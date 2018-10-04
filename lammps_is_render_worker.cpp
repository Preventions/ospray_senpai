// ======================================================================== //
// Copyright 2018 Intel Corporation                                         //
//                                                                          //
// Licensed under the Apache License, Version 2.0 (the "License");          //
// you may not use this file except in compliance with the License.         //
// You may obtain a copy of the License at                                  //
//                                                                          //
//     http://www.apache.org/licenses/LICENSE-2.0                           //
//                                                                          //
// Unless required by applicable law or agreed to in writing, software      //
// distributed under the License is distributed on an "AS IS" BASIS,        //
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. //
// See the License for the specific language governing permissions and      //
// limitations under the License.                                           //
// ======================================================================== //

#include <random>
#include <set>
#include <memory>
#include <algorithm>
#include <array>
#include <chrono>
#include <mpi.h>
#include <unistd.h>
#include <mpiCommon/MPICommon.h>
#include "ospray/ospray_cpp/Camera.h"
#include "ospray/ospray_cpp/Data.h"
#include "ospray/ospray_cpp/Device.h"
#include "ospray/ospray_cpp/FrameBuffer.h"
#include "ospray/ospray_cpp/Geometry.h"
#include "ospray/ospray_cpp/Renderer.h"
#include "ospray/ospray_cpp/TransferFunction.h"
#include "ospray/ospray_cpp/Model.h"
#include "libIS/is_render.h"
#include "colormap.h"
#include "util.h"
#include "image_util.h"
#include "client_server.h"
#include "query_task.h"
#include "kd_tree.h"

using namespace ospcommon;
using namespace ospray::cpp;

float radius = 1.0f;
std::string sim_host;
int sim_port = -1;
int client_port = -1;
bool mpi_multilaunch = false;
bool render_splatter = false;
std::string frame_output_prefix;
std::string benchmark_log_file;
size_t bench_frames = 0;
float bond_dist = 0.f;
int bond_atom_type = -1;
int ao_samples = 1;
bool allow_resize = false;

MPI_Comm worker_comm = MPI_COMM_NULL;
MPI_Comm ospray_comm = MPI_COMM_NULL;

void run_renderer();

box3f make_box(const libISBox3f &box) {
	return box3f(vec3f(box.min.x, box.min.y, box.min.z),
			vec3f(box.max.x, box.max.y, box.max.z));
}

int main(int argc, char **argv) {
	int provided = 0;
	MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
	if (provided != MPI_THREAD_MULTIPLE) {
		std::cerr << "MPI_THREAD_MULTIPLE is required" << std::endl;
		throw std::runtime_error("MPI Thread Multiple Required");
	}

	std::vector<std::string> args(argv, argv + argc);

	for (size_t i = 1; i < args.size(); ++i) {
		if (args[i] == "-sim-host" ) {
			sim_host = argv[++i];
		} else if (args[i] == "-sim-port" ) {
			sim_port = std::stoi(argv[++i]);
		} else if (args[i] == "-port") {
			client_port = std::stoi(argv[++i]);
		} else if (args[i] == "-mpi-multi") {
			mpi_multilaunch = true;
		} else if (args[i] == "-splatter") {
			render_splatter = true;
		} else if (args[i] == "-radius") {
			radius = std::stof(args[++i]);
		} else if (args[i] == "-o") {
			frame_output_prefix = args[++i];
		} else if (args[i] == "-bench") {
			benchmark_log_file = args[++i];
			bench_frames = std::stoull(args[++i]);
		} else if (args[i] == "-bond") {
			bond_atom_type = std::stoi(args[++i]);
			bond_dist = std::stof(args[++i]);
		} else if (args[i] == "-ao") {
			ao_samples = std::stoi(args[++i]);
		} else if (args[i] == "-allow-resize") {
			allow_resize = true;
		}
	}
	if (((sim_host.empty() || sim_port == -1) && !mpi_multilaunch) || client_port == -1) {
		std::cout << "Usage: mpirun -np <N> ./lammps_is_render_worker [options]\n"
			<< "Options:\n"
			<< "  -sim-host <string>     The hostname of the simulation\n"
			<< "  -sim-port <port>       The port libIS on the sim is listening on\n"
			<< "  -port <port>           The port to listen for the viewer client on\n"
			<< "  -mpi-multi             Pass this flag to indicate the run is with MPI's\n"
			<< "                         multi-program launch mode. In this case sim host and\n"
			<< "                         port are not needed\n"
			<< "  -splatter              Render using the distributed splatter instead (requires module_splatter)\n"
			<< "  -no-bonds              Don't compute bonds\n"
			<< "  -radius <float>        Set the radius for atoms\n"
			<< "  -o <file pattern>      Dump frames to files with the prefix 'file_pattern'\n"
			<< "  -bench <log file> <N>  Run in benchmark mode, don't send/recv from a client\n"
			<< "                         and save out perf timings to the log file. In this case\n"
			<< "                         the viewer won't wait for a client to connect and will render N frames\n"
			<< "  -bond <type> <dist>    Compute bonds between atoms of type <type> and distance <dist>\n"
			<< "  -ao <number>           Set the number of AO samples to compute (default 1)\n"
			<< "  -allow-resize          Allow resizing the framebuffer. Note that the number of\n"
			<< "                         rendering clients must evenly divide the number of sim\n"
			<< "                         ranks in this case\n";
		return 1;
	}

	const int ospray_tag = 0x4f535052;
	MPI_Comm_split(MPI_COMM_WORLD, ospray_tag, 0, &worker_comm);
	MPI_Comm_dup(worker_comm, &ospray_comm);

	run_renderer();
	ospShutdown();

	MPI_Comm_free(&worker_comm);
	MPI_Comm_free(&ospray_comm);
	MPI_Finalize();
	return 0;
}
void run_renderer() {
	ospLoadModule("mpi");
	if (render_splatter) {
		ospLoadModule("splatter");
	}
	Device device("mpi_distributed");

	ospDeviceSetVoidPtr(device.handle(), "worldCommunicator",
			static_cast<void*>(&ospray_comm));

	device.set("masterRank", 0);
	//device.set("logLevel", 100);
	ospDeviceSetErrorFunc(device.handle(),
		[](OSPError err, const char *details) {
			std::cout << "OSPRay Error: " << details << std::endl;
		});
	ospDeviceSetStatusFunc(device.handle(),
		[](const char *status) {
			std::cout << "OSPRay Status: " << status << std::endl;
		});
	device.commit();
	device.setCurrent();

	const int rank = mpicommon::world.rank;
	const int world_size = mpicommon::world.size;
	std::cout << "OSPRay with rank " << rank << ", world size: " << world_size << std::endl;

	AppState app;
	AppData appdata;
	Model model;

	if (mpi_multilaunch) {
		std::cout << "Connecting with existing comm for renderer" << std::endl;
		is::render::connectWithExisting(worker_comm, MPI_COMM_WORLD);
	} else {
		std::cout << "Connecting over the network to the simulation" << std::endl;
		is::render::connect(sim_host, sim_port, MPI_COMM_WORLD);
	}

	const bool benchmarking = !benchmark_log_file.empty();
	std::unique_ptr<ClientConnection> client = nullptr;
	char hostname[1024] = {0};
	gethostname(hostname, 1023);
	std::cout << "rank " << rank << " running on " << hostname << "\n";
	if (rank == 0 && !benchmarking) {
		std::cout << "Now listening for client on " << hostname << ":" << client_port << std::endl;
		client = ospcommon::make_unique<ClientConnection>(client_port);
	}

	// Spawn a background thread to query the sim for future timesteps
	auto query_task = std::unique_ptr<QueryTask>(new QueryTask(worker_comm));

	FrameBuffer fb(app.fbSize, OSP_FB_SRGBA, OSP_FB_COLOR | OSP_FB_DEPTH | OSP_FB_ACCUM);
	fb.clear(OSP_FB_COLOR | OSP_FB_DEPTH | OSP_FB_ACCUM);

	Camera camera("perspective");
	camera.set("pos", vec3f(0, 0, -500));
	camera.set("dir", vec3f(0, 0, 1));
	camera.set("up", vec3f(0, 1, 0));
	camera.set("aspect", static_cast<float>(app.fbSize.x) / app.fbSize.y);
	camera.commit();

	Renderer renderer(render_splatter ? "mpi_splatter" : "mpi_raycast");
	renderer.set("camera", camera);
	renderer.set("aoSamples", ao_samples);
	if (render_splatter) {
		renderer.set("bgColor", vec3f(0.0));
	} else {
		renderer.set("bgColor", vec3f(1.0));
	}

	JPGCompressor output_compressor(90, true);

	std::shared_ptr<std::ostream> benchmark_out;
	if (rank == 0 && benchmarking) {
		benchmark_out = std::make_shared<std::ofstream>(benchmark_log_file.c_str());
		*benchmark_out << "Render worker rendering perf on " << world_size << " nodes\n";
	}

	std::shared_ptr<SimulationState> current_state, new_state;
	bool sent_world_bounds = false;
	int frame_number = 0;
	float radians_per_frame = 0.03;
	while (!app.quit && !(bench_frames > 0 && frame_number >= bench_frames)) {
		using namespace std::chrono;

		if (current_state) {
			auto startFrame = high_resolution_clock::now();
			renderer.renderFrame(fb, OSP_FB_COLOR | OSP_FB_DEPTH);
			auto endFrame = high_resolution_clock::now();

			if (rank == 0) {
				const int frameTime = duration_cast<milliseconds>(endFrame - startFrame).count();

				uint32_t *img = (uint32_t*)fb.map(OSP_FB_COLOR);
				float *depth = (float*)fb.map(OSP_FB_DEPTH);

				if (!frame_output_prefix.empty() && frame_number > 0) {
					std::string fname(frame_output_prefix.size() + 4, '_');
					std::snprintf(&fname[0], fname.size() + 1, "%s%04d",
							frame_output_prefix.c_str(), frame_number);
					fname += ".jpg";
					auto jpg = output_compressor.compress(img, app.fbSize.x, app.fbSize.y);
					std::ofstream fout(fname.c_str(), std::ios::binary);
					fout.write((const char*)jpg.first, jpg.second);
				}

				if (!benchmarking) {
					client->send_frame(img, depth, app.fbSize.x, app.fbSize.y, frameTime);
					client->recieve_app_state(app, appdata);
				} else if (frame_number > 0) {
					*benchmark_out << frameTime << "ms" << std::endl;
				}
				fb.unmap(depth);
				fb.unmap(img);
			}
			++frame_number;
		}
		if (rank == 0 && current_state && benchmarking) {
			const vec3f diag = current_state->world_bounds.upper - current_state->world_bounds.lower;
			const float orbit_dist = 1.1 * length(diag);
			const vec3f center = current_state->world_bounds.center();
			vec3f p(orbit_dist * std::sin(frame_number * radians_per_frame),
					0.f,
					orbit_dist * std::cos(frame_number * radians_per_frame));
			app.v[0] = p + center;
			app.v[1] = center - app.v[0];
			app.v[2] = vec3f(0, 1, 0);
			app.cameraChanged = true;
		}

		// Send out the shared app state that the workers need to know, e.g. camera
		// position, if we should be quitting.
		MPI_Bcast(&app, sizeof(AppState), MPI_BYTE, 0, worker_comm);

		if (app.cameraChanged) {
			camera.set("pos", app.v[0]);
			camera.set("dir", app.v[1]);
			camera.set("up", app.v[2]);
			camera.commit();

			fb.clear(OSP_FB_COLOR | OSP_FB_DEPTH | OSP_FB_ACCUM);
			app.cameraChanged = false;
		}
		
		// TODO: This is buggy because if we have different numbers of
		// regiosn on the nodes the new
		// the framebuffer object ID will be different across nodes and we won't
		// be able to send tiles anymore.
		if (allow_resize && app.fbSizeChanged) {
			fb = FrameBuffer(app.fbSize, OSP_FB_SRGBA, OSP_FB_COLOR | OSP_FB_DEPTH | OSP_FB_ACCUM);
			fb.clear(OSP_FB_COLOR | OSP_FB_DEPTH | OSP_FB_ACCUM);
			camera.set("aspect", static_cast<float>(app.fbSize.x) / app.fbSize.y);
			camera.commit();
			app.fbSizeChanged = false;
		}

		if (!new_state) {
			new_state = query_task->take();
		}
		int got_new_state = new_state ? 1 : 0;
		int any_not_updated = 0;
		// Everyone has to get their next set of regions before we can switch
		MPI_Allreduce(&got_new_state, &any_not_updated, 1, MPI_INT, MPI_MIN, worker_comm);
		if (any_not_updated != 0) {
			current_state = new_state;
#if 1
			std::vector<DistributedRegion> region_bounds;
			std::vector<box3f> ghost_regions;
			for (const auto &r : current_state->regions) {
				region_bounds.emplace_back(make_box(r.local), r.simRank);

				box3f ghost_box = make_box(r.ghost);
				if (!ghost_box.empty()) {
					ghost_regions.push_back(ghost_box);
				} else {
					// TODO: Get the actual ghost bounds, seems like there's some bug with
					// using the region bounds?
					//ghost_regions.push_back(region_bounds.back().bounds);
					ghost_regions.push_back(current_state->world_bounds);
				}
			}
			Data region_data(region_bounds.size() * sizeof(DistributedRegion),
					OSP_CHAR, region_bounds.data());

			current_state->model.set("regions", region_data);

			if (!ghost_regions.empty()) {
				Data ghost_region_data(ghost_regions.size() * 2,
						OSP_FLOAT3, ghost_regions.data());
				current_state->model.set("ghostRegions", ghost_region_data);
			}
#endif
			current_state->model.commit();
			renderer.set("model", current_state->model);
			renderer.commit();
			fb.clear(OSP_FB_COLOR | OSP_FB_DEPTH | OSP_FB_ACCUM);
			MPI_Barrier(worker_comm);
			new_state = nullptr;
		}

		if (rank == 0 && current_state && !sent_world_bounds) {
			if (!benchmarking) {
				client->send_metadata(current_state->world_bounds);
			}
			sent_world_bounds = true;
		}
	}

	benchmark_out = nullptr;
	query_task = nullptr;
	is::render::disconnect();
}

