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
#include <mutex>
#include <memory>
#include <thread>
#include <atomic>
#include <cstring>
#include <chrono>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <thread>

#include <glm/glm.hpp>
#include <glm/ext.hpp>
#include <GL/gl3w.h>
#include <GLFW/glfw3.h>
#include <mpi.h>
#include <mpiCommon/MPICommon.h>
#include <ospray_cpp/Camera.h>
#include <ospray_cpp/Data.h>
#include <ospray_cpp/Device.h>
#include <ospray_cpp/FrameBuffer.h>
#include <ospray_cpp/Geometry.h>
#include <ospray_cpp/Renderer.h>
#include <ospray_cpp/TransferFunction.h>
#include <ospray_cpp/Volume.h>
#include <ospray_cpp/Model.h>
#include "imgui.h"
#include "imgui_impl_glfw_gl3.h"

#include "libIS/is_render.h"
#include "arcball_camera.h"
#include "save_bmp.h"
#include "glshaders.h"
#include "gldebug.h"
#include "colormap.h"
#include "app_util.h"
#include "lammps_driver.h"

using namespace ospray::cpp;
using namespace ospcommon;

bool show_ui = true;
bool take_screenshot = false;
int rank, world_size;
float sphere_radius = 1;
bool query_in_background = false;
std::shared_ptr<LammpsDriver> lammps_driver;

// Commandline params
vec3f cam_pos(0, 0, 2);
vec3f cam_target(0, 0, 0);
vec3f cam_up(0, 1, 0);
size_t NUM_SPHERES = 50;
size_t lammps_timesteps = 0;
std::string sim_server;
int sim_port = -1;
bool open_viewer_window = true;
int wait_period = 0;
std::string output_prefix;
std::string lammps_file;

struct Sphere {
	vec3f org;
	int colorID{0};
};

struct LammpsAtom {
	double x, y, z;
};

struct AOTweakParams {
	int n_samples{27};
	int n_turns{16};
	float ball_radius{1.0};
	float sigma{1.2};
	float kappa{1.1};
	int filter_scale{1};
	float edge_sharpness{10};
	int num_blur_passes{2};

	void draw_ui() {
		if (ImGui::Begin("AO Params")) {
			ImGui::SliderInt("n_samples", &n_samples, 1, 100);
			ImGui::SliderInt("n_turns", &n_turns, 1, 100);
			ImGui::SliderFloat("ball_radius", &ball_radius, 0.4, 8);
			ImGui::SliderFloat("sigma", &sigma, 0.1, 5.0);
			ImGui::SliderFloat("kappa", &kappa, 0.1, 5.0);
			ImGui::SliderInt("filter_scale", &filter_scale, 1, 6);
			ImGui::SliderFloat("edge_sharpness", &edge_sharpness, 1, 100);
			ImGui::InputInt("blur passes", &num_blur_passes);
			num_blur_passes = std::max(0, num_blur_passes);
		}
		ImGui::End();
	}
};

void key_callback(GLFWwindow *window, int key, int sc, int action, int mods) {
	WindowState *state = static_cast<WindowState*>(glfwGetWindowUserPointer(window));
	if (action == GLFW_PRESS) {
		switch (key) {
			case GLFW_KEY_ESCAPE:
				glfwSetWindowShouldClose(window, true);
				break;
			case GLFW_KEY_R:
				state->camera.reset();
				break;
			case GLFW_KEY_H:
				show_ui = !show_ui;
				break;
			case GLFW_KEY_P:
				take_screenshot = true;
				break;
			default:
				break;
		}
	}
	ImGui_ImplGlfwGL3_KeyCallback(window, key, sc, action, mods);
}
void cursor_pos_callback(GLFWwindow *window, double x, double y) {
	WindowState *state = static_cast<WindowState*>(glfwGetWindowUserPointer(window));
	const glm::vec2 mouse(x, y);
	if (state->prev_mouse != glm::vec2(-1)) {
		const bool left_down = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS;
		const bool right_down = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT) == GLFW_PRESS;
		const bool middle_down = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_MIDDLE) == GLFW_PRESS;
		if (!ImGui::IsMouseHoveringAnyWindow()) {
			if (left_down) {
				state->camera.rotate(mouse, mouse - state->prev_mouse, state->frame_delta);
				state->camera_changed = true;
			} else if (right_down) {
				state->camera.zoom(mouse.y - state->prev_mouse.y, state->frame_delta);
				state->camera_changed = true;
			} else if (middle_down) {
				state->camera.pan(mouse - state->prev_mouse, state->frame_delta);
				state->camera_changed = true;
			}
		}
	}
	state->prev_mouse = mouse;
}
void framebuffer_size_callback(GLFWwindow *window, int width, int height) {
	WindowState *state = static_cast<WindowState*>(glfwGetWindowUserPointer(window));
	state->app.fb_size = vec2i(width, height);
	state->app.fb_size_changed = true;
}

// Compute new colormap for the updated atom type attributes. Also re-scales so atom type min = 0
Data update_atom_types(std::vector<is::render::Region> &regions) {
	std::array<int, 2> atomTypeRange = {
		std::numeric_limits<int>::max(),
		std::numeric_limits<int>::min()
	};

	for (const auto &region : regions) {
		const Sphere *atoms = reinterpret_cast<const Sphere*>(region.particles.data());
		// TODO: We need to do a reduction to compute the global min/mix ids and generate
		// colors for the unique atom types
		auto min_max_type = std::minmax_element(atoms, atoms + region.particles.size() / sizeof(Sphere),
			[](const Sphere &a, const Sphere &b) {
				return a.colorID < b.colorID;
			});
		const int localMin = std::min(min_max_type.first->colorID, atomTypeRange[0]);
		const int localMax = std::max(min_max_type.second->colorID, atomTypeRange[1]);
		MPI_Allreduce(&localMin, &atomTypeRange[0], 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
		MPI_Allreduce(&localMax, &atomTypeRange[1], 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
	}
	if (atomTypeRange[0] != 0) {
		for (auto &region : regions) {
			Sphere *atoms = reinterpret_cast<Sphere*>(region.particles.data());
			const size_t n_atoms = region.particles.size() / sizeof(Sphere);
			for (size_t i = 0; i < n_atoms; ++i) {
				atoms[i].colorID -= atomTypeRange[0];
			}
		}
		atomTypeRange[1] -= atomTypeRange[0];
		atomTypeRange[0] = 0;
	}
	std::cout << "Got atom type range {" << atomTypeRange[0] << ", " << atomTypeRange[1] << "}\n";
	std::vector<vec4f> atomColors;
	Colormap colormap;
	for (int i = atomTypeRange[0]; i <= atomTypeRange[1]; ++i) {
		const float x = (i - atomTypeRange[0]) / static_cast<float>(atomTypeRange[1] - atomTypeRange[0]);
		atomColors.push_back(vec4f(colormap.map(x), 1));
	}

	Data color_data = Data(atomColors.size(), OSP_FLOAT4, atomColors.data());
	color_data.commit();
	return color_data;
}
void update_geometry(std::vector<is::render::Region> &regions,
		std::vector<Geometry> &geometries, Model &model)
{
	Data color_data = update_atom_types(regions);
	for (size_t i = 0; i < regions.size(); ++i) {
		is::render::Region &region = regions[i];
		std::cout << "Queried region with " << region.numParticles << " particles\n"
			<< "bounds = [(" << region.bounds.min.x << ", " << region.bounds.min.y
			<< ", " << region.bounds.min.z << "), (" << region.bounds.max.x
			<< ", " << region.bounds.max.y << ", " << region.bounds.max.z << ")]\n";
		if (i >= geometries.size()) {
			geometries.emplace_back("spheres");
			model.addGeometry(geometries[i]);
		}
		Data sphere_data = Data(region.particles.size(), OSP_UCHAR, region.particles.data());
		sphere_data.commit();
		Geometry &geom = geometries[i];
		geom.set("spheres", sphere_data);
		geom.set("color", color_data);
		geom.set("offset_colorID", int(sizeof(vec3f)));
		geom.set("radius", sphere_radius);
		geom.commit();
	}
}

void run_ospray_rank(MPI_Comm partition_comm, MPI_Comm insitu_comm);
void run_lammps_rank();

int main(int argc, char **argv) {
	int provided = 0;
	MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
	query_in_background = provided == MPI_THREAD_MULTIPLE;
	std::vector<std::string> args(argv, argv + argc);

	for (size_t i = 0; i < args.size(); ++i) {
		if (args[i] == "-sim") {
			sim_server = args[++i];
		} else if (args[i] == "-port") {
			sim_port = std::stoi(args[++i]);
		} else if (args[i] == "-n") {
			NUM_SPHERES = std::stol(args[++i]);
		} else if (args[i] == "-vp") {
			for (size_t j = 0; j < 3; ++j) {
				cam_pos[j] = std::stof(args[++i]);
			}
		} else if (args[i] == "-vi") {
			for (size_t j = 0; j < 3; ++j) {
				cam_target[j] = std::stof(args[++i]);
			}
		} else if (args[i] == "-vu") {
			for (size_t j = 0; j < 3; ++j) {
				cam_up[j] = std::stof(args[++i]);
			}
		} else if (args[i] == "-o") {
			output_prefix = args[++i];
		} else if (args[i] == "-no-viewer") {
			open_viewer_window = false;
		} else if (args[i] == "-radius") {
			sphere_radius = std::stof(args[++i]);
		} else if (args[i] == "-wait") {
			wait_period = std::stoi(args[++i]);
		} else if (args[i] == "-lammps") {
			lammps_file = args[++i];
		} else if (args[i] == "-lammps-steps") {
			lammps_timesteps = std::stoi(args[++i]);
		}
	}
	if (!lammps_file.empty() && provided != MPI_THREAD_MULTIPLE) {
		throw std::runtime_error("Thread multiple is required to run LAMMPS in app!");
	}
	if (wait_period > 0) {
		std::cout << "waiting for " << wait_period << "s\n";
		std::this_thread::sleep_for(std::chrono::seconds(wait_period));
	}

	MPI_Comm insitu_comm;
	MPI_Comm_dup(MPI_COMM_WORLD, &insitu_comm);
	if (!lammps_file.empty()) {
		lammps_driver = std::make_shared<LammpsDriver>(lammps_file, insitu_comm);
	}

	run_ospray_rank(MPI_COMM_WORLD, insitu_comm);

	if (lammps_driver) {
		is::render::disconnect();
		lammps_driver = nullptr;
	} else if (sim_port > 0) {
		is::render::disconnect();
	}
	MPI_Comm_free(&insitu_comm);

	MPI_Finalize();

	return 0;
}
void run_ospray_rank(MPI_Comm partition_comm, MPI_Comm insitu_comm) {
	ospLoadModule("mpi");

	Device device("mpi_distributed");
	device.set("masterRank", 0);
	device.commit();
	device.setCurrent();

	rank = mpicommon::world.rank;
	world_size = mpicommon::world.size;

	AppState app;
	Model model;
	std::vector<Geometry> geometries;
	std::unique_ptr<QueryTask> query_task;
	if (sim_port > 0 || !lammps_file.empty()) {
		if (lammps_file.empty()) {
			is::render::connect(sim_server, sim_port, insitu_comm);
		} else {
			is::render::connectWithExisting(insitu_comm, insitu_comm);
		}
		std::vector<is::render::Region> regions = is::render::query();

		update_geometry(regions, geometries, model);

		query_task = std::unique_ptr<QueryTask>(new QueryTask(query_in_background));
	} else {
		std::mt19937 rng;
		rng.seed(std::random_device()());

		std::uniform_real_distribution<float> dist(-10, 10);
		std::uniform_int_distribution<int> colorIdDistrib(0, 100);
		std::vector<Sphere> spheres(NUM_SPHERES);
		for (auto &s : spheres) {
			s.org = vec3f(dist(rng), dist(rng), dist(rng));
			s.colorID = colorIdDistrib(rng);
		}
		Data sphere_data = Data(NUM_SPHERES * sizeof(Sphere), OSP_UCHAR, spheres.data());
		sphere_data.commit();

		Colormap colormap;
		std::vector<vec4f> atomColors;
		for (int i = colorIdDistrib.min(); i < colorIdDistrib.max(); ++i) {
			const float x = (i - colorIdDistrib.min()) / static_cast<float>(colorIdDistrib.max() - colorIdDistrib.min());
			atomColors.push_back(vec4f(colormap.map(x), 1));
		}

		Data color_data = Data(atomColors.size(), OSP_FLOAT4, atomColors.data());
		color_data.commit();

		Geometry geom("spheres");
		geom.set("spheres", sphere_data);
		geom.set("color", color_data);
		geom.set("offset_colorID", int(sizeof(vec3f)));
		geom.set("radius", sphere_radius);
		geom.commit();
		model.addGeometry(geom);
	}
	model.commit();

	Camera camera("perspective");
	camera.set("pos", cam_pos);
	camera.set("dir", cam_target - cam_pos);
	camera.set("up", cam_up);
	camera.set("aspect", static_cast<float>(app.fb_size.x) / app.fb_size.y);
	camera.commit();

	Renderer renderer("mpi_raycast");
	// Should just do 1 set here, which is read?
	renderer.set("world", model);
	renderer.set("model", model);
	renderer.set("camera", camera);
	renderer.set("bgColor", vec3f(1.0));
	renderer.commit();
	assert(renderer);

	FrameBuffer fb(app.fb_size, OSP_FB_SRGBA, OSP_FB_COLOR | OSP_FB_DEPTH | OSP_FB_ACCUM);
	fb.clear(OSP_FB_COLOR | OSP_FB_DEPTH | OSP_FB_ACCUM);

	mpicommon::world.barrier();

	std::shared_ptr<WindowState> window_state =
		std::make_shared<WindowState>(app, cam_pos, cam_target, cam_up);
	window_state->frame_delta = 0.16;

	AOTweakParams ao_params;
	GLuint vao, ao_params_buffer;
	// Color, depth textures for OSPRay rendered results
	std::array<GLuint, 2> osp_textures;
	// For blur pass intermediate targets
	std::array<GLuint, 2> blur_textures, blur_fbos;

	std::unique_ptr<GLShader> ao_shader, blur_shader, composite_shader;

	GLFWwindow *window = nullptr;
	if (rank == 0 && open_viewer_window) {
		if (!glfwInit()) {
			std::exit(1);
		}
		glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
		glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
		glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
		//glfwWindowHint(GLFW_OPENGL_DEBUG_CONTEXT, GL_TRUE);
		window = glfwCreateWindow(app.fb_size.x, app.fb_size.y,
				"LAMMPS + SENSEI Viewer", nullptr, nullptr);
		if (!window) {
			glfwTerminate();
			std::exit(1);
		}
		glfwMakeContextCurrent(window);

		if (gl3wInit()) {
			std::cout << "Failed to init gl3w\n";
			glfwTerminate();
			std::exit(1);
		}
		glDisable(GL_DEPTH_TEST);

		// Install callbacks for the ones I don't need to handle them,
		// in the ones I do use I'll forward the call on
		ImGui_ImplGlfwGL3_Init(window, true);

		glfwSetKeyCallback(window, key_callback);
		glfwSetCursorPosCallback(window, cursor_pos_callback);
		glfwSetWindowUserPointer(window, window_state.get());
		glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);

		glGenVertexArrays(1, &vao);
		glBindVertexArray(vao);
		ao_shader = ospcommon::make_unique<GLShader>(vsrc, fAOComputeSrc);
		blur_shader = ospcommon::make_unique<GLShader>(vsrc, fBlurSrc);
		composite_shader = ospcommon::make_unique<GLShader>(vsrc, fFinalCompositeSrc);

		glUseProgram(ao_shader->shader);
		glUniform2f(ao_shader->uniform_location("viewport_dim"), app.fb_size.x, app.fb_size.y);
		glUniform1i(ao_shader->uniform_location("camera_depths"), 1);
		// Mark the camera changed to init the uniforms
		window_state->camera_changed = true;

		glGenBuffers(1, &ao_params_buffer);
		glBindBuffer(GL_UNIFORM_BUFFER, ao_params_buffer);
		glBufferData(GL_UNIFORM_BUFFER, sizeof(AOTweakParams), &ao_params, GL_STREAM_DRAW);
		glBindBufferBase(GL_UNIFORM_BUFFER, 0, ao_params_buffer);

		glUniformBlockBinding(ao_shader->shader,
				glGetUniformBlockIndex(ao_shader->shader, "AOParams"), 0);
		glUniformBlockBinding(blur_shader->shader,
				glGetUniformBlockIndex(blur_shader->shader, "AOParams"), 0);

		glGenTextures(osp_textures.size(), osp_textures.data());
		glBindTexture(GL_TEXTURE_2D, osp_textures[0]);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, app.fb_size.x, app.fb_size.y,
				0, GL_RGBA, GL_UNSIGNED_BYTE, nullptr);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);

		glActiveTexture(GL_TEXTURE1);
		glBindTexture(GL_TEXTURE_2D, osp_textures[1]);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, app.fb_size.x, app.fb_size.y,
				0, GL_RED, GL_FLOAT, nullptr);

		glGenTextures(blur_textures.size(), blur_textures.data());
		glActiveTexture(GL_TEXTURE2);
		for (auto &t : blur_textures) {
			glBindTexture(GL_TEXTURE_2D, t);
			glTexImage2D(GL_TEXTURE_2D, 0, GL_RG32F, app.fb_size.x, app.fb_size.y,
					0, GL_RG, GL_FLOAT, nullptr);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);
		}

		glGenFramebuffers(blur_fbos.size(), blur_fbos.data());
		for (size_t i = 0; i < blur_fbos.size(); ++i) {
			glBindFramebuffer(GL_FRAMEBUFFER, blur_fbos[i]);
			glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0,
					GL_TEXTURE_2D, blur_textures[i], 0);
			assert(glCheckFramebufferStatus(GL_FRAMEBUFFER) == GL_FRAMEBUFFER_COMPLETE);
		}
		glBindFramebuffer(GL_FRAMEBUFFER, 0);
	}

	float ospray_frame_time_avg = 0;
	size_t frame_counter = 0;
	size_t timestep = 0;
	bool new_timestep = true;
	std::vector<is::render::Region> regions;
	while (!app.quit) {
		if (app.camera_changed) {
			camera.set("pos", app.v[0]);
			camera.set("dir", app.v[1]);
			camera.set("up", app.v[2]);
			camera.commit();
			fb.clear(OSP_FB_COLOR | OSP_FB_DEPTH | OSP_FB_ACCUM);
			app.camera_changed = false;
			if (rank == 0) {
				window_state->camera_changed = false;
				glm::vec3 cam_du, cam_dv;
				app.compute_cam_dirs(cam_du, cam_dv);

				if (open_viewer_window) {
					glUseProgram(ao_shader->shader);
					glUniform3fv(ao_shader->uniform_location("cam_du"), 1, glm::value_ptr(cam_du));
					glUniform3fv(ao_shader->uniform_location("cam_dv"), 1, glm::value_ptr(cam_dv));
					glUniform3fv(ao_shader->uniform_location("cam_pos"), 1, &app.v[0].x);
					glUniform3fv(ao_shader->uniform_location("cam_dir"), 1, &app.v[1].x);
				}
			}
		}

		using namespace std::chrono;
		auto start_frame = high_resolution_clock::now();
		renderer.renderFrame(fb, OSP_FB_COLOR | OSP_FB_DEPTH);
		auto end_frame = high_resolution_clock::now();
		const int osprayRenderTime = duration_cast<milliseconds>(end_frame - start_frame).count();
		if (!open_viewer_window && rank == 0) {
			std::cout << "OSPRay Render Frame: " << osprayRenderTime << "ms\n";
		}
		ospray_frame_time_avg = (osprayRenderTime + frame_counter * ospray_frame_time_avg) / (frame_counter + 1);
		++frame_counter;

		if (rank == 0) {
			uint32_t *img = (uint32_t*)fb.map(OSP_FB_COLOR);
			float *imgDepth = (float*)fb.map(OSP_FB_DEPTH);
			if (open_viewer_window) {
				glActiveTexture(GL_TEXTURE0);
				glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, app.fb_size.x, app.fb_size.y,
						0, GL_RGBA, GL_UNSIGNED_BYTE, img);
				glActiveTexture(GL_TEXTURE1);
				glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, app.fb_size.x, app.fb_size.y,
						0, GL_RED, GL_FLOAT, imgDepth);
				glGenerateMipmap(GL_TEXTURE_2D);
				glClear(GL_COLOR_BUFFER_BIT);

				auto start_ao = high_resolution_clock::now();
				// AO compute pass
				glBindFramebuffer(GL_FRAMEBUFFER, blur_fbos[0]);
				glUseProgram(ao_shader->shader);
				glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);

				// Blur passes
				glActiveTexture(GL_TEXTURE2);
				glUseProgram(blur_shader->shader);
				glUniform1i(blur_shader->uniform_location("ao_in"), 2);

				for (size_t i = 0; i < ao_params.num_blur_passes; ++i) {
					// Horizontal blur pass
					glBindFramebuffer(GL_FRAMEBUFFER, blur_fbos[1]);
					glBindTexture(GL_TEXTURE_2D, blur_textures[0]);
					glUniform2i(blur_shader->uniform_location("axis"), 1, 0);
					glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);

					// Vertical blur pass
					glBindFramebuffer(GL_FRAMEBUFFER, blur_fbos[0]);
					glBindTexture(GL_TEXTURE_2D, blur_textures[1]);
					glUniform2i(blur_shader->uniform_location("axis"), 0, 1);
					glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
				}

				// Final AO composite pass
				glBindTexture(GL_TEXTURE_2D, blur_textures[0]);
				glBindFramebuffer(GL_FRAMEBUFFER, 0);
				glUseProgram(composite_shader->shader);
				glUniform1i(composite_shader->uniform_location("camera_colors"), 0);
				glUniform1i(composite_shader->uniform_location("ao_vals"), 2);
				glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
				glFinish();
				auto end_ao = high_resolution_clock::now();
				const int aoComputeTime = duration_cast<milliseconds>(end_ao - start_ao).count();

				// Draw ImGui on top
				ImGui_ImplGlfwGL3_NewFrame();
				if (show_ui) {
					ao_params.draw_ui();

					ImGui::Text("OSPRay Render Time: %dms\nAO Computation Time: %dms",
							osprayRenderTime, aoComputeTime);
				}
				ImGui::Render();

				glBufferData(GL_UNIFORM_BUFFER, sizeof(AOTweakParams), &ao_params, GL_STREAM_DRAW);

				glfwSwapBuffers(window);

				if (take_screenshot) {
					std::vector<uint8_t> img(app.fb_size.x * app.fb_size.y * 4, 0);
					glReadPixels(0, 0, app.fb_size.x, app.fb_size.y, GL_BGRA, GL_UNSIGNED_BYTE, img.data());

					std::stringstream fname;
					fname << output_prefix << "_" << std::setfill('0') << std::setw(4) << timestep << ".bmp";
					if (!save_bmp(fname.str(), app.fb_size.x, app.fb_size.y, img.data())) {
						throw std::runtime_error("failed to save image");
					}
					std::cout << "Saved screenshot to " << timestep << " to '" << fname.str() << "'\n";
					take_screenshot = false;
				}

				glfwPollEvents();
				if (glfwWindowShouldClose(window)) {
					app.quit = true;
				}
			} else if (new_timestep) {
				std::stringstream fname;
				fname << output_prefix << "_" << std::setfill('0') << std::setw(4) << timestep << ".bmp";
				// Go through and swap RGBA format to BGRA
				for (size_t i = 0; i < app.fb_size.x * app.fb_size.y; ++i) {
					uint8_t *px = reinterpret_cast<uint8_t*>(&img[i]);
					std::swap(px[0], px[2]);
				}
				if (!save_bmp(fname.str(), app.fb_size.x, app.fb_size.y, reinterpret_cast<uint8_t*>(img))) {
					throw std::runtime_error("failed to save image");
				}
				std::cout << "Saved timestep " << timestep << " to '" << fname.str() << "'\n";
				new_timestep = false;
			}
			fb.unmap(img);
			fb.unmap(imgDepth);

			const glm::vec3 eye = window_state->camera.eye_pos();
			const glm::vec3 look = window_state->camera.look_dir();
			const glm::vec3 up = window_state->camera.up_dir();
			app.v[0] = vec3f(eye.x, eye.y, eye.z);
			app.v[1] = vec3f(look.x, look.y, look.z);
			app.v[2] = vec3f(up.x, up.y, up.z);
			app.camera_changed = window_state->camera_changed;
		}

		if (lammps_driver && lammps_timesteps > 0 && timestep >= lammps_timesteps) {
			app.quit = true;
		}
		// Send out the shared app state that the workers need to know, e.g. camera
		// position, if we should be quitting.
		MPI_Bcast(&app, sizeof(AppState), MPI_BYTE, 0, MPI_COMM_WORLD);

		if (app.fb_size_changed) {
			fb = FrameBuffer(app.fb_size, OSP_FB_SRGBA, OSP_FB_COLOR | OSP_FB_DEPTH | OSP_FB_ACCUM);
			fb.clear(OSP_FB_COLOR | OSP_FB_DEPTH | OSP_FB_ACCUM);
			camera.set("aspect", static_cast<float>(app.fb_size.x) / app.fb_size.y);
			camera.commit();
			app.fb_size_changed = false;

			if (rank == 0) {
				window_state->camera.update_screen(app.fb_size.x, app.fb_size.y);
				if (open_viewer_window) {
					glViewport(0, 0, app.fb_size.x, app.fb_size.y);
					glUseProgram(ao_shader->shader);
					glUniform2f(ao_shader->uniform_location("viewport_dim"), app.fb_size.x, app.fb_size.y);

					// Need to re-size our intermediate targets
					glActiveTexture(GL_TEXTURE2);
					for (auto &t : blur_textures) {
						glBindTexture(GL_TEXTURE_2D, t);
						glTexImage2D(GL_TEXTURE_2D, 0, GL_RG32F, app.fb_size.x, app.fb_size.y,
								0, GL_RG, GL_FLOAT, nullptr);
					}
				}
			}
		}

		if (sim_port > 0 || lammps_driver) {
			if (!query_in_background) {
				query_task->queryThread();
			}
			if (regions.empty()) {
				regions = query_task->take();
			}
			const int num_regions = regions.size();
			int min_num_regions = 0;
			// Everyone has to get their next set of regions before we can switch
			MPI_Allreduce(&num_regions, &min_num_regions, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
			if (min_num_regions != 0) {
				++timestep;
				new_timestep = true;

				update_geometry(regions, geometries, model);
				model.commit();
				regions.clear();

				fb.clear(OSP_FB_COLOR | OSP_FB_DEPTH | OSP_FB_ACCUM);
			}
		}
	}
	if (rank == 0 && open_viewer_window) {
		ImGui_ImplGlfwGL3_Shutdown();
		glDeleteVertexArrays(1, &vao);
		glDeleteBuffers(1, &ao_params_buffer);
		glDeleteTextures(osp_textures.size(), osp_textures.data());
		glDeleteTextures(blur_textures.size(), blur_textures.data());
		glDeleteFramebuffers(blur_fbos.size(), blur_fbos.data());
		ao_shader = nullptr;
		blur_shader = nullptr;
		composite_shader = nullptr;
	}
}
void run_lammps_rank() {
}

