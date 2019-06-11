#include <iostream>
#include <memory>
#include <vector>
#include <array>
#include <mpi.h>
#include <ospray.h>
#include <SDL.h>
#include "libIS/is_client.h"
#include "gl_core_4_5.h"
#include "imgui.h"
#include "imgui_impl_sdl.h"
#include "imgui_impl_opengl3.h"
#include "arcball_camera.h"
#include "shader.h"
#include "util.h"
#include "transfer_function_widget.h"

const std::string fullscreen_quad_vs = R"(
#version 330 core

const vec4 pos[4] = vec4[4](
	vec4(-1, 1, 0.5, 1),
	vec4(-1, -1, 0.5, 1),
	vec4(1, 1, 0.5, 1),
	vec4(1, -1, 0.5, 1)
);

void main(void){
	gl_Position = pos[gl_VertexID];
}
)";

const std::string display_texture_fs = R"(
#version 330 core

uniform sampler2D img;

out vec4 color;

void main(void){ 
	ivec2 uv = ivec2(gl_FragCoord.xy);
	color = texelFetch(img, uv, 0);
})";

int win_width = 1280;
int win_height = 720;
int rank = -1;
int world_size = -1;
std::string display_field;
glm::vec2 value_range(0.f, 1.f);

struct AppState {
	glm::vec3 cam_pos, cam_dir, cam_up;
	glm::ivec2 fb_dims;
	bool camera_changed = false;
	bool fbsize_changed = false;
	bool done = false;
	bool tfcn_changed = false;
};

void run_viewer(const std::vector<std::string> &args);
void run_worker(const std::vector<std::string> &args);
OSPModel build_regions(const std::vector<is::SimState> &regions, OSPTransferFunction tfcn);
void log_regions(const std::vector<is::SimState> &regions);
void update_transfer_fcn(OSPTransferFunction tfcn, const std::vector<uint8_t> &colormap);

std::ostream& operator<<(std::ostream &os, const libISBox3f &b) {
	os << "{(" << b.min.x << ", " << b.min.y << ", " << b.min.z
		<< "), (" << b.max.x << ", " << b.max.y << ", " << b.max.z << ")}";
	return os;
}

glm::vec2 transform_mouse(glm::vec2 in) {
	return glm::vec2(in.x * 2.f / win_width - 1.f, 1.f - 2.f * in.y / win_height);
}

int main(int argc, const char **argv) {
	if (argc < 3) {
		std::cout << "Usage: " << argv[0] << " (simhost) (simport)\n";
		return 1;
	}

	int provided = 0;
	MPI_Init_thread(&argc, (char***)&argv, MPI_THREAD_MULTIPLE, &provided);
	if (provided != MPI_THREAD_MULTIPLE && provided != MPI_THREAD_SERIALIZED) {
		std::cout << "MPI Thread Multiple or Serialized is required!\n";
		return 1;
	}
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	ospLoadModule("ispc");
	if (ospLoadModule("mpi") != OSP_NO_ERROR) {
		std::cout << "Failed to load MPI module\n";
		return 1;
	}

	OSPDevice device = ospNewDevice("mpi_distributed");
	ospDeviceCommit(device);
	ospSetCurrentDevice(device);

	std::cout << "Connecting to sim at " << argv[1] << ":" << argv[2] << "\n";
	is::client::connect(argv[1], std::stoi(argv[2]), MPI_COMM_WORLD);

	std::vector<std::string> args(argv, argv + argc);

	for (size_t i = 1; i < args.size(); ++i) {
		if (args[i] == "-field") {
			display_field = args[++i];
		} else if (args[i] == "-range") {
			value_range.x = std::stof(args[++i]);
			value_range.y = std::stof(args[++i]);
		}
	}

	if (rank == 0) {
		run_viewer(args);
	} else {
		run_worker(args);
	}
	ospShutdown();

	is::client::disconnect();
	MPI_Finalize();

	return 0;
}
void run_viewer(const std::vector<std::string> &args) {
	if (SDL_Init(SDL_INIT_EVERYTHING) != 0) {
		std::cerr << "Failed to init SDL: " << SDL_GetError() << "\n";
		throw std::runtime_error("Failed to init SDL");
	}
	const char* glsl_version = "#version 330 core";
	SDL_GL_SetAttribute(SDL_GL_CONTEXT_FLAGS, SDL_GL_CONTEXT_FORWARD_COMPATIBLE_FLAG);
	SDL_GL_SetAttribute(SDL_GL_CONTEXT_PROFILE_MASK, SDL_GL_CONTEXT_PROFILE_CORE);
	SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 3);
	SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 3);

	// Create window with graphics context
	SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);
	SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE, 24);
	SDL_GL_SetAttribute(SDL_GL_STENCIL_SIZE, 8);

	SDL_Window* window = SDL_CreateWindow("OSPRay Senpai",
			SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, win_width, win_height,
			SDL_WINDOW_OPENGL | SDL_WINDOW_RESIZABLE);

	SDL_GLContext gl_context = SDL_GL_CreateContext(window);
	SDL_GL_SetSwapInterval(1);
	SDL_GL_MakeCurrent(window, gl_context);

	if (ogl_LoadFunctions() == ogl_LOAD_FAILED) {
		std::cerr << "Failed to initialize OpenGL\n";
		throw std::runtime_error("Failed to init OpenGL");
	}

	// Setup Dear ImGui context
	ImGui::CreateContext();

	// Setup Dear ImGui style
	ImGui::StyleColorsDark();

	// Setup Platform/Renderer bindings
	ImGui_ImplSDL2_InitForOpenGL(window, gl_context);
	ImGui_ImplOpenGL3_Init(glsl_version);

	ImGuiIO& io = ImGui::GetIO();

	glm::vec3 eye(0, 0, 5);
	glm::vec3 center(0);
	glm::vec3 up(0, 1, 0);
	for (size_t i = 1; i < args.size(); ++i) {
		if (args[i] == "-eye") {
			eye.x = std::stof(args[++i]);
			eye.y = std::stof(args[++i]);
			eye.z = std::stof(args[++i]);
		} else if (args[i] == "-center") {
			center.x = std::stof(args[++i]);
			center.y = std::stof(args[++i]);
			center.z = std::stof(args[++i]);
		} else if (args[i] == "-up") {
			up.x = std::stof(args[++i]);
			up.y = std::stof(args[++i]);
			up.z = std::stof(args[++i]);
		}
	}

	ArcballCamera camera(eye, center, up);

	auto display_render = std::make_unique<Shader>(fullscreen_quad_vs, display_texture_fs);

	GLuint render_texture;
	glGenTextures(1, &render_texture);
	glBindTexture(GL_TEXTURE_2D, render_texture);
	glTexStorage2D(GL_TEXTURE_2D, 1, GL_RGBA8, win_width, win_height);

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

	GLuint vao;
	glCreateVertexArrays(1, &vao);
	glBindVertexArray(vao);

	glClearColor(0.0, 0.0, 0.0, 0.0);
	glDisable(GL_DEPTH_TEST);

	OSPTransferFunction tfcn = ospNewTransferFunction("piecewise_linear");

	OSPFrameBuffer fb = ospNewFrameBuffer(osp::vec2i{win_width, win_height},
			OSP_FB_SRGBA, OSP_FB_COLOR | OSP_FB_ACCUM);

	// Query the first timestep from the simulation
	auto regions = is::client::query();
	log_regions(regions);
	OSPModel world = build_regions(regions, tfcn);

	// Start querying for the next timestep asynchronously
	auto region_future = is::client::query_async();

	OSPCamera ospcamera = ospNewCamera("perspective");
	OSPRenderer renderer = ospNewRenderer("mpi_raycast");
	ospSetObject(renderer, "model", world);
	ospSetObject(renderer, "camera", ospcamera);
	// Note: we commit for the first time in the render loop as well
	// b/c we mark the camera as changed

	TransferFunctionWidget widget;

	AppState app_state;
	size_t frame_id = 0;
	glm::vec2 prev_mouse(-2.f);
	bool camera_changed = true;
	bool fbsize_changed = false;
	bool data_changed = false;
	while (!app_state.done) {
		SDL_Event event;
		while (SDL_PollEvent(&event)) {
			ImGui_ImplSDL2_ProcessEvent(&event);
			if (event.type == SDL_QUIT) {
				app_state.done = true;
			}
			if (!io.WantCaptureKeyboard && event.type == SDL_KEYDOWN) {
				if (event.key.keysym.sym == SDLK_ESCAPE) {
					app_state.done = true;
				} else if (event.key.keysym.sym == SDLK_p) {
					auto eye = camera.eye();
					auto center = camera.center();
					auto up = camera.up();
					std::cout << "-eye " << eye.x << " " << eye.y << " " << eye.z
						<< " -center " << center.x << " " << center.y << " " << center.z
						<< " -up " << up.x << " " << up.y << " " << up.z << "\n";
				}
			}
			if (event.type == SDL_WINDOWEVENT && event.window.event == SDL_WINDOWEVENT_CLOSE
					&& event.window.windowID == SDL_GetWindowID(window)) {
				app_state.done = true;
			}
			if (!io.WantCaptureMouse) {
				if (event.type == SDL_MOUSEMOTION) {
					const glm::vec2 cur_mouse = transform_mouse(glm::vec2(event.motion.x, event.motion.y));
					if (prev_mouse != glm::vec2(-2.f)) {
						if (event.motion.state & SDL_BUTTON_LMASK) {
							camera.rotate(prev_mouse, cur_mouse);
							camera_changed = true;
						} else if (event.motion.state & SDL_BUTTON_RMASK) {
							camera.pan(cur_mouse - prev_mouse);
							camera_changed = true;
						}
					}
					prev_mouse = cur_mouse;
				} else if (event.type == SDL_MOUSEWHEEL) {
					camera.zoom(event.wheel.y * 0.1);
					camera_changed = true;
				}
			}
			if (event.type == SDL_WINDOWEVENT && event.window.event == SDL_WINDOWEVENT_RESIZED) {
				fbsize_changed = true;

				win_width = event.window.data1;
				win_height = event.window.data2;
				io.DisplaySize.x = win_width;
				io.DisplaySize.y = win_height;

				glDeleteTextures(1, &render_texture);
				glGenTextures(1, &render_texture);
				// Setup the render textures for color and normals
				glBindTexture(GL_TEXTURE_2D, render_texture);
				glTexStorage2D(GL_TEXTURE_2D, 1, GL_RGBA8, win_width, win_height);

				glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
				glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
				glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
				glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
			}
		}

		int ready = 0;
		if (region_future.wait_for(std::chrono::milliseconds(0)) == std::future_status::ready) {
			ready = 1;
		}

		int all_new_state = 0;
		MPI_Allreduce(&ready, &all_new_state, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

		if (all_new_state) {
			data_changed = true;
			regions = region_future.get();
			OSPModel newworld = build_regions(regions, tfcn);
			ospSetObject(renderer, "model", newworld);
			ospCommit(renderer);
			// TODO: Releasing the world here shouldn't segfault
			//ospRelease(world);
		
			region_future = is::client::query_async();
		}

		app_state.cam_pos = camera.eye();
		app_state.cam_dir = camera.dir();
		app_state.cam_up = camera.up();
		app_state.fb_dims = glm::ivec2(win_width, win_height);
		app_state.camera_changed = camera_changed;
		app_state.fbsize_changed = fbsize_changed;
		app_state.tfcn_changed = widget.changed();
		MPI_Bcast(&app_state, sizeof(app_state), MPI_BYTE, 0, MPI_COMM_WORLD);

		if (app_state.tfcn_changed) {
			auto colormap = widget.get_colormap();
			size_t colormap_size = colormap.size();
			MPI_Bcast(&colormap_size, sizeof(colormap_size), MPI_BYTE, 0, MPI_COMM_WORLD);
			MPI_Bcast(colormap.data(), colormap.size(), MPI_BYTE, 0, MPI_COMM_WORLD);
			update_transfer_fcn(tfcn, colormap);
		}

		if (camera_changed || fbsize_changed) {
			if (fbsize_changed) {
				ospRelease(fb);
				fb = ospNewFrameBuffer(osp::vec2i{app_state.fb_dims.x, app_state.fb_dims.y},
						OSP_FB_SRGBA, OSP_FB_COLOR | OSP_FB_ACCUM);
			}

			ospSet3fv(ospcamera, "pos", &app_state.cam_pos.x);
			ospSet3fv(ospcamera, "dir", &app_state.cam_dir.x);
			ospSet3fv(ospcamera, "up", &app_state.cam_up.x);
			ospSet1f(ospcamera, "fovy", 65.f);
			ospSet1f(ospcamera, "aspect",
					static_cast<float>(app_state.fb_dims.x) / app_state.fb_dims.y);
			ospCommit(ospcamera);
			ospCommit(renderer);
		}
		if (data_changed || camera_changed || fbsize_changed) {
			ospFrameBufferClear(fb, OSP_FB_COLOR | OSP_FB_ACCUM);
		}

		ImGui_ImplOpenGL3_NewFrame();
		ImGui_ImplSDL2_NewFrame(window);
		ImGui::NewFrame();

		if (ImGui::Begin("Transfer Function")) {
			widget.draw_ui();
		}
		ImGui::End();

		ospRenderFrame(fb, renderer, OSP_FB_COLOR);

		// Rendering
		ImGui::Render();
		glViewport(0, 0, (int)io.DisplaySize.x, (int)io.DisplaySize.y);

		{
			const uint32_t *mapped = static_cast<const uint32_t*>(ospMapFrameBuffer(fb, OSP_FB_COLOR));
			glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, win_width, win_height, GL_RGBA,
					GL_UNSIGNED_BYTE, mapped);
			ospUnmapFrameBuffer(mapped, fb);
		}

		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		glUseProgram(display_render->program);
		glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);

		ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

		SDL_GL_SwapWindow(window);

		++frame_id;
		camera_changed = false;
		data_changed = false;
	}

	// Release the shader before we release the GL context
	display_render = nullptr;

	ImGui_ImplOpenGL3_Shutdown();
	ImGui_ImplSDL2_Shutdown();
	ImGui::DestroyContext();

	SDL_GL_DeleteContext(gl_context);
	SDL_DestroyWindow(window);
	SDL_Quit();
}
void run_worker(const std::vector<std::string>&) {
	OSPTransferFunction tfcn = ospNewTransferFunction("piecewise_linear");
	OSPFrameBuffer fb = ospNewFrameBuffer(osp::vec2i{win_width, win_height},
			OSP_FB_SRGBA, OSP_FB_COLOR | OSP_FB_ACCUM);

	// Query the first timestep from the simulation
	auto regions = is::client::query();
	log_regions(regions);
	OSPModel world = build_regions(regions, tfcn);

	// Start querying for the next timestep asynchronously
	auto region_future = is::client::query_async();

	OSPCamera ospcamera = ospNewCamera("perspective");
	OSPRenderer renderer = ospNewRenderer("mpi_raycast");
	ospSetObject(renderer, "model", world);
	ospSetObject(renderer, "camera", ospcamera);
	// Note: we commit for the first time in the render loop as well
	// b/c we mark the camera as changed

	AppState app_state;
	size_t frame_id = 0;
	glm::vec2 prev_mouse(-2.f);
	bool data_changed = false;
	while (!app_state.done) {
		int ready = 0;
		if (region_future.wait_for(std::chrono::milliseconds(0)) == std::future_status::ready) {
			ready = 1;
		}

		int all_new_state = 0;
		MPI_Allreduce(&ready, &all_new_state, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

		if (all_new_state) {
			data_changed = true;
			regions = region_future.get();
			OSPModel newworld = build_regions(regions, tfcn);
			ospSetObject(renderer, "model", newworld);
			ospCommit(renderer);
			// TODO: Releasing the world here shouldn't segfault
			//ospRelease(world);
		
			region_future = is::client::query_async();
		}

		MPI_Bcast(&app_state, sizeof(app_state), MPI_BYTE, 0, MPI_COMM_WORLD);

		if (app_state.tfcn_changed) {
			std::vector<uint8_t> colormap;
			size_t colormap_size = 0;
			MPI_Bcast(&colormap_size, sizeof(colormap_size), MPI_BYTE, 0, MPI_COMM_WORLD);
			colormap.resize(colormap_size);
			MPI_Bcast(colormap.data(), colormap.size(), MPI_BYTE, 0, MPI_COMM_WORLD);
			update_transfer_fcn(tfcn, colormap);
		}

		if (app_state.camera_changed || app_state.fbsize_changed) {
			if (app_state.fbsize_changed) {
				ospRelease(fb);
				fb = ospNewFrameBuffer(osp::vec2i{app_state.fb_dims.x, app_state.fb_dims.y},
						OSP_FB_SRGBA, OSP_FB_COLOR | OSP_FB_ACCUM);
			}

			ospSet3fv(ospcamera, "pos", &app_state.cam_pos.x);
			ospSet3fv(ospcamera, "dir", &app_state.cam_dir.x);
			ospSet3fv(ospcamera, "up", &app_state.cam_up.x);
			ospSet1f(ospcamera, "fovy", 65.f);
			ospSet1f(ospcamera, "aspect",
					static_cast<float>(app_state.fb_dims.x) / app_state.fb_dims.y);
			ospCommit(ospcamera);
			ospCommit(renderer);
		}
		if (data_changed || app_state.camera_changed || app_state.fbsize_changed) {
			ospFrameBufferClear(fb, OSP_FB_COLOR | OSP_FB_ACCUM);
		}

		ospRenderFrame(fb, renderer, OSP_FB_COLOR);

		++frame_id;
		data_changed = false;
	}
}
OSPModel build_regions(const std::vector<is::SimState> &regions, OSPTransferFunction tfcn) {
	OSPModel world = ospNewModel();
	if (regions.size() > 1) {
		std::cout << "WARNING: Multiple regions per-rank are not handled right now\n";
	}
	auto region = regions[0];
	ospSet1i(world, "id", rank);
	ospSet3fv(world, "region.lower", &region.local.min.x);
	ospSet3fv(world, "region.upper", &region.local.max.x);

	if (region.particles.numParticles > 0) {
		OSPData data = ospNewData(region.particles.array->numBytes(),
				OSP_UCHAR, region.particles.array->data(),
				OSP_DATA_SHARED_BUFFER);
		ospCommit(data);

		OSPGeometry spheres = ospNewGeometry("spheres");
		ospSetData(spheres, "spheres", data);
		ospSet1i(spheres, "bytes_per_sphere", region.particles.array->stride());
		
		OSPMaterial mat = ospNewMaterial2("scivis", "OBJMaterial");
		ospSet3f(mat, "Kd", 1.f, 1.f, 1.f / rank);
		ospCommit(mat);
		ospSetMaterial(spheres, mat);

		ospCommit(spheres);
		ospAddGeometry(world, spheres);
	}

	for (const auto &f : region.fields) {
		if (f.first != display_field) {
			continue;
		}

		auto &field = f.second;
		OSPVolume vol = ospNewVolume("shared_structured_volume");
		ospSet3i(vol, "dimensions", field.dims[0], field.dims[1], field.dims[2]);
		switch (field.dataType) {
			case UINT8: ospSetString(vol, "voxelType", "uchar"); break;
			case FLOAT: ospSetString(vol, "voxelType", "float"); break;
			case DOUBLE: ospSetString(vol, "voxelType", "double"); break;
			default: throw std::runtime_error("Invalid voxel type");
		}

		OSPData vol_data = ospNewData(field.array->numBytes(), OSP_UCHAR,
				field.array->data(), OSP_DATA_SHARED_BUFFER);
		ospCommit(vol_data);
		ospSetData(vol, "voxelData", vol_data);
		ospSet1f(vol, "samplingRate", 1.f);
		ospSet1i(vol, "adaptiveSampling", 0);

		ospSetObject(vol, "transferFunction", tfcn);
		ospCommit(vol);
		ospAddVolume(world, vol);
	}

	ospCommit(world);
	return world;
}
void log_regions(const std::vector<is::SimState> &regions) {
	for (int i = 0; i < world_size; ++i) {
		if (rank == 0) {
			std::cout << "Rank " << rank << " has " << regions.size() << " regions\n";
			// For each region we received, print out its data
			for (const auto &r : regions) {
				std::cout << "Region:\n"
					<< "world bounds: " << r.world << "\n"
					<< "local bounds: " << r.local << "\n"
					<< "ghost bounds: " << r.ghost << "\n"
					<< "Region has " << r.particles.numParticles
					<< " particles and " << r.fields.size() << " fields\n";
				std::cout << "Fields: {";
				for (const auto &f : r.fields) {
					std::cout << "(" << f.first << ", ";
					if (f.second.dataType == UINT8) {
						std::cout << "uint8";
					} else if (f.second.dataType == FLOAT) {
						std::cout << "float";
					} else if (f.second.dataType == DOUBLE) {
						std::cout << "double";
					} else {
						std::cout << "INVALID!";
					}
					std::cout << ", [" << f.second.dims[0] << ", "
						<< f.second.dims[1] << ", " << f.second.dims[2]
						<< "]), ";
				}
				std::cout << "}\n";
			}
			std::cout << "--------" << std::endl;
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
}
void update_transfer_fcn(OSPTransferFunction tfcn, const std::vector<uint8_t> &colormap) {
	std::vector<float> colors;
	std::vector<float> opacities;
	for (size_t i = 0; i < colormap.size() / 4; ++i) {
		colors.push_back(colormap[i * 4] / 255.f);
		colors.push_back(colormap[i * 4 + 1] / 255.f);
		colors.push_back(colormap[i * 4 + 2] / 255.f);
		opacities.push_back(colormap[i * 4 + 3] / 255.f);
	}
	OSPData colors_data = ospNewData(colors.size() / 3, OSP_FLOAT3, colors.data());
	ospCommit(colors_data);
	OSPData opacity_data = ospNewData(opacities.size(), OSP_FLOAT, opacities.data());
	ospCommit(opacity_data);

	//The value range here will be different from Will's code. It will need to match Timo's data.
	ospSetData(tfcn, "colors", colors_data);
	ospSetData(tfcn, "opacities", opacity_data);
	ospSet2f(tfcn, "valueRange", value_range.x, value_range.y);
	ospCommit(tfcn);
}

