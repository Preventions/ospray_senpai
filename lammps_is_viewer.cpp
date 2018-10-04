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
#include <iostream>
#include <algorithm>
#include <vector>
#include <array>
#include <chrono>
#include <functional>
#ifdef _WIN32
#include <winsock2.h>
#include <windows.h>
#endif

#include <turbojpeg.h>
#include <GL/gl3w.h>
#include <GLFW/glfw3.h>

#include "imgui/imgui.h"
#include "imgui_impl_glfw_gl3.h"

#include "gldebug.h"
#include "glshaders.h"
#include "arcball.h"
#include "util.h"
#include "image_util.h"
#include "client_server.h"

using namespace ospcommon;
using namespace std::chrono;

bool save_screenshot = false;

// Extra stuff we need in GLFW callbacks
struct WindowState {
	Arcball &camera;
	vec2f prev_mouse;
	bool camera_changed;
	AppState &app;

	WindowState(AppState &app, Arcball &camera)
		: camera(camera), prev_mouse(-1), camera_changed(true), app(app)
	{}
};

void keyCallback(GLFWwindow *window, int key, int scancode, int action, int mods) {
	if (ImGui::GetIO().WantCaptureKeyboard) {
		return;
	}

	WindowState *state = static_cast<WindowState*>(glfwGetWindowUserPointer(window));
	if (action == GLFW_PRESS) {
		switch (key) {
			case GLFW_KEY_ESCAPE:
				glfwSetWindowShouldClose(window, true);
				break;
			case 'P':
			case 'p':
				save_screenshot = true;
				break;
			default:
				break;
		}
	}
	// Forward on to ImGui
	ImGui_ImplGlfwGL3_KeyCallback(window, key, scancode, action, mods);
}

void cursorPosCallback(GLFWwindow *window, double x, double y) {
	if (ImGui::GetIO().WantCaptureMouse) {
		return;
	}

	WindowState *state = static_cast<WindowState*>(glfwGetWindowUserPointer(window));

	const vec2f mouse(x, y);
	if (state->prev_mouse != vec2f(-1)) {
		const bool leftDown =
			glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS;
		const bool rightDown =
			glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT) == GLFW_PRESS;
		const bool middleDown =
			glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_MIDDLE) == GLFW_PRESS;
		const vec2f prev = state->prev_mouse;

		if (leftDown) {
			const vec2f mouseFrom(clamp(prev.x * 2.f / state->app.fbSize.x - 1.f,  -1.f, 1.f),
					clamp(1.f - 2.f * prev.y / state->app.fbSize.y, -1.f, 1.f));
			const vec2f mouseTo  (clamp(mouse.x * 2.f / state->app.fbSize.x - 1.f,  -1.f, 1.f),
					clamp(1.f - 2.f * mouse.y / state->app.fbSize.y, -1.f, 1.f));
			state->camera.rotate(mouseFrom, mouseTo);
			state->camera_changed = true;
		} else if (middleDown || (rightDown && glfwGetKey(window, GLFW_KEY_W))) {
			const vec2f mouseFrom(clamp(prev.x * 2.f / state->app.fbSize.x - 1.f,  -1.f, 1.f),
					clamp(1.f - 2.f * prev.y / state->app.fbSize.y, -1.f, 1.f));
			const vec2f mouseTo   (clamp(mouse.x * 2.f / state->app.fbSize.x - 1.f,  -1.f, 1.f),
					clamp(1.f - 2.f * mouse.y / state->app.fbSize.y, -1.f, 1.f));
			const vec2f mouseDelta = mouseTo - mouseFrom;
			state->camera.pan(mouseDelta);
			state->camera_changed = true;
		} else if (rightDown) {
			state->camera.zoom(mouse.y - prev.y);
			state->camera_changed = true;
		}
	}
	state->prev_mouse = mouse;
}

void framebufferSizeCallback(GLFWwindow *window, int width, int height) {
	WindowState *state = static_cast<WindowState*>(glfwGetWindowUserPointer(window));
	state->app.fbSize = vec2i(width, height);
	state->app.fbSizeChanged = true;
}

void charCallback(GLFWwindow*, unsigned int c) {
	ImGuiIO& io = ImGui::GetIO();
	if (c > 0 && c < 0x10000) {
		io.AddInputCharacter((unsigned short)c);
	}
}

void app_loop(const std::string &serverhost, const int port, GLFWwindow *window);
int main(int argc, const char **argv) {
	std::vector<std::string> args(argv, argv + argc);

	std::string serverhost;
	int port = -1;
	for (int i = 1; i < argc; ++i) {
		if (args[i] == "-server") {
			serverhost = argv[++i];
		} else if (args[i] == "-port") {
			port = std::atoi(argv[++i]);
		}
	}
	if (serverhost.empty() || port < 0) {
		throw std::runtime_error("Usage: ./lammps_is_viewer -server <server host> -port <port>");
	}


	if (!glfwInit()) {
		std::cout << "Failed to init GLFW\n";
		return 1;
	}

	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	//glfwWindowHint(GLFW_OPENGL_DEBUG_CONTEXT, GL_TRUE);

	GLFWwindow *window = glfwCreateWindow(1024, 1024,
			"In Situ LAMMPS OSPRay Remote Viewer", nullptr, nullptr);

	if (!window) {
		glfwTerminate();
		std::cout << "Failed to open window\n";
		return 1;
	}
	glfwMakeContextCurrent(window);
	glfwSwapInterval(1);

	if (gl3wInit()) {
		std::cout << "Failed to init gl3w\n";
		glfwTerminate();
		std::exit(1);
	}
	//register_debug_callback();

	ImGui_ImplGlfwGL3_Init(window, false);

	glfwSetKeyCallback(window, keyCallback);
	glfwSetCursorPosCallback(window, cursorPosCallback);
	glfwSetFramebufferSizeCallback(window, framebufferSizeCallback);
	glfwSetMouseButtonCallback(window, ImGui_ImplGlfwGL3_MouseButtonCallback);
	glfwSetScrollCallback(window, ImGui_ImplGlfwGL3_ScrollCallback);
	glfwSetCharCallback(window, charCallback);

#ifdef _WIN32
	WSADATA wsa_data;
	if (WSAStartup(MAKEWORD(2, 2), &wsa_data) != 0) {
		std::cout << "Failed to start winsock" << std::endl;
		throw std::runtime_error("Failed to start winsock");
	}
#endif

	app_loop(serverhost, port, window);

#ifdef _WIN32
	WSACleanup();
#endif

	ImGui_ImplGlfwGL3_Shutdown();
	glfwDestroyWindow(window);
	return 0;
}
void app_loop(const std::string &serverhost, const int port, GLFWwindow *window) {
	AppState app;
	//app.fbSize = vec2i(1920, 1080);
	//app.fbSizeChanged = true;

	AppData appdata;
	bool got_world_bounds = false;
	box3f world_bounds(vec3f(-1), vec3f(1));
	Arcball arcball_camera(world_bounds);
	auto window_state = std::make_shared<WindowState>(app, arcball_camera);
	glfwSetWindowUserPointer(window, window_state.get());

	std::vector<uint32_t> imgBuf;
	std::vector<float> depth_buf;
	int frameTime = 0;

	std::array<GLuint, 2> osp_upload_tex;
	glGenTextures(osp_upload_tex.size(), osp_upload_tex.data());
	for (const auto &c : osp_upload_tex) {
		glBindTexture(GL_TEXTURE_2D, c);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);
	}
	glBindTexture(GL_TEXTURE_2D, osp_upload_tex[0]);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, app.fbSize.x, app.fbSize.y, 0,
			GL_RGBA, GL_UNSIGNED_BYTE, nullptr);
	glBindTexture(GL_TEXTURE_2D, osp_upload_tex[1]);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, app.fbSize.x, app.fbSize.y, 0,
			GL_RED, GL_FLOAT, nullptr);

	GLuint blit_from_fbo;
	glGenFramebuffers(1, &blit_from_fbo);
	glBindFramebuffer(GL_FRAMEBUFFER, blit_from_fbo);
	glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0,
			GL_TEXTURE_2D, osp_upload_tex[0], 0);

	// Mark the camera changed to init the uniforms
	window_state->camera_changed = true;
	app.cameraChanged = true;

	// Compressor for saving screenshots
	JPGCompressor compressor(100, true);

	std::vector<unsigned char> jpgBuf;
	JPGDecompressor decompressor;
	ServerConnection server(serverhost, port, app);
	while (!app.quit)  {
		ImGui_ImplGlfwGL3_NewFrame();

		if (server.get_new_frame(jpgBuf, depth_buf, frameTime)) {
			decompressor.decompress(jpgBuf.data(), jpgBuf.size(), app.fbSize.x,
					app.fbSize.y, imgBuf);

			glActiveTexture(GL_TEXTURE0);
			glBindTexture(GL_TEXTURE_2D, osp_upload_tex[0]);
			glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, app.fbSize.x, app.fbSize.y, 0,
					GL_RGBA, GL_UNSIGNED_BYTE, imgBuf.data());

			glActiveTexture(GL_TEXTURE1);
			glBindTexture(GL_TEXTURE_2D, osp_upload_tex[1]);
			glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, app.fbSize.x, app.fbSize.y, 0,
					GL_RED, GL_FLOAT, depth_buf.data());
		}

		glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
		glBindFramebuffer(GL_READ_FRAMEBUFFER, blit_from_fbo);
		glBlitFramebuffer(0, 0, app.fbSize.x, app.fbSize.y,
				0, 0, app.fbSize.x, app.fbSize.y,
				GL_COLOR_BUFFER_BIT, GL_NEAREST);

		if (save_screenshot) {
			save_screenshot = false;
			auto jpg = compressor.compress(imgBuf.data(), app.fbSize.x, app.fbSize.y);
			std::ofstream fout("screenshot.jpg", std::ios::binary);
			fout.write(reinterpret_cast<const char*>(jpg.first), jpg.second);

			std::cout << "Screenshot saved to 'screenshot.jpg'" << std::endl;
		}

		if (!got_world_bounds) {
			ImGui::Text("Waiting for server to load");
		}
		ImGui::Text("Render Worker: frame took %dms", frameTime);
		ImGui::Text("Client Application: %.3f ms/frame (%.1f FPS)",
				1000.0f / ImGui::GetIO().Framerate,
				ImGui::GetIO().Framerate);

		ImGui::Render();

		glfwSwapBuffers(window);

		glfwPollEvents();
		if (glfwWindowShouldClose(window)) {
			app.quit = true;
		}

		if (!got_world_bounds && server.get_world_bounds(world_bounds)) {
			std::cout << "Got world bounds = " << world_bounds << std::endl;
			got_world_bounds = true;
			arcball_camera = Arcball(world_bounds);
			window_state->camera_changed = true;
		}

		const vec3f eye = window_state->camera.eyePos();
		const vec3f look = window_state->camera.lookDir();
		const vec3f up = window_state->camera.upDir();
		app.v[0] = vec3f(eye.x, eye.y, eye.z);
		app.v[1] = vec3f(look.x, look.y, look.z);
		app.v[2] = vec3f(up.x, up.y, up.z);
		app.cameraChanged = window_state->camera_changed;
		window_state->camera_changed = false;

		server.update_app_state(app, appdata);

		if (app.fbSizeChanged) {
			imgBuf.resize(app.fbSize.x * app.fbSize.y, 0);
			depth_buf.resize(app.fbSize.x * app.fbSize.y, 0.f);
			app.fbSizeChanged = false;
			glBindTexture(GL_TEXTURE_2D, osp_upload_tex[0]);
			glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, app.fbSize.x, app.fbSize.y, 0,
					GL_RGBA, GL_UNSIGNED_BYTE, nullptr);
			glBindTexture(GL_TEXTURE_2D, osp_upload_tex[1]);
			glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, app.fbSize.x, app.fbSize.y, 0,
					GL_RED, GL_FLOAT, nullptr);

			glViewport(0, 0, app.fbSize.x, app.fbSize.y);
		}
	}

	glDeleteTextures(osp_upload_tex.size(), osp_upload_tex.data());
	glDeleteFramebuffers(1, &blit_from_fbo);
}

