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

#include <iostream>
#include <chrono>
#include "client_server.h"

ServerConnection::ServerConnection(const std::string &server, const int port,
		const AppState &app_state)
	: server_host(server), server_port(port), new_frame(false), app_state(app_state),
	have_world_bounds(false)
{
	server_thread = std::thread([&](){ connection_thread(); });
}
ServerConnection::~ServerConnection() {
	{
		std::lock_guard<std::mutex> lock(state_mutex);
		app_state.quit = true;
	}
	server_thread.join();
}
bool ServerConnection::get_new_frame(std::vector<unsigned char> &buf, std::vector<float> &depth, int &time) {
	std::lock_guard<std::mutex> lock(frame_mutex);
	if (new_frame) {
		buf = jpg_buf;
		depth = depth_buf;
		time = frame_time;
		new_frame = false;
		return true;
	}
	return false;
}
bool ServerConnection::get_world_bounds(ospcommon::box3f &bounds) {
	if (have_world_bounds) {
		bounds = world_bounds;
		return true;
	}
	return false;
}
void ServerConnection::update_app_state(const AppState &state, const AppData &data) {
	std::lock_guard<std::mutex> lock(state_mutex);

	app_state.v = state.v;
	app_state.fbSize = state.fbSize;
	if (state.timestepChanged) {
		app_state.currentTimestep = state.currentTimestep;
	}
	app_state.cameraChanged = state.cameraChanged ? state.cameraChanged : app_state.cameraChanged;
	app_state.fbSizeChanged = state.fbSizeChanged ? state.fbSizeChanged : app_state.fbSizeChanged;
	app_state.tfcnChanged = state.tfcnChanged ? state.tfcnChanged : app_state.tfcnChanged;
	app_state.timestepChanged = state.timestepChanged ? state.timestepChanged : app_state.timestepChanged;
	app_state.fieldChanged = state.fieldChanged ? state.fieldChanged : app_state.fieldChanged;

	if (state.fieldChanged) {
		app_data.currentVariable = data.currentVariable;
	}
	app_data.tfcn_colors = data.tfcn_colors;
	app_data.tfcn_alphas = data.tfcn_alphas;
}
void ServerConnection::connection_thread() {
	SocketFabric fabric(server_host, server_port);
	ospcommon::networking::BufferedReadStream read_stream(fabric);
	ospcommon::networking::BufferedWriteStream write_stream(fabric);

	// Get the world bounds from the server
	read_stream >> world_bounds;
	have_world_bounds = true;

	AppState state;
	AppData data;
	std::vector<unsigned char> recv_jpg;
	std::vector<float> recv_depth;
	while (true) {
		// Receive a frame from the server
		uint64_t jpg_size = 0;
		read_stream >> jpg_size;
		recv_jpg.resize(jpg_size, 0);
		read_stream.read(recv_jpg.data(), jpg_size);

		uint64_t depth_size = 0;
		read_stream >> depth_size;
		recv_depth.resize(depth_size, std::numeric_limits<float>::infinity());
		read_stream.read(recv_depth.data(), depth_size * sizeof(float));

		int recv_frame_time = 0;
		read_stream >> recv_frame_time;
		{
			std::lock_guard<std::mutex> lock(frame_mutex);
			jpg_buf = recv_jpg;
			depth_buf = recv_depth;
			frame_time = recv_frame_time;
			new_frame = true;
		}

		// Send over the latest app state
		{
			std::lock_guard<std::mutex> lock(state_mutex);
			state = app_state;
			data = app_data;
		}
		write_stream.write(&app_state, sizeof(AppState));
		if (app_state.fieldChanged) {
			write_stream << app_data.currentVariable;
		}
		if (app_state.tfcnChanged) {
			write_stream << app_data.tfcn_colors << app_data.tfcn_alphas;
		}
		write_stream.flush();

		if (app_state.quit) {
			return;
		}
		app_state.tfcnChanged = false;
		app_state.cameraChanged = false;
		app_state.fbSizeChanged = false;
		app_state.timestepChanged = false;
		app_state.fieldChanged = false;
	}
}

ClientConnection::ClientConnection(const int port)
	: compressor(90), fabric(port), read_stream(fabric), write_stream(fabric)
{}
void ClientConnection::send_metadata(const ospcommon::box3f &world_bounds) {
	write_stream << world_bounds;
}
void ClientConnection::send_frame(uint32_t *img, float *depth, int width, int height, int frame_time) {
	auto jpg = compressor.compress(img, width, height);
	write_stream << uint64_t(jpg.second);
	write_stream.write(jpg.first, jpg.second);

	const uint64_t depth_size = width * height;
	write_stream << depth_size;
	write_stream.write(depth, depth_size * sizeof(float));

	write_stream << frame_time;
	write_stream.flush();
}
void ClientConnection::recieve_app_state(AppState &app, AppData &data) {
	read_stream.read(&app, sizeof(AppState));
	if (app.fieldChanged) {
		read_stream >> data.currentVariable;
	}
	if (app.tfcnChanged) {
		read_stream >> data.tfcn_colors >> data.tfcn_alphas;
	}
}

