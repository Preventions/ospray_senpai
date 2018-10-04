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

#pragma once

#include <vector>
#include <atomic>
#include <thread>
#include <set>
#include <mutex>
#include "ospcommon/networking/Socket.h"
#include "ospcommon/box.h"
#include "ospcommon/networking/BufferedDataStreaming.h"
#include "util.h"
#include "image_util.h"
#include "socket_fabric.h"

// A connection to the render worker server
class ServerConnection {
	std::string server_host;
	int server_port;

	std::vector<unsigned char> jpg_buf;
	std::vector<float> depth_buf;
	int frame_time;
	bool new_frame;
	std::mutex frame_mutex;

	AppState app_state;
	AppData app_data;
	std::mutex state_mutex;

	std::atomic<bool> have_world_bounds;
	ospcommon::box3f world_bounds;

	std::thread server_thread;

public:
	ServerConnection(const std::string &server, const int port,
			const AppState &app_state);
	~ServerConnection();
	/* Get the new JPG recieved from the network, if we've got a new one,
	 * otherwise the buf is unchanged.
	 */
	bool get_new_frame(std::vector<unsigned char> &buf, std::vector<float> &depth, int &frame_time);
	// Get the world bounds of the dataset, if we've received it from the server
	bool get_world_bounds(ospcommon::box3f &bounds);
	// Update the app state to be sent over the network for the next frame
	void update_app_state(const AppState &state, const AppData &data);

private:
	void connection_thread();
};

// A client connecting to the render worker server
class ClientConnection {
	JPGCompressor compressor;
	SocketFabric fabric;
	ospcommon::networking::BufferedReadStream read_stream;
	ospcommon::networking::BufferedWriteStream write_stream;

public:
	ClientConnection(const int port);
	void send_metadata(const ospcommon::box3f &world_bounds);
	void send_frame(uint32_t *color, float *depth, int width, int height, int frame_time);
	void recieve_app_state(AppState &app, AppData &data);
};

