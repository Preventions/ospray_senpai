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
#include <cstring>
#include <stdexcept>
#ifndef _WIN32
#include <unistd.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netinet/ip.h>
#include <netdb.h>
#else
#include <winsock2.h>
#include <WS2tcpip.h>
#include <windows.h>
#endif
#include <atomic>
#include "socket_fabric.h"

#ifdef _WIN32
#include <winsock2.h>
#include <ws2tcpip.h>
#define SHUT_RDWR SD_BOTH
#endif

#ifdef _WIN32

struct WinSockContext {
	WinSockContext() {
		WSADATA wsa_data;
		if (WSAStartup(MAKEWORD(2, 2), &wsa_data) != 0) {
			throw std::runtime_error("Failed to initialize WinSock!");
		}
	}
	~WinSockContext() {
		WSACleanup();
	}
};

static std::unique_ptr<WinSockContext> winsock = nullptr;

void initialize() {
	static std::atomic<bool> initialized = false;
	if (!initialized) {
		static std::mutex mutex;
		std::lock_guard<std::mutex> lock(mutex);
		if (!winsock) {
			winsock = std::make_unique<WinSockContext>();
		}
		initialized = true;
	}
}
#else
void initialize() {}
#endif

static void close_socket(int socket) {
	if (socket > -1) {
		shutdown(socket, SHUT_RDWR);
#ifdef _WIN32
		closesocket(socket);
#else
		close(socket);
#endif
	}
}

SocketFabric::SocketFabric(const uint16_t port) : buffer(4096, 0) {
	initialize();
	int listen_socket = socket(AF_INET, SOCK_STREAM, 0);
#ifdef SO_REUSEADDR
	{
		int flag = 1;
		setsockopt(listen_socket, SOL_SOCKET, SO_REUSEADDR, (char*)&flag, sizeof(int));
	}
#endif
	struct sockaddr_in serv_addr = {0};
	serv_addr.sin_family = AF_INET;
	serv_addr.sin_port = htons(port);
	serv_addr.sin_addr.s_addr = INADDR_ANY;

	if (bind(listen_socket, (struct sockaddr*)&serv_addr, sizeof(serv_addr)) < 0) {
		throw std::runtime_error("Failed to bind to socket");
	}

	if (listen(listen_socket, 1) < 0) {
		throw std::runtime_error("Failed to listen on socket");
	}

	// Accept an incoming connection
	struct sockaddr_in addr;
	socklen_t len = sizeof(addr);
	sock_fd = accept(listen_socket, (struct sockaddr*)&addr, &len);
	if (sock_fd == -1) {
		throw std::runtime_error("Failed to accept connection");
	}

#ifdef TCP_NODELAY
	{
		int flag = 1; 
		setsockopt(sock_fd, IPPROTO_TCP, TCP_NODELAY, (char*)&flag, sizeof(int));
	}
#endif

#ifdef SO_NOSIGPIPE
	{
		int flag = 1;
		setsockopt(sock_fd, SOL_SOCKET, SO_NOSIGPIPE, (void*)&flag, sizeof(int));
	}
#endif
	close_socket(listen_socket);
}
SocketFabric::SocketFabric(const std::string &hostname, const uint16_t port) : buffer(4096, 0) {
	initialize();
	sock_fd = socket(AF_INET, SOCK_STREAM, 0);

	struct hostent* server = gethostbyname(hostname.c_str());
	if (server == nullptr) {
		throw std::runtime_error("Hostname " + hostname + " not found!");
	}

	struct sockaddr_in serv_addr = {0};
	serv_addr.sin_family = AF_INET;
	serv_addr.sin_port = htons(port);
	std::memcpy((char*)&serv_addr.sin_addr.s_addr, (char*)server->h_addr, server->h_length);

	if (connect(sock_fd, (struct sockaddr*)&serv_addr, sizeof(serv_addr)) < 0) {
		throw std::runtime_error("Connection to " + hostname + " on port " + std::to_string(port) + " failed");
	}

#ifdef TCP_NODELAY
	{
		int flag = 1; 
		setsockopt(sock_fd, IPPROTO_TCP, TCP_NODELAY, (char*)&flag, sizeof(int));
	}
#endif

#ifdef SO_NOSIGPIPE
	{
		int flag = 1;
		setsockopt(sock_fd, SOL_SOCKET, SO_NOSIGPIPE, (void*)&flag, sizeof(int));
	}
#endif
}
SocketFabric::~SocketFabric() {
	close_socket(sock_fd);
}

void SocketFabric::send(void *mem, size_t s) {
#ifndef _WIN32
	if (write(sock_fd, mem, s) != s) {
		throw std::runtime_error("Failed to write all bytes");
	}
#else
	if (::send(sock_fd, (const char*)mem, s, 0) == SOCKET_ERROR) {
		throw std::runtime_error("Failed to write all bytes");
	}
#endif
}
size_t SocketFabric::read(void *&mem) {
	if (buffer.empty()) {
		buffer.resize(4096, 0);
	}
#ifndef _WIN32
	const size_t s = ::read(sock_fd, buffer.data(), buffer.size());
#else
	const size_t s = recv(sock_fd, buffer.data(), buffer.size(), 0);
#endif
	mem = buffer.data();
	return s;
}

