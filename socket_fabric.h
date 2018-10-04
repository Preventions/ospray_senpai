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

#include <cstdint>
#include <string>
#include <vector>
#include "ospcommon/networking/Fabric.h"

struct SocketFabric : public ospcommon::networking::Fabric {
  int sock_fd;
  std::vector<char> buffer;

public:
  SocketFabric(const uint16_t port);
  SocketFabric(const std::string &hostname, const uint16_t port);
  ~SocketFabric();

  SocketFabric(const SocketFabric&) = delete;
  SocketFabric& operator=(const SocketFabric&) = delete;

  virtual void send(void *mem, size_t s) override;
  virtual size_t read(void *&mem) override;
};

