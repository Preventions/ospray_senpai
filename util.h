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

#include <array>
#include <vector>
#include "ospcommon/vec.h"
#include "ospcommon/box.h"

using vec3sz = ospcommon::vec_t<size_t, 3>;

// Struct for bcasting out the camera change info and general app state
#pragma pack(1)
struct AppState {
  // eye pos, look dir, up dir
  std::array<ospcommon::vec3f, 3> v;
  ospcommon::vec2i fbSize;
  uint64_t currentTimestep;
  bool cameraChanged, quit, fbSizeChanged,
       tfcnChanged, timestepChanged, fieldChanged;

  AppState();
  void compute_cam_dirs(ospcommon::vec3f &cam_du, ospcommon::vec3f &cam_dv) const;
};

// Struct for holding the other app data buffers and info that
// we can't bcast directly.
// TODO: These we don't really need for the lammps in situ stuff
struct AppData {
  std::string currentVariable;
  std::vector<ospcommon::vec3f> tfcn_colors;
  std::vector<float> tfcn_alphas;
};

struct DistributedRegion {
	ospcommon::box3f bounds;
	int id;

	DistributedRegion(ospcommon::box3f bounds = ospcommon::empty, int id = -1);
};

