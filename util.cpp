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

#include <cmath>
#include <vector>
#include "ospcommon/vec.h"
#include "util.h"

using namespace ospcommon;

AppState::AppState() : fbSize(1024), cameraChanged(false), quit(false),
  fbSizeChanged(false), tfcnChanged(false), timestepChanged(false),
  fieldChanged(false)
{}
void AppState::compute_cam_dirs(vec3f &cam_du, vec3f &cam_dv) const {
	const float fovy = 60.0;
	vec3f dir = normalize(vec3f(v[1].x, v[1].y, v[1].z));
	cam_du = normalize(cross(dir, vec3f(v[2].x, v[2].y, v[2].z)));
	cam_dv = normalize(cross(cam_du, dir));
	const float img_y = 2.0 * std::tan((M_PI / 180.f) * 0.5 * fovy);
	const float img_x = img_y * fbSize.x / static_cast<float>(fbSize.y);
	cam_du *= img_x;
	cam_dv *= img_y;
}

DistributedRegion::DistributedRegion(box3f bounds, int id)
	: bounds(bounds), id(id)
{}

