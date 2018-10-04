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

#include "colormap.h"

using namespace ospcommon;

Colormap::Colormap() {
	// Jet
	colors.push_back(vec3f(0       , 0, 0.562493));
	colors.push_back(vec3f(0       , 0, 1       ));
	colors.push_back(vec3f(0       , 1, 1       ));
	colors.push_back(vec3f(0.500008, 1, 0.500008));
	colors.push_back(vec3f(1       , 1, 0       ));
	colors.push_back(vec3f(1       , 0, 0       ));
	colors.push_back(vec3f(0.500008, 0, 0       ));

	/*
	// Ice/Fire
	colors.push_back(vec3f(0        , 0          , 0          ));
	colors.push_back(vec3f(0        , 0.120394   , 0.302678   ));
	colors.push_back(vec3f(0        , 0.216587   , 0.524575   ));
	colors.push_back(vec3f(0.0552529, 0.345022   , 0.659495   ));
	colors.push_back(vec3f(0.128054 , 0.492592   , 0.720287   ));
	colors.push_back(vec3f(0.188952 , 0.641306   , 0.792096   ));
	colors.push_back(vec3f(0.327672 , 0.784939   , 0.873426   ));
	colors.push_back(vec3f(0.60824  , 0.892164   , 0.935546   ));
	colors.push_back(vec3f(0.881376 , 0.912184   , 0.818097   ));
	colors.push_back(vec3f(0.9514   , 0.835615   , 0.449271   ));
	colors.push_back(vec3f(0.904479 , 0.690486   , 0          ));
	colors.push_back(vec3f(0.854063 , 0.510857   , 0          ));
	colors.push_back(vec3f(0.777096 , 0.330175   , 0.000885023));
	colors.push_back(vec3f(0.672862 , 0.139086   , 0.00270085 ));
	colors.push_back(vec3f(0.508812 , 0          , 0          ));
	colors.push_back(vec3f(0.299413 , 0.000366217, 0.000549325));
	colors.push_back(vec3f(0.0157473, 0.00332647 , 0          ));
	*/

	/*
	// Cool/Warm
	colors.push_back(vec3f(0.231373, 0.298039 , 0.752941));
	colors.push_back(vec3f(0.865003, 0.865003 , 0.865003));
	colors.push_back(vec3f(0.705882, 0.0156863, 0.14902));

	// Blue/Red
	colors.push_back(vec3f(0, 0, 1));
	colors.push_back(vec3f(1, 0, 0));

	// Grayscale
	colors.push_back(vec3f(0));
	colors.push_back(vec3f(1));
	*/
}
vec3f Colormap::map(const float x) const {
	if (!std::isnormal(x)) {
		return colors[0];
	}
	const int lo = std::floor(x * colors.size());
	const int hi = std::min(lo + 1, int(colors.size() - 1));
	const float rem = x - lo;
	if (hi == 0) {
		return colors[0];
	} else if (lo == colors.size()) {
		return colors.back();
	}
	return (1.f - rem) * colors[lo] + rem * colors[hi];
}

