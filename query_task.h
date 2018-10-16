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

#include <thread>
#include <atomic>
#include <mutex>
#include <vector>
#include "ospray/ospray_cpp/Model.h"
#include "ospray/ospray_cpp/Geometry.h"

#include "libIS/is_render.h"

struct Sphere {
	ospcommon::vec3f pos;
	int atom_type;

	Sphere();
	Sphere(float x, float y, float z, int type);
	Sphere(ospcommon::vec3f v, int type);
};

struct SimulationState {
	std::vector<is::SimState> regions;
	std::vector<ospray::cpp::Geometry> geom;
	std::vector<ospray::cpp::Model> models, ghost_models;
	ospcommon::box3f world_bounds;

	SimulationState() = default;
	// TODO: This should be ok to do the cleanup here even on a
	// background thread since it's a local op. But will it be
	// thread safe with the local operations being done? It's
	// an independent model so should be ok.
	~SimulationState();
	SimulationState(const SimulationState &) = delete;
	SimulationState& operator=(const SimulationState &) = delete;
};

// Struct responsible for polling the sim data in the background and
// building the OSPRay geometry and model for the new timestep received
class QueryTask {
	std::mutex mutex;
	MPI_Comm query_comm;
	std::shared_ptr<SimulationState> available_state;
	std::atomic<bool> quitThread;
	std::thread thread;

public:
	QueryTask(MPI_Comm worker_comm);
	~QueryTask();
	// Take the loaded region, returns nullptr if there's no new
	// data from the simulation since the last time take was called.
	std::shared_ptr<SimulationState> take();
	void queryThread();

private:
	// Compute new colormap for the updated atom type attributes. Also re-scales so atom type min = 0
	std::vector<ospcommon::vec4f> colormap_atoms(std::vector<is::SimState> &regions,
			std::array<int, 2> &atom_type_range);
	ospray::cpp::Geometry make_atom_geometry(const is::SimState &region,
			const std::vector<ospcommon::vec4f> &atom_colors);
	ospray::cpp::Geometry make_bond_geometry(const is::SimState &region,
			const std::vector<ospcommon::vec4f> &atom_colors,
			const ospcommon::box3f &world_bounds,
			const std::array<int, 2> &atom_type_range);
};

