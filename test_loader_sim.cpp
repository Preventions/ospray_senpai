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

#include <chrono>
#include <iostream>
#include <cstring>
#include <thread>
#include <random>
#include <vector>
#include <mpi.h>
#include <libIS/is_sim.h>
#include <ospcommon/FileName.h>
#include <ospcommon/box.h>

using namespace ospcommon;

struct Particle {
	vec3f pos;
	float radius;

	Particle(const vec3f &p, float radius) : pos(p), radius(radius) {}
};

int rank, world_size;

void load_cosmic_web(const FileName &filename, std::vector<Particle> &particles,
    box3f &local_bounds);

int main(int ac, char **av) {
	MPI_Init(&ac, &av);
	std::vector<std::string> args(av, av + ac);

	std::vector<FileName> cosmic_web_bricks;
	int n_steps = 100;
	bool mpi_multilaunch = false;
	for (int i = 1; i < ac; ++i) {
		if (args[i] == "-step") {
			n_steps = std::atoi(av[++i]);
		} else if (args[i] == "-mpi-multi") {
			mpi_multilaunch = true;
		} else {
			cosmic_web_bricks.push_back(args[i]);
		}
	}

	MPI_Comm sim_comm = MPI_COMM_WORLD;
	if (mpi_multilaunch) {
		const int test_tag = 0x54455354;
		MPI_Comm_split(MPI_COMM_WORLD, test_tag, 0, &sim_comm);
	}
	MPI_Comm_rank(sim_comm, &rank);
	MPI_Comm_size(sim_comm, &world_size);

	std::cout << "#sim rank " << rank << "/" << world_size << "\n";

	if (mpi_multilaunch) {
		std::cout << "Connecting with existing comm" << std::endl;
		libISInitWithExisting(sim_comm, MPI_COMM_WORLD);
	} else {
		std::cout << "Waiting for network connection" << std::endl;
		libISInit(MPI_COMM_WORLD, 29374);
	}

	box3f bounds;
	const size_t bricks_per_rank = cosmic_web_bricks.size() / world_size;
	std::vector<Particle> particles;
	for (size_t i = 0; i < bricks_per_rank; ++i) {
		load_cosmic_web(cosmic_web_bricks[rank * bricks_per_rank + i],
				particles, bounds);
	}
	std::cout << "Loaded " << particles.size() << " particles in total\n";
	std::cout << "Data size to send " << particles.size() * sizeof(Particle) << " bytes\n";
	libISBox3f local;
	local.min.x = bounds.lower.x;
	local.min.y = bounds.lower.y;
	local.min.z = bounds.lower.z;
	local.max.x = bounds.upper.x;
	local.max.y = bounds.upper.y;
	local.max.z = bounds.upper.z;

	// TODO: Need to collect the world bounds
	libISSimState *state = libISMakeSimState();
	//libISSetWorldBounds(state, world);
	libISSetLocalBounds(state, local);
	libISSetGhostBounds(state, local);
	libISSetParticles(state, particles.size(), 0, sizeof(Particle), particles.data());

	for (int i = 0; i < n_steps; ++i) {
		MPI_Barrier(sim_comm);
		std::this_thread::sleep_for(std::chrono::seconds(2));
		if (rank == 0) {
			std::cout << "Timestep " << i << "\n";
		}
		libISProcess(state);
	}

	libISFreeSimState(state);
	libISFinalize();

	if (mpi_multilaunch) {
		MPI_Comm_free(&sim_comm);
	}
	MPI_Finalize();
	return 0;
}

#pragma pack(1)
struct CosmicWebHeader {
  // number of particles in this dat file
  int np_local;
  float a, t, tau;
  int nts;
  float dt_f_acc, dt_pp_acc, dt_c_acc;
  int cur_checkpoint, cur_projection, cur_halofind;
  float massp;
};
std::ostream& operator<<(std::ostream &os, const CosmicWebHeader &h) {
  os << "{\n\tnp_local = " << h.np_local
    << "\n\ta = " << h.a
    << "\n\tt = " << h.t
    << "\n\ttau = " << h.tau
    << "\n\tnts = " << h.nts
    << "\n\tdt_f_acc = " << h.dt_f_acc
    << "\n\tdt_pp_acc = " << h.dt_pp_acc
    << "\n\tdt_c_acc = " << h.dt_c_acc
    << "\n\tcur_checkpoint = " << h.cur_checkpoint
    << "\n\tcur_halofind = " << h.cur_halofind
    << "\n\tmassp = " << h.massp
    << "\n}";
  return os;
}

void load_cosmic_web(const FileName &filename, std::vector<Particle> &particles,
    box3f &local_bounds)
{
  std::ifstream fin(filename.c_str(), std::ios::binary);

  if (!fin.good()) {
    throw std::runtime_error("could not open particle data file " + filename.str());
  }

  CosmicWebHeader header;
  if (!fin.read(reinterpret_cast<char*>(&header), sizeof(CosmicWebHeader))) {
    throw std::runtime_error("Failed to read header");
  }

  std::cout << "Cosmic Web Header: " << header << "\n";

  // Compute the brick offset for this file, given in the last 3 numbers of the name
  std::string brick_name = filename.name();
  brick_name = brick_name.substr(brick_name.size() - 3, 3);
  const int brick_number = std::stoi(brick_name);
  // The cosmic web bricking is 8^3
  const int brick_z = brick_number / 64;
  const int brick_y = (brick_number / 8) % 8;
  const int brick_x = brick_number % 8;
  std::cout << "Brick position = { " << brick_x << ", " << brick_y
    << ", " << brick_z << " }\n";
  // Each cell is 768x768x768 units
  const float step = 768.f;
  const vec3f offset(step * brick_x, step * brick_y, step * brick_z);

  particles.reserve(particles.size() + header.np_local);
  for (int i = 0; i < header.np_local; ++i) { 
    vec3f position, velocity;

    if (!fin.read(reinterpret_cast<char*>(&position), sizeof(vec3f))) {
      throw std::runtime_error("Failed to read position for particle");
    }
    if (!fin.read(reinterpret_cast<char*>(&velocity), sizeof(vec3f))) {
      throw std::runtime_error("Failed to read velocity for particle");
    }
    position += offset;

    local_bounds.lower = min(position, local_bounds.lower);
    local_bounds.upper = max(position, local_bounds.upper);
    particles.emplace_back(position, 1.0);
  }
}


