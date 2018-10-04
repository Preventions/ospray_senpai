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

#include <fstream>
#include <iostream>
#include <vector>

#include <mpi.h>
#include "library.h"
#include "modify.h"
#include "fix.h"
#include "fix_external.h"
#include "lammps_driver.h"
#include "libIS/is_sim.h"

using namespace LAMMPS_NS;

void proxy_callback(void *ptr, bigint ntimestep, int nlocal, int *id,
		double **x, double **f);

LammpsDriver::LammpsDriver(const std::string &lammps_input, MPI_Comm world)
	: world_comm(world), lammps_input(lammps_input), lammps(nullptr), quit_thread(false)
{
	MPI_Comm_dup(world_comm, &lammps_comm);
	MPI_Comm_rank(world_comm, &rank);
	lammps_open(0, nullptr, lammps_comm, (void**)&lammps);
	thread = std::thread([&](){ runLammps(); });
}
LammpsDriver::~LammpsDriver() {
	quit_thread = true;
	if (thread.joinable()) {
		thread.join();
	}
	lammps_close(lammps);
	MPI_Comm_free(&lammps_comm);
}
void LammpsDriver::runLammps() {
	std::cout << "Loading lammps input: '" << lammps_input << "'\n";
	// Read the input file and setup the problem by Bcasting out to other ranks
	if (rank == 0) {
		std::ifstream input(lammps_input.c_str());
		for (std::string line; std::getline(input, line);) {
			int len = line.size();
			// Skip empty lines
			if (len == 0) {
				continue;
			}
			MPI_Bcast(&len, 1, MPI_INT, 0, lammps_comm);
			MPI_Bcast(&line[0], len, MPI_CHAR, 0, lammps_comm);
			lammps_command(lammps, &line[0]);
		}
		// Bcast out we're done with the file
		int len = 0;
		MPI_Bcast(&len, 1, MPI_INT, 0, lammps_comm);
	} else {
		while (true) {
			int len = 0;
			MPI_Bcast(&len, 1, MPI_INT, 0, lammps_comm);
			if (len == 0) {
				break;
			} else {
				std::vector<char> line(len + 1, '\0');
				MPI_Bcast(line.data(), len, MPI_CHAR, 0, lammps_comm);
				lammps_command(lammps, line.data());
			}
		}
	}

	// Setup the fix external callback
	int ifix = lammps->modify->find_fix_by_style("external");
	FixExternal *fix = (FixExternal*)lammps->modify->fix[ifix];
	fix->set_callback(proxy_callback, this);

	libISInitWithExisting(lammps_comm, world_comm);
	while (!quit_thread) {
		lammps_command(lammps, "run 10");
	}
	libISFinalize();
}
void LammpsDriver::lammps_callback(int64_t ntimestep, int nlocal,
		int *id, double **x, double **f)
{
	libISBox3f bounds = libISMakeBox3f();
	const double *lmp_lo = static_cast<const double*>(lammps_extract_global(lammps, "boxlo"));
	const double *lmp_hi = static_cast<const double*>(lammps_extract_global(lammps, "boxhi"));
	libISVec3f islo{lmp_lo[0], lmp_lo[1], lmp_lo[2]};
	libISVec3f ishi{lmp_hi[0], lmp_hi[1], lmp_hi[2]};
	libISBoxExtend(&bounds, &islo);
	libISBoxExtend(&bounds, &ishi);

	const int *type = static_cast<const int*>(lammps_extract_atom(lammps, "type"));
	// TODO: Which side is better to convert on? We have
	// to convert here anyway since LAMMPS is SoA and we're
	// expecting AoS data coming in to libIS.
	struct LammpsAtom {
		float x, y, z;
		int type;
	};
	std::vector<LammpsAtom> atoms;
	atoms.reserve(nlocal);
	for (int i = 0; i < nlocal; ++i) {
		atoms.push_back(LammpsAtom{(*x)[i * 3], (*x)[i * 3 + 1],
				(*x)[i * 3 + 2], type[i]});
	}
	libISProcess(nlocal, reinterpret_cast<char*>(atoms.data()),
			sizeof(LammpsAtom), bounds); 
}

void proxy_callback(void *ptr, bigint ntimestep, int nlocal, int *id,
		double **x, double **f)
{
	LammpsDriver *driver = static_cast<LammpsDriver*>(ptr);
	driver->lammps_callback(ntimestep, nlocal, id, x, f);
}

