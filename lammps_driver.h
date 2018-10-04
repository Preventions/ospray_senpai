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
#include <string>
#include <mpi.h>

#include "lammps.h"
#include "fix.h"

class LammpsDriver {
	int rank;
	std::string lammps_input;
	MPI_Comm world_comm, lammps_comm;
	LAMMPS_NS::LAMMPS *lammps;
	std::atomic<bool> quit_thread;
	std::thread thread;

public:
	LammpsDriver(const std::string &lammps_input, MPI_Comm world);
	~LammpsDriver();

private:
	void runLammps();
	void lammps_callback(int64_t ntimestep, int nlocal,
			int *id, double **x, double **f);
	friend void proxy_callback(void *ptr, LAMMPS_NS::bigint ntimestep,
			int nlocal, int *id, double **x, double **f);
};

