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

#include <set>
#include <vtkSmartPointer.h>
#include <vtkMolecule.h>
#include <vtkPoints.h>
#include <vtkIntArray.h>
#include <vtkDataSetAttributes.h>
#include <vtkSimpleBondPerceiver.h>
#include <mpiCommon/MPICommon.h>
#include "ospray/ospray_cpp/Data.h"
#include "ospray/ospray_cpp/Model.h"
#include "colormap.h"
#include "kd_tree.h"
#include "query_task.h"

using namespace ospcommon;
using namespace ospray::cpp;

extern float radius;
extern float bond_dist;
extern int bond_atom_type;

Sphere::Sphere() : pos(0), atom_type(0) {}
Sphere::Sphere(float x, float y, float z, int type)
	: pos(x, y, z), atom_type(type)
{}
Sphere::Sphere(vec3f v, int type)
	: pos(v), atom_type(type)
{}

SimulationState::~SimulationState() {
	for (auto &m : models) {
		m.release();
	}
	for (auto &gm : ghost_models) {
		gm.release();
	}
	for (auto &g : geom) {
		g.release();
	}
}

QueryTask::QueryTask(MPI_Comm worker_comm) : quitThread(false) {
	MPI_Comm_dup(worker_comm, &query_comm);
	thread = std::thread([&](){ queryThread(); });
}
QueryTask::~QueryTask() {
	quitThread = true;
	if (thread.joinable()) {
		thread.join();
	}
	MPI_Comm_free(&query_comm);
}
std::shared_ptr<SimulationState> QueryTask::take() {
	std::lock_guard<std::mutex> lock(mutex);
	if (available_state) {
		std::shared_ptr<SimulationState> ret = available_state;
		available_state = nullptr;
		return ret;
	}
	return nullptr;
}
void QueryTask::queryThread() {
	do {
		std::shared_ptr<SimulationState> state = std::make_shared<SimulationState>();
        bool simDisconnected = false;
		state->regions = is::client::query(&simDisconnected);
		if (quitThread) {
			break;
		}
        if (!is::client::sim_connected()) {
            std::cout << "Simulation disconnected, quitting\n";
            std::exit(0);
        }

		struct ISSphere { vec3f pos; int type; };
		// Transform the particles so we can give them radii based on the atom type
		for (auto &r : state->regions) {
			auto new_spheres = std::make_shared<is::OwnedArray>(sizeof(Sphere) * r.particles.array->size(),
					sizeof(Sphere));

			const ISSphere *isspheres = reinterpret_cast<const ISSphere*>(r.particles.array->data());
			Sphere *spheres = reinterpret_cast<Sphere*>(new_spheres->data());
			for (size_t i = 0; i < r.particles.array->size(); ++i) {
				spheres[i].pos = isspheres[i].pos;
				spheres[i].atom_type = isspheres[i].type;
				spheres[i].radius = isspheres[i].type == 1 ? 1.7 : 0.8;
			}
			r.particles.array = new_spheres;
		}

		// TODO: We probably don't need to recompute the colormap each time,
		// unless the sim is adding new atoms in.
		std::array<int, 2> atom_type_range;
		auto colormap = colormap_atoms(state->regions, atom_type_range);
		colormap[0] = vec4f(5.f, 113.f, 176.f, 255.f) / 255.f;
		colormap[1] = vec4f(202.f, 0.f, 32.f, 255.f) / 255.f;

		// TODO: We need to compute the world bounds ourselves, since the LAMMPS
		// driver doesn't send this right now
		box3f local_bounds;
		for (const auto &r : state->regions) {
			const vec3f region_lower = vec3f(r.local.min.x, r.local.min.y, r.local.min.z) - vec3f(radius); 
			const vec3f region_upper = vec3f(r.local.max.x, r.local.max.y, r.local.max.z) + vec3f(radius);

			local_bounds.extend(region_lower);
			local_bounds.extend(region_upper);

			Model model;
			Model ghost_model;
			model.set("id", r.simRank);
			ghost_model.set("id", r.simRank);
			model.set("region.lower", region_lower);
			model.set("region.upper", region_upper);

			auto atoms = make_atom_geometry(r, colormap);
			model.addGeometry(atoms);
			ghost_model.addGeometry(atoms);
			atoms.release();
			if (bond_atom_type > -1) {
				auto bonds = make_bond_geometry(r, colormap, state->world_bounds, atom_type_range);
				if (bonds) {
					model.addGeometry(bonds);
					ghost_model.addGeometry(bonds);
					bonds.release();
				}
			}
			model.commit();
			ghost_model.commit();
			state->models.push_back(model);
			state->ghost_models.push_back(ghost_model);
		}
		MPI_Allreduce(&local_bounds.lower, &state->world_bounds.lower, 3, MPI_FLOAT, MPI_MIN,
				query_comm);
		MPI_Allreduce(&local_bounds.upper, &state->world_bounds.upper, 3, MPI_FLOAT, MPI_MAX,
				query_comm);

		std::lock_guard<std::mutex> lock(mutex);
		available_state = state;
	} while (!quitThread);
}
std::vector<vec4f> QueryTask::colormap_atoms(std::vector<is::SimState> &regions,
		std::array<int, 2> &atom_type_range)
{
	atom_type_range[0] = std::numeric_limits<int>::max();
	atom_type_range[1] = std::numeric_limits<int>::min();

	for (const auto &region : regions) {
		const Sphere *atoms = reinterpret_cast<const Sphere*>(region.particles.array->data());
		// TODO: We need to do a reduction to compute the global min/mix ids and generate
		// colors for the unique atom types
		auto min_max_type = std::minmax_element(atoms, atoms + region.particles.numParticles,
			[](const Sphere &a, const Sphere &b) {
				return a.atom_type < b.atom_type;
			});
		atom_type_range[0] = std::min(min_max_type.first->atom_type, atom_type_range[0]);
		atom_type_range[1] = std::max(min_max_type.second->atom_type, atom_type_range[1]);
	}

	const int localMin = atom_type_range[0];
	const int localMax = atom_type_range[1];
	MPI_Allreduce(&localMin, &atom_type_range[0], 1, MPI_INT, MPI_MIN, query_comm);
	MPI_Allreduce(&localMax, &atom_type_range[1], 1, MPI_INT, MPI_MAX, query_comm);

	if (atom_type_range[0] != 0) {
		for (auto &region : regions) {
			Sphere *atoms = reinterpret_cast<Sphere*>(region.particles.array->data());
			const size_t n_atoms = region.particles.numParticles + region.particles.numGhost;
			for (size_t i = 0; i < n_atoms; ++i) {
				atoms[i].atom_type -= atom_type_range[0];
			}
		}
		atom_type_range[1] -= atom_type_range[0];
		atom_type_range[0] = 0;
	}

	std::vector<vec4f> atomColors;
	Colormap colormap;
	for (int i = atom_type_range[0]; i <= atom_type_range[1]; ++i) {
		const float x = (i - atom_type_range[0]) / static_cast<float>(atom_type_range[1] - atom_type_range[0]);
		atomColors.push_back(vec4f(colormap.map(x), 1));
	}
	return atomColors;
}
ospray::cpp::Geometry QueryTask::make_atom_geometry(const is::SimState &region,
		const std::vector<vec4f> &atom_colors)
{
	Data color_data = Data(atom_colors.size(), OSP_FLOAT4, atom_colors.data());
	color_data.commit();
	// The geometries are interleaved with the spheres representing the atoms
	// and the cylinders representing the bonds
	Geometry geom("spheres");

	// Setup the local sphere data for this region
	Data sphere_data = Data(region.particles.array->numBytes(),
			OSP_UCHAR, region.particles.array->data(), OSP_DATA_SHARED_BUFFER);
	sphere_data.commit();
	geom.set("spheres", sphere_data);
	geom.set("color", color_data);
	geom.set("offset_colorID", int(sizeof(vec3f)));
	geom.set("offset_radius", int(sizeof(vec3f) + sizeof(int)));
	geom.set("bytes_per_sphere", int(sizeof(Sphere)));
	geom.set("radius", radius);
	geom.commit();
	return geom;
}
// Map the simulation atom type to the periodic table atom type for the
// rhodopsin simulation
uint16_t map_rhodo_atom_types(int type) {
	switch (type) {
		case 1: case 2: case 3: case 4: case 5: case 6: case 7: case 8: case 9:
		case 36: case 37: case 38: case 39: case 40: case 54: case 55: case 56:
			return 1;

		case 10: case 11: case 12: case 13: case 14: case 15: case 16: case 17:
		case 18: case 19: case 20: case 21: case 22: case 41: case 42: case 43:
		case 44: case 45: case 46: case 58: case 59: case 60: case 61: case 62:
		case 63: case 64: case 65: case 66: case 67: case 68:
			return 6;

		case 23: case 24: case 25: case 26: case 27: case 28: case 29: case 47: case 57:
			return 7;

		case 30: case 31: case 32: case 33: case 48: case 49: case 50:
			return 8;

		case 34: case 35:
			return 16;

		case 51:
			return 15;

		case 52:
			return 16;

		case 53:
			return 11;
	}
	return -1;
}
uint16_t map_silicine_atom_types(int type) {
	switch (type) {
		case 0: return 77; // atom type 1 is Iridium - rescaled to 0
		case 1: return 14; // atom type 2 is Silicon - rescaled to 1
	}
	return -1;
}
ospray::cpp::Geometry QueryTask::make_bond_geometry(const is::SimState &region,
		const std::vector<vec4f> &atom_colors, const box3f &world_bounds,
		const std::array<int, 2> &atom_type_range)
{
	std::cout << "Computing bonds" << std::endl;
	const int bond_type = bond_atom_type;
	std::vector<vec3f> cylinder_pts;
	std::cout << "Computing bonds for region with "
		<< region.particles.numParticles << " atoms and "
		<< region.particles.numGhost << " ghost, total of "
		<< region.particles.numParticles + region.particles.numGhost
		<< " atoms\n";

	auto molecule = vtkSmartPointer<vtkMolecule>::New();
	const Sphere *spheres = reinterpret_cast<const Sphere*>(region.particles.array->data());
	size_t firstGhostAtom = 0;
	for (size_t i = 0; i < region.particles.numParticles + region.particles.numGhost; ++i) {
		if (spheres[i].atom_type == bond_atom_type) {
			if (i < region.particles.numParticles) {
				++firstGhostAtom;
			}
			const uint16_t type = map_silicine_atom_types(spheres[i].atom_type);
			molecule->AppendAtom(type, spheres[i].pos.x, spheres[i].pos.y, spheres[i].pos.z);
		}
	}

	auto bonds_molecule = vtkSmartPointer<vtkMolecule>::New();
	auto bond_filter = vtkSmartPointer<vtkSimpleBondPerceiver>::New();
	bond_filter->SetInputData(molecule.Get());
	bond_filter->SetOutput(bonds_molecule.Get());
	bond_filter->SetTolerance(bond_dist);
	bond_filter->Update();
	std::cout << "Computed " << bonds_molecule->GetNumberOfBonds() << " bonds" << std::endl;

	if (bonds_molecule->GetNumberOfBonds() > 0) {
		cylinder_pts.reserve(bonds_molecule->GetNumberOfBonds() * 2);
		for (size_t i = 0; i < bonds_molecule->GetNumberOfBonds(); ++i) {
			vtkBond bond = bonds_molecule->GetBond(i);
			if (bond.GetBeginAtomId() < firstGhostAtom || bond.GetEndAtomId() < firstGhostAtom) {
				vec3f p;
				bonds_molecule->GetAtomPosition(bond.GetBeginAtomId(), &p.x);
				cylinder_pts.push_back(p);

				bonds_molecule->GetAtomPosition(bond.GetEndAtomId(), &p.x);
				cylinder_pts.push_back(p);
			}
		}

		Data cylinder_data = Data(cylinder_pts.size() * sizeof(vec3f), OSP_UCHAR,
				cylinder_pts.data());

		const std::vector<vec4f> cylinder_colors(cylinder_pts.size() / 2, atom_colors[bond_type]);
		Data cylinder_color_data = Data(cylinder_colors.size(), OSP_FLOAT4,
				cylinder_colors.data());

		Geometry geom("cylinders");
		geom.set("cylinders", cylinder_data);
		geom.set("radius", 0.3);
		geom.set("color", cylinder_color_data);
		geom.commit();
		return geom;
	}
	return nullptr;
}

