#include <list>

using namespace std;

class SphFluidSolver;
struct Particle;
struct GridElement;

#ifndef _SPH_H_
#define _SPH_H_
#include "Vector.h"
struct Particle {
	double density;
	Vector2d position;
	Vector2d velocity;
	Vector2d force;
	double red, green, blue;
	Vector2d viscosity_force;
	Vector2d pressure_force;
	Particle() { }
};

struct GridElement {
	list<Particle> particles;
};

const double gas_constant = 1000.0;
const double mu =  0.001; 
const double rest_density = 1.3; 
const double point_damping = 3;
const double core_radius = 1.5;
const double timestep = 0.01;
const double mass = 1.0;

class SphFluidSolver {
	const int grid_width;
	const int grid_height;
	GridElement *grid_elements;
	GridElement *sleeping_grid_elements;
public:
	SphFluidSolver(double domain_width, double domain_height) 
		: grid_width((int) (domain_width / core_radius) + 1), grid_height((int) (domain_height / core_radius) + 1) {}

	void update(void(*inter_hook)() = NULL, void(*post_hook)() = NULL);
	void init_particles(Particle *particles, int count);
	template <typename Function>
	void foreach_particle(Function function) {
		for (int j = 0; j < grid_height; j++) {
			for (int i = 0; i < grid_width; i++) {
				GridElement &grid_element = grid_elements[grid_height*j + i];
				list<Particle> &plist = grid_element.particles;
				for (list<Particle>::iterator piter = plist.begin(); piter != plist.end(); piter++) {
					function(*piter);
				}
			}
		}
	}

private:

	double kernel(const Vector2d &r);

	Vector2d gradient_kernel(const Vector2d &r);

	double laplacian_kernel(const Vector2d &r);

	Vector2d gradient_pressure_kernel(const Vector2d &r);

	double laplacian_viscosity_kernel(const Vector2d &r);

	void add_density(Particle &particle, Particle &neighbour);

	void sum_density(GridElement &grid_element, Particle &particle);

	void sum_all_density(int i, int j, Particle &particle);

	void update_densities(int i, int j);

	void add_forces(Particle &particle, Particle &neighbour);

	void sum_forces(GridElement &grid_element, Particle &particle);

	void sum_all_forces(int i, int j, Particle &particle);

	void update_forces(int i, int j);

	void update_particle(Particle &particle);

	void update_particles(int i, int j);

	void reset_particle(Particle &particle);

	void reset_particles();

	void insert_into_grid(int i, int j);

	void update_grid();

	void update_densities();

	void update_forces();

	void update_particles();

	GridElement &grid(int i, int j);

	GridElement &sleeping_grid(int i, int j);

	int grid_index(int i, int j);

	void add_to_grid(GridElement *target_grid, Particle &particle);
};

#endif /* _SPH_H_ */

