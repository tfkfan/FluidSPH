#include "sph.h"
#include <time.h>
#include <cstdio>

#include <iostream>

using namespace std;

#define PI_FLOAT				3.141592653589793f
#define SQR(x)					((x) * (x))
#define CUBE(x)					((x) * (x) * (x))
#define POW6(x)					(CUBE(x) * CUBE(x))
#define POW9(x)					(POW6(x) * CUBE(x))
//1)  315.0f / (64.0f * PI_FLOAT * POW9(h)) * CUBE(SQR(h) - dot(r, r))
//2)  Vector3f common = 0.5f * material.gas_constant
//			* ((particle.density - material.rest_density) + (neighbour.density - material.rest_density)) 
//if (dot(r, r) < SQR(0.001f)) {
//		return Vector3f(0.0f);
//	}
//	return -45.0f / (PI_FLOAT * POW6(h)) * SQR(h - length(r)) * normalize(r);
//3) return   945.0f / (32.0f * PI_FLOAT * POW9(h)) * (SQR(h) - dot(r, r)) * (7.0f * dot(r, r) - 3.0f * SQR(h));
inline double SphFluidSolver::kernel(const Vector2d &r) {
	return 315.0f / (64.0f * PI_FLOAT * POW9(core_radius)) * CUBE(SQR(core_radius) - dot(r, r));
}
inline Vector2d SphFluidSolver::gradient_kernel(const Vector2d &r) {
	return -945.0f / (32.0f * PI_FLOAT * POW9(core_radius)) * SQR(SQR(core_radius) - dot(r, r)) * r;
}
inline double SphFluidSolver::laplacian_kernel(const Vector2d &r) {
	return   945.0f / (32.0f * PI_FLOAT * POW9(core_radius))
	       * (SQR(core_radius) - dot(r, r)) * (7.0f * dot(r, r) - 3.0f * SQR(core_radius));
}
inline Vector2d SphFluidSolver::gradient_pressure_kernel(const Vector2d &r) {
	if (dot(r, r) < SQR(0.001f)) 
		return Vector2d(0.0f);
	return -45.0f / (PI_FLOAT * POW6(core_radius)) * SQR(core_radius - length(r)) * normalize(r);
}
inline double SphFluidSolver::laplacian_viscosity_kernel(const Vector2d &r) {
	return 45.0f / (PI_FLOAT * POW6(core_radius)) * (core_radius - length(r));
}
inline void SphFluidSolver::add_density(Particle &particle, Particle &neighbour) {
	if (particle.id > neighbour.id) 
		return;
	Vector2d r = particle.position - neighbour.position;
	if (dot(r, r) > SQR(core_radius)) 
		return;
    float common = kernel(r);
    particle.density += neighbour.mass * common;
	neighbour.density += particle.mass * common;
}
void SphFluidSolver::sum_density(GridElement &grid_element, Particle &particle) {
	list<Particle> &plist = grid_element.particles;
	for (list<Particle>::iterator piter = plist.begin(); piter != plist.end(); piter++) 
		add_density(particle, *piter);
}
inline void SphFluidSolver::sum_all_density(int i, int j, Particle &particle) {
	for (int y = j - 1; y <= j + 1; y++) {
		for (int x = i - 1; x <= i + 1; x++) {
			if (   (x < 0) || (x >= grid_width)
				|| (y < 0) || (y >= grid_height)) {
				continue;
			}

			sum_density(grid(x, y), particle);
		}
	}
}
void SphFluidSolver::update_densities(int i, int j) {
	GridElement &grid_element = grid(i, j);
	list<Particle> &plist = grid_element.particles;
	for (list<Particle>::iterator piter = plist.begin(); piter != plist.end(); piter++) 
		sum_all_density(i, j, *piter);
}
inline void SphFluidSolver::add_forces(Particle &particle, Particle &neighbour) {
	if (particle.id >= neighbour.id) 
		return;
	Vector2d r = particle.position - neighbour.position;
	if (dot(r, r) > SQR(core_radius)) 
		return;
	/* Compute the pressure force. */
	Vector2d common = 0.5f * gas_constant
			* ((particle.density - rest_density) + (neighbour.density - rest_density))
	        * gradient_pressure_kernel(r);
	particle.force += -neighbour.mass / neighbour.density * common;
	particle.pressure_force += -neighbour.mass / neighbour.density * common;
	neighbour.force -= -particle.mass / particle.density * common;
	neighbour.pressure_force -= -particle.mass / particle.density * common;
	/* Compute the viscosity force. */
	common = mu * (neighbour.velocity - particle.velocity)
	         * laplacian_viscosity_kernel(r);
	particle.force += neighbour.mass / neighbour.density * common;
	particle.viscosity_force += neighbour.mass / neighbour.density * common;
	neighbour.force -= particle.mass / particle.density * common;
	neighbour.viscosity_force -= particle.mass / particle.density * common;
	/* Compute the gradient of the color field. */
	common = gradient_kernel(r);
	particle.color_gradient += neighbour.mass / neighbour.density * common;
	neighbour.color_gradient -= particle.mass / particle.density * common;
	/* Compute the laplacian of the color field. */
	float value = laplacian_kernel(r);
	particle.color_laplacian += neighbour.mass / neighbour.density * value;
	neighbour.color_laplacian += particle.mass / particle.density * value;
}
void SphFluidSolver::sum_forces(GridElement &grid_element, Particle &particle) {
	list<Particle>  &plist = grid_element.particles;
	for (list<Particle>::iterator piter = plist.begin(); piter != plist.end(); piter++) 
		add_forces(particle, *piter);
}
void SphFluidSolver::sum_all_forces(int i, int j, Particle &particle) {
	for (int y = j - 1; y <= j + 1; y++) {
		for (int x = i - 1; x <= i + 1; x++) {
			if ((x < 0) || (x >= grid_width) || (y < 0) || (y >= grid_height)) 
				continue;
			sum_forces(grid(x, y), particle);
		}
	}
}
void SphFluidSolver::update_forces(int i, int j) {
	GridElement &grid_element = grid(i, j);
	list<Particle>&plist = grid_element.particles;
	for (list<Particle>::iterator piter = plist.begin(); piter != plist.end(); piter++) 
		sum_all_forces(i, j, *piter);
}
inline void SphFluidSolver::update_particle(Particle &particle) {
	if (length(particle.color_gradient) > 0.001f) 
		particle.force +=-sigma * particle.color_laplacian*normalize(particle.color_gradient);

	Vector2d acceleration =   particle.force / particle.density - point_damping * particle.velocity / particle.mass;
	particle.velocity += timestep * acceleration;
	particle.position += timestep * particle.velocity;
}
void SphFluidSolver::update_particles(int i, int j) {
	GridElement &grid_element = grid(i, j);
	list<Particle> &plist = grid_element.particles;
	for (list<Particle>::iterator piter = plist.begin(); piter != plist.end(); piter++) 
		update_particle(*piter);
}
inline void SphFluidSolver::reset_particle(Particle &particle) {
	particle.density = 0.0;
	particle.force = Vector2d(0.0);
	particle.viscosity_force = Vector2d(0.0);
	particle.pressure_force = Vector2d(0.0);
	particle.color_gradient = Vector2d(0.0);
	particle.color_laplacian = 0.0;
}
void SphFluidSolver::reset_particles() {
	for (int j = 0; j < grid_height; j++) {
		for (int i = 0; i < grid_width; i++) {
			GridElement &grid_element = grid(i, j);
			list<Particle> &plist = grid_element.particles;
			for (list<Particle>::iterator piter = plist.begin(); piter != plist.end(); piter++) 
				reset_particle(*piter);
		}
	}
}
inline void SphFluidSolver::insert_into_grid(int i, int j) {
	GridElement &grid_element = grid(i, j);
	list<Particle> &plist = grid_element.particles;
	for (list<Particle>::iterator piter = plist.begin(); piter != plist.end(); piter++) 
		add_to_grid(sleeping_grid_elements, *piter);
}
void SphFluidSolver::update_grid() {
	for (int j = 0; j < grid_height; j++) {
		for (int i = 0; i < grid_width; i++) {
			insert_into_grid(i, j);
			grid(i, j).particles.clear();
		}
	}
	swap(grid_elements, sleeping_grid_elements);
}
void SphFluidSolver::update_densities() {
	#pragma omp parallel for
	for (int j = 0; j < grid_height; j++) 
		for (int i = 0; i < grid_width; i++) 
			update_densities(i, j);
}
void SphFluidSolver::update_forces() {
	#pragma omp parallel for
	for (int j = 0; j < grid_height; j++) 
		for (int i = 0; i < grid_width; i++) 
			update_forces(i, j);
}
void SphFluidSolver::update_particles() {
	#pragma omp parallel for
	for (int j = 0; j < grid_height; j++) 
		for (int i = 0; i < grid_width; i++) 
			update_particles(i, j);
}
void SphFluidSolver::update(void(*inter_hook)(), void(*post_hook)()) {
	reset_particles();
    update_densities();
    update_forces();
	if (inter_hook != NULL) 
		inter_hook();
    update_particles();
	if (post_hook != NULL) 
		post_hook();
	update_grid();
}
void SphFluidSolver::init_particles(Particle *particles, int count) {
	grid_elements = new GridElement[grid_width * grid_height];
	sleeping_grid_elements = new GridElement[grid_width * grid_height];
	for (int x = 0; x < count; x++) {
		particles[x].id = x;
		add_to_grid(grid_elements, particles[x]);
	}
}
inline GridElement &SphFluidSolver::grid(int i, int j) {
	return grid_elements[grid_index(i, j)];
}
inline GridElement &SphFluidSolver::sleeping_grid(int i, int j) {
	return sleeping_grid_elements[grid_index(i, j)];
}
inline int SphFluidSolver::grid_index(int i, int j) {
	return grid_height*j + i;
}
inline void SphFluidSolver::add_to_grid(GridElement *target_grid, Particle &particle) {
	int i = (int) (particle.position.x / core_radius);
	int j = (int) (particle.position.y / core_radius);
	target_grid[grid_index(i, j)].particles.push_back(particle);
}

