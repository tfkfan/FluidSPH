#ifndef __SPHSYSTEM_H__
#define __SPHSYSTEM_H__
#include "Vector2D.h"
#define PI 3.141592f
#define INF 1E-12f
#define sqr(x) x*x
class Particle
{
public:
	Vec2f pos;
	Vec2f vel;
	Vec2f acc;
	Vec2f ev;
	float dens;
	float pres;
	Particle* next;
};
class Cell
{
public:
	Particle* head;
};
class SPHSystem
{
public:
	SPHSystem();
	~SPHSystem();
	void initFluid();
	void addSingleParticle(Vec2f pos, Vec2f vel);
	Vec2i calcCellPos(Vec2f pos);
	uint calcCellHash(Vec2i pos);
	//kernel function
	float poly6(float r2) {
		return 315.0f / (64.0f * PI * pow(kernel_radius, 9)) *
			pow(kernel_radius * kernel_radius - r2, 3);
	}
	float spiky(float r) {
		return -45.0f / (PI * pow(kernel_radius, 6)) * (kernel_radius - r) *
			(kernel_radius - r);
	}
	float visco(float r) {
		return 45.0f / (PI * pow(kernel_radius, 6)) * (kernel_radius - r);
	}
	//animation
	void compTimeStep();
	void buildGrid();
	void compDensPressure();
	void compForce();
	void advection();
	void animation();
	//getters
	uint getNumParticle() { return numParticle; }
	Vec2f getWorldSize() { return worldSize; }
	Particle* getParticles() { return particles; }
	Cell* getCells() { return cells; }
private:
	float kernel_radius;
	float mass;
	uint maxParticle;
	uint numParticle;
	Vec2i gridSize;
	Vec2f worldSize;
	float cellSize;
	uint totCell;
	//params
	Vec2f gravity;
	float stiffness;
	float restDensity;
	float timeStep;
	float wallDamping;
	float viscosity;
	Particle* particles;
	Cell* cells;
	float dx, dy;
	float ddx, ddy;
	Vec2f tmpPos;
};
#endif