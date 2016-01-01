#ifndef __SPHSYSTEM_H__
#define __SPHSYSTEM_H__

#include "Vector2D.h"
#include<list>

#define PI 3.141592f
#define INF 1E-12f

using namespace std;

class Particle{
public:
	Vec2f pos;
	Vec2f vel;
	Vec2f acc;
	Vec2f ev;

	float dens;
	float pres;
	Particle(){	}
	Particle(const Particle& p){
		this->acc = p.acc;
		this->pos = p.pos;
		this->vel = p.vel;
		this->ev = p.ev;
	}
	~Particle(){}
};

class Cell{
public:
	list<Particle> *pList;
};

class SPHSystem{
public:
	SPHSystem();
	~SPHSystem();
	void initFluid();
	Vec2i calcCellPos(Vec2f pos);

	//kernel function
	float poly6(float r2){ return 315.0f/(64.0f * PI * pow(kernel_radius, 9)) * pow(kernel_radius*kernel_radius-r2, 3); }
	float spiky(float r){ return -45.0f/(PI * pow(kernel_radius, 6)) * (kernel_radius-r) * (kernel_radius-r); }
	float visco(float r){ return 45.0f/(PI * pow(kernel_radius, 6)) * (kernel_radius-r); }

	//animation
	void compTimeStep();
	void buildGrid();
	void compDensPressure();
	void compForce();
	void advection();
	void animation();

	//getters
	uint getNumParticle(){ return numParticle; }
	Vec2f getWorldSize(){ return worldSize; }
	Particle* getParticles(){ return particles; }
	Cell** getCells(){ return cells; }
	
	Vec2i gridSize;

private:
	float kernel_radius;
	float mass;

	uint maxParticle;
	uint numParticle;

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

	Particle *particles;
	Cell **cells;

	void compNearDensPressure(Particle& p, Vec2i cellPos);
	void compNearForce(Particle& p,Vec2i cellPos);
};

#endif
