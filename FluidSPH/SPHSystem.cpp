#include "SPHSystem.h"
#include <cmath>
#include <string>
#include <iostream>

using namespace std;

SPHSystem::SPHSystem(){
	kernel_radius = 0.04f;
	mass = 0.02f;

	maxParticle = 10;
	numParticle = 0;

	worldSize.x = 2;
	worldSize.y = 2;
	cellSize = 1.00f;
	gridSize.x = (int)(worldSize.x/cellSize)+1;
	gridSize.y = (int)(worldSize.y/cellSize)+1;
	totCell = (uint)(gridSize.x) * (uint)(gridSize.y);
	
	//params
	gravity.x = 0.0f;
	gravity.y = -9.8f;
	stiffness = 1000.0f;
	restDensity = 1000.0f;
	timeStep = 0.0005f;
	wallDamping = 0.0f;
	viscosity = 8.1f;

	particles = new Particle[maxParticle];

	cells = new Cell*[gridSize.x];
	for(int i=0;i<gridSize.x;i++)
		cells[i] = new Cell[gridSize.y];

	cout << "SPHSystem:" << endl;
	cout << "GridSizeX:" << gridSize.x << endl;
	cout << "GridSizeY:" << gridSize.y << endl;
	//cout << "TotalCell:" << totCell << endl;
}

SPHSystem::~SPHSystem(){
	free(particles);
	free(cells);
}

void SPHSystem::initFluid(){
	
	float x = 0 , y = worldSize.y;
	float dx = kernel_radius*0.8,dy = kernel_radius*0.8;

	for(int i = 0; i < maxParticle; i++, x+=dx){
		if(x>=worldSize.x){
			x=0;
			y-=dy;
		}
		
		Vec2f pos;
		pos.x = x;
		pos.y = y;

		Particle p;
		p.pos.x = pos.x;
		p.pos.y = pos.y;

		particles[numParticle] = p;
		numParticle++;
	}
	cout << "NUM Particle:" <<  numParticle << endl;
}


Vec2i SPHSystem::calcCellPos(Vec2f pos){
	Vec2i res;
	res.x = (int)(pos.x/cellSize);
	res.y = (int)(pos.y/cellSize);
	return res;
}
	

void SPHSystem::compTimeStep(){
	Particle *p;
	float maxAcc = 0.0f;
	float curAcc;
	for(uint i=0; i<numParticle; i++){
		p = &(particles[i]);
		curAcc = p->acc.Length();
		if(curAcc > maxAcc) maxAcc=curAcc;
	}

	if(maxAcc > 0.0f)
		timeStep = kernel_radius/maxAcc*0.4f;
	else
		timeStep = 0.002f;
}

void SPHSystem::buildGrid(){
	for(int i=0;i<gridSize.x;i++)
		for(int j=0;j<gridSize.y;j++)
			cells[i][j].pList.clear();
	
	for(int i=0;i<maxParticle;i++){
			Vec2i index = calcCellPos(particles[i].pos);
			
			Cell& c = cells[index.x][index.y];
			list<Particle>& l = c.pList;
			l.push_back(particles[i]);
		}
}
void SPHSystem::compNearDensPressure(Particle& p, Vec2i cellPos){
	Vec2i nearPos;
	p.dens = 0.0f;
	p.pres = 0.0f;

	for(int m = -1; m <= 1; m++)
		for(int n = -1; n <= 1; n++){
			nearPos.x = cellPos.x + m;
			nearPos.y = cellPos.y + n;

			if(nearPos.x<0||nearPos.x>=gridSize.x||nearPos.y<0||nearPos.y>=gridSize.y)
				continue;

			list<Particle>& np = cells[nearPos.x][nearPos.y].pList;

			if(np.empty())
				continue;

			for (list<Particle>::iterator it = np.begin(); it != np.end(); it++){
				Vec2f distVec = it->pos - p.pos;
				float dist = distVec.LengthSquared();
					
				//Сомнительный код
				if(dist>=kernel_radius*kernel_radius)
					continue;

				p.dens = p.dens + mass * poly6(dist);
			}	
		}
	p.pres = (pow(p.dens / restDensity, 7) - 1) * stiffness;
}
void SPHSystem::compDensPressure(){
	Vec2i cellPos;

	for(int i = 0 ; i < gridSize.x; i++)
		for(int j = 0;j < gridSize.y; j++){

			Cell cell = cells[i][j];

			cellPos.x = i;
			cellPos.y = j;

			for (list<Particle>::iterator p_it = cell.pList.begin(); p_it != cell.pList.end(); p_it++)
				compNearDensPressure(*p_it, cellPos);
		}
}

void SPHSystem::compForce(){
	Vec2i cellPos;
	Vec2i nearPos;
	for(int i = 0 ; i < gridSize.x; i++)
		for(int j = 0;j < gridSize.y; j++){

			Cell cell = cells[i][j];

			cellPos.x = i;
			cellPos.y = j;

			for (list<Particle>::iterator p_it = cell.pList.begin(); p_it != cell.pList.end(); p_it++){
				Particle& p = *p_it;

				for(int m = -1; m <= 1; m++)
					for(int n = -1; n <= 1; n++){
						nearPos.x = cellPos.x + n;
						nearPos.y = cellPos.y + m;

						if(nearPos.x<0||nearPos.x>=gridSize.x||nearPos.y<0||nearPos.y>=gridSize.y)
							continue;

						list<Particle>& np = cells[nearPos.x][nearPos.y].pList;
						if(!np.empty()){
							for (list<Particle>::iterator it = np.begin(); it != np.end(); it++){
								Vec2f distVec = p.pos - it->pos;
								float dist2 = distVec.LengthSquared();

								if(dist2 < kernel_radius*kernel_radius && dist2 > INF){
									float dist = sqrt(dist2);
									float V = mass / it->dens;

									float tempForce = V * (p.pres + it->pres) * spiky(dist);
									p.acc = p.acc - distVec*tempForce/dist;

									Vec2f relVel;
									relVel = it->ev-p.ev;
									tempForce = V * viscosity * visco(dist);
									p.acc = p.acc + relVel*tempForce; 
								}
							}
						}
					}
			}
		}
}

void SPHSystem::advection(){
	for(int i = 0 ; i < gridSize.x; i++)
		for(int j = 0;j < gridSize.y; j++){
			list<Particle>& l = cells[i][j].pList;
			for (list<Particle>::iterator it = l.begin(); it != l.end(); it++){
				Particle& p = *it;
				p.vel = p.vel + p.acc*timeStep;
				p.pos = p.pos + p.vel*timeStep;

				if(p.pos.x < 0.0f){
					p.vel.x = p.vel.x * wallDamping;
					p.pos.x = 0.0f;
				}
				if(p.pos.x >= worldSize.x){
					p.vel.x = p.vel.x * wallDamping;
					p.pos.x = worldSize.x - 0.0001f;
				}
				if(p.pos.y < 0.0f){
					p.vel.y = p.vel.y * wallDamping;
					p.pos.y = 0.0f;
				}
				if(p.pos.y >= worldSize.y){
					p.vel.y = p.vel.y * wallDamping;
					p.pos.y = worldSize.y - 0.0001f;
				}

				p.ev=(p.ev+p.vel)/2;
			}
	}
}

void SPHSystem::animation(){
	buildGrid();
	compDensPressure();
	compForce();
	//compTimeStep();
	advection();
}
