#ifndef _SPHMATHHANDLER_HPP_
#define _SPHMATHHANDLER_HPP_
#include "SphParticle2d.hpp"
#include <math.h>
namespace sph
{
	class SphMathHandler
	{
	protected:
		int cntAllParticles;
		SphParticle2d *particles;
		double mass,viscosity,h,density0,k;
	    double PI,g,dt;
		int i,j;
		void initParameters(const int count)
		{
			cntAllParticles=count;
			PI=3.14159265;
			g = 9.81;
			dt=0.00001;
			mass=0.0033;
			viscosity=0.0001;
			h=3;
			density0=1000;
			k=100;
		}
		void createParticles()
		{	
			particles = new SphParticle2d[cntAllParticles];
			for(int l=0,l2=0, i=0;i<cntAllParticles;i++,l2+=5)
			{
					particles[i].cr=fRand(0,1);
					particles[i].cg=fRand(0,1);
					particles[i].cb=fRand(0,1);
					particles[i].x = 50+l2;
					if(l2>100)
					{
						l2=0;
						l+=5;
					}
					particles[i].y=50+l;
			}
		}
		double fRand(double fMin, double fMax)
		{
			double f = (double)rand() / RAND_MAX;
			return fMin + f * (fMax - fMin);
		}
		double getDistance(int _i)
		{
			return sqrt(particles[_i].x*particles[_i].x+particles[_i].y*particles[_i].y);
		}
		double getVelocity(int _i)
		{
			return sqrt(particles[_i].vx*particles[_i].vx+particles[_i].vy*particles[_i].vy);
		}
		double kernel(double d)
		{		
			
			if(abs(d)<=h)
			{
				
				return pow(h*h-d*d,3);
			}
			else
				return 0;
		}
		double gradkernel(double d)
		{	
		
			if(abs(d)<=h)
			{
			
				return pow(h-d,2);
			}
			else 
				return 0;
		}
		double dblgradkernel(double d)
		{
	
			if(abs(d)<=h)
			{
				
				return h-d;
			}
			else 
				return 0;
		}
		void solveSystem()
		{
			//coeffs for kernels
			double cpoly6=315/(64*PI*pow(h,9));
			double cpress=-45/(PI*pow(h,6));
			double cvisc=45/(PI*pow(h,6));
			//densities and pressure
			for(i=0;i<cntAllParticles;i++)
			{
				particles[i].density=0;
				for(j=0;j<cntAllParticles;j++)
				{
					double q=getDistance(i)-getDistance(j);
					particles[i].density += mass*kernel(q)*cpoly6;		
				}
				particles[i].pressure=k*(particles[i].density-density0);
			}
			
			//forces and leapfrog
			for(i=0;i<cntAllParticles;i++)
			{
				double fpress=0;
				for(j=0;j<cntAllParticles;j++)
				{
					if(particles[j].density==0)
					{
						continue;
					}
					double q = getDistance(i)-getDistance(j);
					fpress += -mass*((particles[j].pressure+particles[i].pressure)/(2*particles[j].density))*gradkernel(q)*cpress;
					//double fvisc = viscosity*mass*(getVelocity(j)-getVelocity(i))*dblgradkernel(q)*cvisc/(particles[i].density*particles[j].density);
					//acceleration (forces are already divided by deensity)
				}
				if(particles[i].density!=0)
				{
					particles[i].ax=fpress/particles[i].density;
					particles[i].ay=fpress/particles[i].density-g ;
				}
				else
				{
					particles[i].ax=0;
					particles[i].ay=-g;
				}
				//integration

				particles[i].vx = particles[i].vx + particles[i].ax*dt;
				particles[i].vy = particles[i].vy + particles[i].ay*dt;

				particles[i].x = particles[i].x + particles[i].vx*dt;
				particles[i].y = particles[i].y + particles[i].vy*dt;
				//box collision 
				if(particles[i].x<=1)
				{
					particles[i].vx=-particles[i].vx;
					particles[i].x+=3;
				}
				if(particles[i].x>=640)
				{
					particles[i].vx=-particles[i].vx;
					particles[i].x-=3;
				}
				if(particles[i].y<=1)
				{
					particles[i].vy=-particles[i].vy;
					particles[i].y+=3;
				}
				if(particles[i].y>=480)
				{
					particles[i].vy=-particles[i].vy;
					particles[i].y-=3;
				}
			}
		}
	};
}
#endif