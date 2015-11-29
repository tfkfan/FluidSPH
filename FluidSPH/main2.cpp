#include <stdlib.h>
#include <iostream>
#include <time.h>
#include <math.h>
#include "glut.h"

using namespace std;


const int cntAllParticles=800;
const double dt=0.01;
const double mass = 1;
const double viscosity=3.1;
const double h = 0.17;
const double density0=1.2;
const double k=10;
const double PI=3.14159265;
const double g=-100;
const double gamma = 7;

const double cpoly6=4/(PI*pow(h,8));
const double cpress=-30/(PI*pow(h,5));
const double cvisc=-20/(PI*pow(h,5));
const int width=640,height=480;
const double xmin=0;
const double xmax=1;
const double ymin=0;
const double ymax=1;
class Particle2d
{
public:
	double x,y;
	double vx,vy;
	double ax,ay;
	double cr,cg,cb;
	double density;
	double fpx,fpy;
	double fv;
	Particle2d()
	{
		x=y=0;
		vy=vx=0;
		ax=ay=0;
		cr=cg=cb;
		density=fpx=fpy=fv=0;
	}
	~Particle2d()
	{
	}
};

Particle2d particles[cntAllParticles];

inline double fRand(double fMin, double fMax)
{
	double f = (double)rand() / RAND_MAX;
	return fMin + f * (fMax - fMin);
}
void initParticles()
{		
	double l=0,l2=0;
	for(int i=0;i<cntAllParticles;i++,l+=0.01)
	{
			particles[i].cr=fRand(0.5,1);
			particles[i].cg=fRand(0.5,1);
			particles[i].cb=fRand(0.5,1);
			
			particles[i].x = l;
			particles[i].y= l2;
			if(l>=0.2)
			{
				l=0;
				l2+=0.01;
			}
	}
	/*
	l=0;
	for(int i=cntAllParticles-300;i<cntAllParticles-200;i++,l+=0.01){
		particles[i].cr=1;
		particles[i].cg=1;
		particles[i].cb=1;

		particles[i].y = l;
		particles[i].x = 0;
	}
	l=0.01;
	for(int i=cntAllParticles-200;i<cntAllParticles-100;i++,l+=0.01){
		particles[i].cr=1;
		particles[i].cg=1;
		particles[i].cb=1;

		particles[i].x = l;
		particles[i].y = 0;
	}
	l=0.01;
	for(int i=cntAllParticles-100;i<cntAllParticles;i++,l+=0.01){
		particles[i].cr=1;
		particles[i].cg=1;
		particles[i].cb=1;

		particles[i].y = l;
		particles[i].x = 1;
	}
	*/
}
inline double getDistance(int i)
{
	return sqrt(particles[i].x*particles[i].x+particles[i].y*particles[i].y);
}
inline double getDistance(int i,int j)
{
	double dx = particles[j].x-particles[i].x;
	double dy = particles[j].y-particles[i].y;
	return sqrt(dx*dx+dy*dy);
}
inline double getVelocity(int i,int j)
{
	double dx = particles[j].vx-particles[i].vx;
	double dy = particles[j].vy-particles[i].vy;
	return sqrt(dx*dx+dy*dy);
}
inline double getVelocity(int i)
{
	return sqrt(particles[i].vx*particles[i].vx+particles[i].vy*particles[i].vy);
}
inline double kernel(double r){
  if (r < h)
    return cpoly6*pow(h*h-r*r,3);
  else
    return 0;
}
inline Particle2d  gradKernel(double r,double x,double y){
  if (r < 0.0001 && x!= 0 && y!=0)
  {
     Particle2d p;
     p.x=cpress*x*pow(h-r,2)/abs(x);
     p.y=cpress*y*pow(h-r,2)/abs(y);
     return p;
  }
  else{
    Particle2d p;
    return p;
  }
}
inline double doubleGradKernel(double r){
  if (r < h)
    return cvisc*(h-r);
  else
    return 0;
}
void solveSystem()
{
	for(int i=0;i<cntAllParticles;i++){
		Particle2d& pi = particles[i];
		pi.density = 0;
		pi.fpx = 0;
		pi.fpy = 0;
		pi.fv = 0;
	}

	for(int i=0;i<cntAllParticles;i++){
		Particle2d& pi = particles[i];
		for(int j=0;j<cntAllParticles;j++){
			Particle2d& pj = particles[j];
			double r = getDistance(i,j);
			double ker = kernel(r);
			pi.density += mass*ker;	
		    pj.density += mass*ker;
		}
	}

	for(int i=0;i<cntAllParticles;i++)
	{
		Particle2d& pi = particles[i];
		pi.ay=g;
		for(int j=0;j<cntAllParticles;j++)
		{
			Particle2d& pj = particles[j];
			double r = getDistance(i,j);

			Particle2d gradP = gradKernel(r,pj.x-pi.x,pj.y-pi.y);
			double temp = mass * k * ((pi.density - density0) + (pj.density - density0));
			if(pi.density!=0 && pj.density!=0){
				pi.fpx -= temp * gradP.x /pj.density;
				pi.fpy -= temp * gradP.y /pj.density;
			
				pj.fpx += temp * gradP.x /pi.density;
				pj.fpy += temp * gradP.y /pi.density;
			
				temp = mass * viscosity * getVelocity(i,j)* doubleGradKernel(r);
				pi.fv += temp / pj.density;
			
				pj.fv -= temp/ pi.density;
			}
		}
		if(pi.density!=0)
		{
			pi.ax+=(pi.fpx+pi.fv)/pi.density;
			pi.ay+=(pi.fpy+pi.fv)/pi.density;
		}
		else
		{
			pi.ax=0;
			pi.ay=g;
		}
		
		
		//integration
		pi.vx += pi.ax*dt;
		pi.vy += pi.ay*dt;

		pi.x += pi.vx*dt;
		pi.y += pi.vy*dt;

		//box collision 
		
		if(particles[i].x<=0)
		{
			particles[i].vx=-particles[i].vx/1.1;
			particles[i].x=0;
		}
		if(particles[i].x>=1)
		{
			particles[i].vx=-particles[i].vx/1.1;
			particles[i].x=1;
		}
		if(particles[i].y<=0)
		{
			particles[i].vy=-particles[i].vy/1.1;
			particles[i].y=0;
		}
		if(particles[i].y>=1)
		{
			particles[i].vy=-particles[i].vy/1.1;
			particles[i].y=1;
		}
		
		
	}
}
inline void DrawCircle(double x, double y, double r, int num_segments) 
{
	glBegin(GL_TRIANGLE_FAN);
	for (int i = 0; i< num_segments; i++)   {
		double theta = 2.0f * 3.1415926f * float(i) / float(num_segments);//get the current angle 
		double cx = r * cosf(theta);//calculate the x component 
		double cy = r * sinf(theta);//calculate the y component 
		glVertex2d(x + cx, y + cy);//output vertex 
	}
	glEnd();
}
void InitWindow(int width, int height)
{
	glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGB);
	glutInitWindowSize(width,height);
	glutInitWindowPosition(100,100);
	glutCreateWindow("TestSPH");
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(0,1,0,1);
	glClearColor(0.0f,0.0f,0.0f,1.0f);
}
void render()
{
	glClear(GL_COLOR_BUFFER_BIT);
	glColor3f(1,1,1);
	for(int i=0;i<cntAllParticles;i++)
	{
		Particle2d& p = particles[i];
		glColor3d(p.cr,p.cg,p.cb);
		DrawCircle(p.x,p.y,0.005,20);
	}
	/*
	glEnableClientState(GL_VERTEX_ARRAY);
	glVertexPointer(2,GL_DOUBLE,sizeof(Particle2d),particles);
	glDrawArrays(GL_POINTS,0, cntAllParticles);
	glDisableClientState(GL_VERTEX_ARRAY); 
	*/
	glutSwapBuffers();
}
void keyBD(unsigned char key,int x,int y)
{
	switch(key)
	{
		case  'w':
			exit(1);
			break;
	}
}
void timer(int=0)
{
	render();
	solveSystem();
	glutTimerFunc(1,timer,0);
}
int main(int argc, char** argv)
{
	 
	glutInitWindowSize(width, height);
    glutInit(&argc, argv);
    InitWindow(width,height);
	glutKeyboardFunc(keyBD);
	glutDisplayFunc(render);
	glutTimerFunc(1,timer,0);
	initParticles();
	glutMainLoop();
	return 0;
}