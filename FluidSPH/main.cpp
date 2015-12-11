#include <iostream>
#include <string>
#include "glut.h"
#include "SPHSystem.h"

using namespace std;

int winX = 600;
int winY = 600;

SPHSystem *sph;

void init(){
	glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
	gluOrtho2D(0.0, sph->getWorldSize().x, 0.0, sph->getWorldSize().y);
	glViewport(0, 0, winX, winY);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glEnable(GL_BLEND);
	glEnable(GL_POINT_SMOOTH);
	glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
	
	glClearColor(1.0f, 1.0f, 1.0f, 1.0);
}

void drawParticles(){
	Particle *p = sph->getParticles();
	glColor3f(0.2f, 0.2f, 1.0f);
	glPointSize(5.0f);

	glBegin(GL_POINTS);
		for(uint i=0; i<sph->getNumParticle(); i++)
		{
			glVertex2f(p[i].pos.x, p[i].pos.y);
		}
	glEnd();
}

void displayFunc(){
	sph->animation();

	glClear(GL_COLOR_BUFFER_BIT);
	drawParticles();
	glutSwapBuffers();
}

void idleFunc(){
	glutPostRedisplay();
}

void reshapeFunc(int width, int height){
	winX = width;
	winY = height;
	glViewport(0, 0, winX, winY);
	glutReshapeWindow(winX, winY);
}

int nParticles = 4000;

int main(int argc, char **argv){
	sph = new SPHSystem();
	sph->initFluid(nParticles);
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE);
	glutInitWindowSize(winX, winY);
	glutCreateWindow("SPH Fluid 2D");
	init();

	glutDisplayFunc(displayFunc);
	glutReshapeFunc(reshapeFunc);
	glutIdleFunc(idleFunc);
    glutMainLoop();

	free(sph);
	return 0;
}
