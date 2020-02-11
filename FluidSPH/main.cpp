#include <iostream>
#include <string>
#include "glut.h"
#include "SPHSystem.h"
using namespace std;
int winX = 600;
int winY = 600;
SPHSystem* sph;
void init() {
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(0.0,
		sph->getWorldSize().x,
		0.0,
		sph->getWorldSize().y);
	glViewport(0, 0, winX, winY);
	glBlendFunc(GL_SRC_ALPHA,
		GL_ONE_MINUS_SRC_ALPHA);
	glEnable(GL_BLEND);
	glEnable(GL_POINT_SMOOTH);

	glHint(GL_POINT_SMOOTH_HINT,
			GL_NICEST);
	glClearColor(1.0f, 1.0f,
		1.0f, 1.0);
}
int divider = 100;
void drawParticles() {
	Particle* p = sph->getParticles();
	glColor3f(0.2f, 0.2f, 1.0f);
	glPointSize(2.0f);
	for (uint i = 0; i < sph->getNumParticle(); i++) {
		glBegin(GL_POINTS);
		glVertex2f(p[i].pos.x, p[i].pos.y);
		glEnd();
		float d = p[i].vel.Length();
		glBegin(GL_LINES);
		glVertex2f(p[i].pos.x, p[i].pos.y);
		if (d <= 50)
			glVertex2f(p[i].pos.x + p[i].vel.x / divider, p[i].pos.y +
				p[i].vel.y / divider);
		else
			glVertex2f(p[i].pos.x + p[i].vel.x / d, p[i].pos.y + p[i].vel.y / d);
		glEnd();
	}
}
void displayFunc() {
	sph->animation();
	glClear(GL_COLOR_BUFFER_BIT);
	drawParticles();
	glutSwapBuffers();
}
void idleFunc() {
	glutPostRedisplay();
}
void reshapeFunc(int width, int height) {
	winX = width;
	winY = height;
	glViewport(0, 0, winX, winY);
	glutReshapeWindow(winX, winY);
}
int main(int argc, char** argv) {
	sph = new SPHSystem();
	sph->initFluid();
	
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
