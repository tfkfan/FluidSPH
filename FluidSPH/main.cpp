#include <iostream>
#include <time.h>
#include <cstdio>
#include "glut.h"
#include "sph.h"

using namespace std;

int wndWidth = 700, wndHeight = 700;

bool paused = false;
bool track_gravity = true;
const double num_segments = 20;
const double r = 0.3;
const double WIDTH = 50;
const double HEIGHT = 50;
int sphereId;

Vector2d gravity_direction(0,-1);

const int particle_count = 3000;

SphFluidSolver solver(WIDTH, HEIGHT);
const float gravity = 100.0f;
const float scale = 1.0f;
float collision_restitution = 1.0f;

void init_liquid() {
	Particle *particles = new Particle[particle_count];
	int count = particle_count;
	Particle *particle_iter = particles;
	while (true) {
		for (int j = 0; j < HEIGHT; j++) {
			for (int i = 0; i < WIDTH; i++) {
				if (count-- == 0) {
					solver.init_particles(particles, particle_count);
					return;
				}
				particle_iter->position.x = i / scale;
				particle_iter->position.y = j / scale;
				particle_iter++;
			}
		}
	}
}
void draw_particle(Particle &particle) {
	glColor3f(0,0,0);
	glBegin(GL_TRIANGLE_FAN);
	for (int i = 0; i< num_segments; i++)   {
		double theta = 2.0f * 3.1415926f * float(i) / float(num_segments);//get the current angle 
		double cx = r * cosf(theta);//calculate the x component 
		double cy = r * sinf(theta);//calculate the y component 
		glVertex2d(particle.position.x + cx, particle.position.y + cy);//output vertex 
	}
	glEnd();
	/*
	Vector3f p = scale * particle.position;
	
	glTranslatef(+p.x, +p.y, +p.z);
	glCallList(sphereId);
	glTranslatef(-p.x, -p.y, -p.z);
	*/
}

void add_gravity_force(Particle &particle) {
	particle.force += gravity * gravity_direction * particle.density;
}

void add_global_forces() {
	solver.foreach_particle(add_gravity_force);
}

void handle_particle_collision_cube(Particle &particle) {
	double &px = particle.position.x;
	double &py = particle.position.y;
	double &vx = particle.velocity.x;
	double &vy = particle.velocity.y;

	if (px < 0 || px > WIDTH / scale) {
		px = min(max(px, 0.0), WIDTH / scale);
		vx *= -collision_restitution;
	}
	if (py < 0 || py > HEIGHT / scale) {
		py = min(max(py, 0.0), HEIGHT / scale);
		vy *= -collision_restitution;
	}
}

void handle_collisions() {
	solver.foreach_particle(handle_particle_collision_cube);
}


void reshape(int width, int height) {
	glViewport(0, 0, width, height);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(0,WIDTH,0,HEIGHT);
	glClearColor(1.0f,1.0f,1.0f,1.0f);
}

void keyboard(unsigned char key, int x, int y) {
	switch (key) {
	case ' ':
		paused = !paused;
		break;
	case 'g':
	case 'G':
		//solver.reset_particles();
		break;
	case 'q':
	case 'Q':
	case 0x1bU: /* ESC */
		exit(0);
	default:
		break;
	}
}

void display() {
	glClear(GL_COLOR_BUFFER_BIT);
	glClearColor(1.0f,1.0f,1.0f,1.0f);

	if (!paused) 
		solver.update(add_global_forces, handle_collisions);
	solver.foreach_particle(draw_particle);

	glutSwapBuffers();
	glPopMatrix();
}

void idle() {
	glutPostRedisplay();
}

void print_usage() {
	cout << endl;
	cout << "KEYSTROKE       ACTION" << endl;
	cout << "=========       ======" << endl << endl;
	cout << "q, Q, <ESC>     exit the program" << endl;
	cout << endl;
}

int main(int argc, char *argv[]) {
	print_usage();

	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE);
	glutInitWindowSize(wndWidth, wndHeight);
	glutInitWindowPosition(0, 0);
	glutCreateWindow("SPH Fluids");

	glutReshapeFunc(reshape);
	glutKeyboardFunc(keyboard);
	glutDisplayFunc(display);
	glutIdleFunc(idle);

	init_liquid();

	glutMainLoop();

	return (0);
}

