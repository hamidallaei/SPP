#ifndef _VISUALPARTICLE_
#define _VISUALPARTICLE_

#include "particle.h"

const int circle_points_num = 100;

GLfloat unitcircle[2*circle_points_num];
GLuint indices[circle_points_num];


void Init_Circle()
{
	for (int i = 0; i < circle_points_num; i++)
	{
		unitcircle[2*i] = cos(2*i*PI/circle_points_num);
		unitcircle[2*i+1] = sin(2*i*PI/circle_points_num);
		indices[i] = i;
	}
}

struct RGB{
	float red,green,blue;
};

class VisualParticle : public BasicParticle {
public:
	static float radius;
	static float tail;
	static RGB color;
	void Draw();
};

float VisualParticle::radius;
float VisualParticle::tail;
RGB VisualParticle::color;

void VisualParticle::Draw()
{
	glEnableClientState (GL_VERTEX_ARRAY);
	glVertexPointer (2, GL_FLOAT, 0, unitcircle);

	glTranslatef(r.x,r.y,0);
	glScalef(radius,radius,1);

	glDrawElements(GL_POLYGON,circle_points_num,GL_UNSIGNED_INT,indices);

	float scale = 20*tail;
	glBegin(GL_LINES);
		glVertex3f(0, 0, 0.0);
		glVertex3f(-scale*v.x,-scale*v.y, 0);
	glEnd();
	
	glDisableClientState( GL_VERTEX_ARRAY );
	glLoadIdentity();
}

#endif
