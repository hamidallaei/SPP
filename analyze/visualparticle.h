#ifndef _VISUALPARTICLE_
#define _VISUALPARTICLE_

#include "../shared/particle.h"

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
	static float thickness;
	static RGB color;
	void Find_Color(RGB&);
	void Draw(float, float, float);
	void Draw_Magnified();
	void Draw_Magnified(C2DVector r0, float d0, C2DVector r1, float d1);
};

float VisualParticle::radius;
float VisualParticle::tail;
float VisualParticle::thickness;
RGB VisualParticle::color;

void HSV_To_RGB(float h, float s, float v,RGB& color)
{
	int index;
	float f, p, q, t;
	if( s == 0 ) {
		// achromatic (grey)
		color.red = color.green = color.blue = v;
		return;
	}

	h /= 2*M_PI;			// sector 0 to 5
	h *= 6;
	index = floor( h );
	f = h - index;			// factorial part of h
	p = v * ( 1 - s );
	q = v * ( 1 - s * f );
	t = v * ( 1 - s * ( 1 - f ) );

	switch( index ) {
		case 0:
			color.red = v;
			color.green = t;
			color.blue = p;
			break;
		case 1:
			color.red = q;
			color.green = v;
			color.blue = p;
			break;
		case 2:
			color.red = p;
			color.green = v;
			color.blue = t;
			break;
		case 3:
			color.red = p;
			color.green = q;
			color.blue = v;
			break;
		case 4:
			color.red = t;
			color.green = p;
			color.blue = v;
			break;
		default:		// case 5:
			color.red = v;
			color.green = p;
			color.blue = q;
			break;
	}
}

void VisualParticle::Find_Color(RGB& input_color)
{
	float f, p, q, t, s = 1, value = 1;



	float theta = atan2(v.y,v.x);
	if (theta < 0)
		theta += 2*M_PI;

	HSV_To_RGB(theta,1,1, input_color);
}

void VisualParticle::Draw(float x0 = 0, float y0 = 0, float scale = 1)
{
	glEnableClientState (GL_VERTEX_ARRAY);

	Find_Color(color);
	glColor4f(color.red, color.green, color.blue,0.5);
	glLineWidth(thickness);
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

	glVertexPointer (2, GL_FLOAT, 0, unitcircle);

	glTranslatef(r.x+x0,r.y+y0,0);
	glScalef(scale*radius,scale*radius,1);

	glDrawElements(GL_POLYGON,circle_points_num,GL_UNSIGNED_INT,indices);

	glBegin(GL_LINES);
		glVertex2f(0, 0);
		glVertex2f(-tail*v.x,-tail*v.y);
	glEnd();
	
	glDisableClientState( GL_VERTEX_ARRAY );
	glLoadIdentity();
	glLineWidth(1);
}

void VisualParticle::Draw_Magnified(C2DVector r0, float d0, C2DVector r1, float d1)
{
	float scale = d1 / d0;
	C2DVector p;
	p.x = r.x - r0.x;
	p.y = r.y - r0.y;
	if (fabs(p.x) < (d0-radius) && fabs(p.y) < (d0-radius))
	{
		p *= scale;
		p += r1;

		const int thickness_factor = 3; // 2000
		glLineWidth(scale*thickness / thickness_factor);
		glEnableClientState (GL_VERTEX_ARRAY);
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

		glVertexPointer (2, GL_FLOAT, 0, unitcircle);

		glTranslatef(p.x,p.y,0);
		glScalef(scale*radius,scale*radius,1);

		glDrawElements(GL_POLYGON,circle_points_num,GL_UNSIGNED_INT,indices);

		glBegin(GL_LINES);
			glVertex2f(0, 0);
			glVertex2f(-tail*v.x,-tail*v.y);
		glEnd();
	
		glDisableClientState( GL_VERTEX_ARRAY );
		glLoadIdentity();
		glLineWidth(1);
	}
}

#endif
