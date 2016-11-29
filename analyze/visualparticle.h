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

void Draw_Circle(C2DVector r, C2DVector tail, float radius, float thickness, RGB color, GLenum mode)
{
	glEnableClientState (GL_VERTEX_ARRAY);

	glColor4f(color.red, color.green, color.blue,0.5);
	glLineWidth(thickness);

	glPolygonMode(GL_FRONT_AND_BACK, mode);
	glPolygonMode(GL_FRONT_AND_BACK, mode);

	glVertexPointer (2, GL_FLOAT, 0, unitcircle);

	glTranslatef(r.x,r.y,0);
	glScalef(radius,radius,1);

	glDrawElements(GL_POLYGON,circle_points_num,GL_UNSIGNED_INT,indices);

	glBegin(GL_LINES);
		glVertex2f(0, 0);
		glVertex2f(tail.x,tail.y);
	glEnd();
		glColor4f(color.red, color.green, color.blue,1.0);
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		glDrawElements(GL_POLYGON,circle_points_num,GL_UNSIGNED_INT,indices);
	
	glDisableClientState( GL_VERTEX_ARRAY );
	glLoadIdentity();
	glLineWidth(1);
}

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

class VisualMembrane : public BasicParticle00 {
public:
	static float radius;
	static float thickness;
	static GLenum mode; // GL_FILL (filled) or GL_LINE (empty) or GL_POINT

	virtual void Draw(float, float, float);
	virtual void Draw_Magnified(C2DVector r0, float d0, C2DVector r1, float d1);
};

float VisualMembrane::radius = 0.5;
float VisualMembrane::thickness = 1;
GLenum VisualMembrane::mode = GL_FILL;

void VisualMembrane::Draw(float x0 = 0, float y0 = 0, float scale = 1)
{
	C2DVector null(0);
	null.Null();

	RGB color;
	color.red = 0.2;
	color.green = 0.2;
	color.blue = 0.2;

	Draw_Circle(r, null, radius*scale, thickness, color, mode);
}

void VisualMembrane::Draw_Magnified(C2DVector r0, float d0, C2DVector r1, float d1)
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

		C2DVector null(0);
		null.Null();

		RGB color;
		color.red = 0.2;
		color.green = 0.2;
		color.blue = 0.2;

		Draw_Circle(p, null, radius*scale, scale*thickness / thickness_factor, color, mode);
	}
}


class VisualParticle : public BasicParticle0 {
public:
	C2DVector v;
	static float radius;
	static float tail;
	static float thickness;
	static GLenum mode; // GL_FILL (filled) or GL_LINE (empty) or GL_POINT

	void Find_Color(RGB&);
	virtual void Draw(float, float, float);
	virtual void Draw_Magnified(C2DVector r0, float d0, C2DVector r1, float d1);
};

float VisualParticle::radius;
float VisualParticle::tail;
float VisualParticle::thickness;
GLenum VisualParticle::mode = GL_FILL;

void VisualParticle::Find_Color(RGB& input_color)
{
	float f, p, q, t, s = 1, value = 1;

	theta = theta - floor(theta/(2*M_PI))*2*M_PI;
	HSV_To_RGB(theta,1,1, input_color);
}

void VisualParticle::Draw(float x0 = 0, float y0 = 0, float scale = 1)
{
	RGB color;
	Find_Color(color);
	Draw_Circle(r, v*(-tail), radius*scale, thickness, color, mode);
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

		RGB color;
		Find_Color(color);
		Draw_Circle(p, v*(-tail), radius*scale, scale*thickness / thickness_factor, color, mode);
	}
}


class VisualChain : public VisualParticle {
public:
	static int chain_length;

	void Draw(float, float, float);
	void Draw_Magnified(C2DVector r0, float d0, C2DVector r1, float d1);
};

int VisualChain::chain_length = 2;

void VisualChain::Draw(float x0 = 0, float y0 = 0, float scale = 1)
{
	C2DVector v_temp, r_temp;
	v_temp.x = cos(theta);
	v_temp.y = sin(theta);

	RGB color;
	Find_Color(color);

	for (int i = 0; i < chain_length; i++)
	{
		C2DVector s_i = v_temp*(((1-chain_length)/2.0 + i));

		r_temp = r + s_i;

		Draw_Circle(r_temp, v_temp*(-tail), radius*scale, thickness, color, mode);
	}
}

void VisualChain::Draw_Magnified(C2DVector r0, float d0, C2DVector r1, float d1)
{
	float scale = d1 / d0;
	C2DVector v_temp, p;
	v_temp.x = cos(theta);
	v_temp.y = sin(theta);

	RGB color;
	Find_Color(color);

	for (int i = 0; i < chain_length; i++)
	{
		C2DVector s_i = v_temp*(((1-chain_length)/2.0 + i));
		C2DVector p;
		p = r - r0 + s_i;

		if (fabs(p.x) < (d0-radius) && fabs(p.y) < (d0-radius))
		{
			p *= scale;
			p += r1;
			const int thickness_factor = 3; // 2000

			Draw_Circle(p, v*(-tail), radius*scale, scale*thickness / thickness_factor, color, mode);
		}
	}
}

#endif
