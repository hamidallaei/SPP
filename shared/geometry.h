#ifndef _GEOMETRY_
#define _GEOMETRY_

#include "c2dvector.h"
#include "wall.h"

class Geometry{
public:
	int wall_num;
	Wall wall[max_wall_num];
	C2DVector center;
	Real total_length;

	Geometry();
	void Reset();
	void Add_Wall(C2DVector point_1, C2DVector point_2);
	void Add_Wall(Real point_1_x, Real point_1_y, Real point_2_x, Real point_2_y);
	void Find_Center();
	void Translate(C2DVector delta);
	void Rotate(Real phi);

	void Interact(Particle* p);
};

Geometry::Geometry()
{
}

void Geometry::Reset()
{
	wall_num = 0;
}

void Geometry::Add_Wall(C2DVector point_1, C2DVector point_2)
{
	wall[wall_num].Init(point_1, point_2);
	total_length += wall[wall_num].length;
	wall_num++;
}


void Geometry::Add_Wall(Real point_1_x, Real point_1_y, Real point_2_x, Real point_2_y)
{
	wall[wall_num].Init(point_1_x, point_1_y, point_2_x, point_2_y);
	total_length += wall[wall_num].length;
	wall_num++;
}


void Geometry::Find_Center()
{
	center.Null();
	for (int i = 0; i < wall_num; i++)
		center += (wall[i].point_1 + wall[i].point_2)*(wall[i].length/2);
	center = center / total_length;
}

void Geometry::Translate(C2DVector delta)
{
	for (int i = 0; i < wall_num; i++)
		wall[i].Translate(delta);
}

void Geometry::Rotate(Real phi)
{
	for (int i = 0; i < wall_num; i++)
		wall[i].Rotate(phi, center);
}

void Geometry::Interact(Particle* p)
{
	for (int i = 0; i < wall_num; i++)
		wall[i].Interact(p);
}


#endif
