#ifndef _GEOMETRY_
#define _GEOMETRY_

#include "c2dvector.h"

class Geometry{
public:
	int wall_num;
	Wall wall[max_wall_num];
	C2DVector center;
	Real total_length;

	Geometry();
	void Reset();
	void Add_Wall(C2DVector point_1, C2DVector point_2);
	void Find_Center();
	void Translate(C2DVector delta);
	void Rotate(Real phi);
};

Geometry::Geometry()
{
}

void Geometry::Reset()
{
	wall_num = 0;
}

void Add_Wall(C2DVector point_1, C2DVector point_2)
{
	wall[wall_num].Init(point_1, point_2);
	total_length += wall[wall_num].length;
	wall_num++;
}

void Find_Center()
{
	center.Null();
	for (int i = 0; i < wall_num; i++)
		center += (wall_num[i].point_1 + wall_num[i].point_2)*(wall_num[i].length/2);
	center = center / total_length;
}

void Translate(C2DVector delta)
{
	for (int i = 0; i < wall_num; i++)
		wall[i].Translate(delta);
}

void Rotate(Real phi)
{
	for (int i = 0; i < wall_num; i++)
		wall[i].Rotate(phi,center);
}

#endif
