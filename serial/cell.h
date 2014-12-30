#ifndef _CELL_
#define _CELL_

#include "c2dvector.h"
#include "parameters.h"
#include <vector>

class Cell{
public:
	vector<Particle*> particle;
	C2DVector r;
	static C2DVector dim;

	Cell();

	void Init();
	void Init(Real x, Real y);
	void Delete();
	void Add(Particle* p);
	void Interact(Cell* c);
	void Self_Interact();
	bool Inside(C2DVector);
};

Cell::Cell()
{
}

void Cell::Init(Real x, Real y)
{
	r.x = x;
	r.y = y;
}

void Cell::Delete()
{
	particle.clear();
}

void Cell::Add(Particle* p)
{
	particle.push_back(p);
}

//void Cell::Add(Wall w)
//{
//	bool b = Inside(w);
//	for (int i = 0; i < wall.size(); i++)
//		b = b && (wall[i] != &w);
//	if (b)
//		wall.push_back(&w);
//}

void Cell::Interact(Cell* c)
{
	for (int i = 0; i < particle.size(); i++)
	{
		for (int j = 0; j < c->particle.size(); j++)
			particle[i]->Interact(c->particle[j]);
	}
}

void Cell::Self_Interact()
{
	for (int i = 0; i < particle.size(); i++)
	{
		for (int j = i+1; j < particle.size(); j++)
			particle[i]->Interact(particle[j]);
	}
}

bool Cell::Inside(C2DVector point)
{
	C2DVector dr = point - r;
	#ifdef PERIODIC_BOUNDARY_CONDITION
		dr.Periodic_Transform();
	#endif
	bool result = ((dr.x < dim.x) && (-dim.x < dr.x)) && ((dr.y < dim.y) && (-dim.y < dr.y));
	return (result);
}


C2DVector Cell::dim;

#endif
