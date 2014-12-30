#ifndef _CELL_
#define _CELL_

#include "c2dvector.h"
#include "parameters.h"
#include <vector>

class Cell{
public:
	vector<int> pid; // particle_id
	C2DVector r; // Center position of the cell in the box
	static C2DVector dim; // Dimension of the cell width and height
	static Particle* particle; // This is a pointer to the original particle array pointer of the box. We need this pointer in some subroutins

	Cell();

	void Init();
	void Init(Real x, Real y);
	void Delete();
	void Add(int p); // Add a particle id to the list of pid of this cell.
	void Interact(Cell* c); // Interact all particles wihtin this cell with the cell c
	void Self_Interact(); // Interact all particles within this cell with themselve
	void Move();
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
	pid.clear();
}

void Cell::Add(int p)
{
	pid.push_back(p);
}

void Cell::Interact(Cell* c)
{
	for (int i = 0; i < pid.size(); i++)
	{
		for (int j = 0; j < c->pid.size(); j++)
			particle[pid[i]].Interact(particle[c->pid[j]]);
	}
}

void Cell::Self_Interact()
{
	for (int i = 0; i < pid.size(); i++)
	{
		for (int j = i+1; j < pid.size(); j++)
			particle[pid[i]].Interact(particle[pid[j]]);
	}
}

void Cell::Move()
{
	for (int i = 0; i < pid.size(); i++)
		particle[pid[i]].Move();
}

C2DVector Cell::dim;
Particle* Cell::particle = NULL; // Be carefull that this pointer be initiated in future

#endif

