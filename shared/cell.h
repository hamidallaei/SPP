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
	void Clear_Neighbor_List(); // This will clean the neighbor list of particles inside this cell
	void Neighbor_List(); // Adding neighboring particles to their list in a cell but each pair of close particles are presented only one time as a member of neighbor list of one of the pair particles.
	void Neighbor_List(Cell* c); // Adding neighboring particles of different cells to the neighbor list of particles.
	void Interact(); // Interacting using nieghbor list. Particles outside of this cell are also considered.
	void Interact(Cell* c); // Interact all particles wihtin this cell with the cell c
	void Self_Interact(); // Interact all particles within this cell with themselve
	void Move();

	#ifdef RUNGE_KUTTA
	void Move_Runge_Kutta_1();
	void Move_Runge_Kutta_2();
	#endif

	C2DVector Compute_Polarization_Sum(); // Compute polarization of the particles inside the cell
};

Cell::Cell()
{
	if (((Lx2 / divisor_x) < rv) || ((Ly2 / divisor_y) < rv))
	{
		cout << "Error, cells are too small that befor a cell update the run particles neighbors change. You have to decrease number of cells (divisor_x or divisor_y)" << endl;
		exit(0);
	}
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

void Cell::Clear_Neighbor_List()
{
	for (int i = 0; i < pid.size(); i++)
		particle[pid[i]].neighbor_id.clear();
}

void Cell::Neighbor_List(Cell* c)
{
	for (int i = 0; i < pid.size(); i++)
	{
		for (int j = 0; j < c->pid.size(); j++)
		{
			C2DVector dr = particle[pid[i]].r - particle[c->pid[j]].r;
			#ifdef PERIODIC_BOUNDARY_CONDITION
				dr.Periodic_Transform();
			#endif
			Real d = sqrt(dr.Square());
			if (d < rv)
				particle[pid[i]].neighbor_id.push_back(c->pid[j]);
		}
	}
}

void Cell::Neighbor_List()
{
	for (int i = 0; i < pid.size(); i++)
	{
		for (int j = i+1; j < pid.size(); j++)
		{
			C2DVector dr = particle[pid[i]].r - particle[pid[j]].r;
			Real d = sqrt(dr.Square());
			if (d < rv)
				particle[pid[i]].neighbor_id.push_back(pid[j]);
		}
	}
}

void Cell::Interact()
{
	for (int i = 0; i < pid.size(); i++)
		for (int j = 0; j < particle[pid[i]].neighbor_id.size(); j++)
			particle[pid[i]].Interact(particle[particle[pid[i]].neighbor_id[j]]);
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

#ifdef RUNGE_KUTTA
void Cell::Move_Runge_Kutta_1()
{
	for (int i = 0; i < pid.size(); i++)
		particle[pid[i]].Move_Runge_Kutta_1();
}

void Cell::Move_Runge_Kutta_2()
{
	for (int i = 0; i < pid.size(); i++)
		particle[pid[i]].Move_Runge_Kutta_2();
}
#endif

// Compute polarization of the particles inside the cell
C2DVector Cell::Compute_Polarization_Sum()
{
	C2DVector vp_sum;
	for (int i = 0; i < pid.size(); i++)
		vp_sum += particle[pid[i]].v;
	return vp_sum;
}

C2DVector Cell::dim;
Particle* Cell::particle = NULL; // Be carefull that this pointer be initiated in future

#endif

