#ifndef _PAIR_
#define _PAIR_

#include "../shared/c2dvector.h"
#include <boost/tuple/tuple.hpp>
#include <vector>

class Pair_Set;

class Cell{
public:
	vector<int> pid; // particle_id

	Cell();
	~Cell();

	void Add(int p); // Add a particle id to the list of pid of this cell.
	void Add_Pairs_Self(Pair_Set* ps);
	void Add_Pairs(Cell* c, Pair_Set* ps);
};

class Pair_Set{
public:
	vector<int> pid1;
	vector<int> pid2;
	int step;
	Real r_cut, r2_cut;
	static SceneSet* sceneset;

	Pair_Set();
	~Pair_Set();

	bool Find_Close_Particle(Real r_cut, int t);
// finding mean square distance of pairs t frames in advance
	Real Find_Mean_Square_Distance(int tau);
};

SceneSet* Pair_Set::sceneset = NULL;

Pair_Set::Pair_Set()
{
}

Pair_Set::~Pair_Set()
{
	pid1.clear();
	pid2.clear();
}

bool Pair_Set::Find_Close_Particle(Real input_r_cut, int t)
{
	if (sceneset == NULL)
	{
		cout << "Error: sceneset is not set to point any varialbe (NULL)" << endl;
		return (false);
	}
	if (t > sceneset->scene.size())
	{
		cout << "Error: The specifide time for finding pairs " << t << " is larger than the whole number of snapshots" << endl;
		return (false);
	}
	step = t;

	int grid_dim = (int) (2*sceneset->L_min / r_cut);
	grid_dim--;

	grid_dim = min(1000,grid_dim);

	Cell** c = new Cell*[grid_dim];
	for (int i = 0; i < grid_dim; i++)
		c[i] = new Cell[grid_dim];
	
	for (int i = 0; i < Scene::number_of_particles; i++)
	{
		int x = (int) floor((sceneset->scene[step].particle[i].r.x + sceneset->L_min)*grid_dim / (2*sceneset->L_min));
		int y = (int) floor((sceneset->scene[step].particle[i].r.y + sceneset->L_min)*grid_dim / (2*sceneset->L_min));
		c[x][y].Add(i);
	}

	r_cut = input_r_cut;
	r2_cut = r_cut*r_cut;

	for (int x = 0; x < grid_dim; x++)
		for (int y = 0; y < grid_dim; y++)
		{
			c[x][y].Add_Pairs_Self(this);
			c[x][y].Add_Pairs(&c[(x+1)%grid_dim][y], this);
			c[x][y].Add_Pairs(&c[(x+1)%grid_dim][(y+1)%grid_dim], this);
			c[x][y].Add_Pairs(&c[(x+1)%grid_dim][(y-1+grid_dim)%grid_dim], this);
			c[x][y].Add_Pairs(&c[x][(y+1)%grid_dim], this);
		}

	for (int i = 0; i < grid_dim; i++)
		delete [] c[i];

	delete [] c;

	return (true);
}

Real Pair_Set::Find_Mean_Square_Distance(int tau)
{
	if ((tau+step) > sceneset->scene.size())
	{
		cout << "Error: The specifide time for finding distances " << tau+step << " is larger than the whole number of snapshots" << endl;
		return (-1);
	}
	Real sum_d2 = 0;
	for (int n = 0; n < pid1.size(); n++)
	{
		C2DVector dr = sceneset->scene[tau+step].particle[pid1[n]].r - sceneset->scene[tau+step].particle[pid2[n]].r;
		dr.x -= 2*sceneset->L*((int) (dr.x / sceneset->L_min));
		dr.y -= 2*sceneset->L*((int) (dr.y / sceneset->L_min));
		Real d2 = dr.Square();
		sum_d2 += d2;
	}
	sum_d2 /= pid1.size();
	return (sum_d2);
}

// Deffinition of cell class


Cell::Cell()
{
}

Cell::~Cell()
{
	pid.clear();
}

void Cell::Add(int p)
{
	pid.push_back(p);
}

void Cell::Add_Pairs_Self(Pair_Set* ps)
{
	for (int i = 0; i < pid.size(); i++)
		for (int j = i+1; j < pid.size(); j++)
		{
			C2DVector dr = ps->sceneset->scene[ps->step].particle[pid[i]].r - ps->sceneset->scene[ps->step].particle[pid[j]].r;
			dr.x -= 2*ps->sceneset->L*((int) (dr.x / ps->sceneset->L_min));
			dr.y -= 2*ps->sceneset->L*((int) (dr.y / ps->sceneset->L_min));
			Real d2 = dr.Square();
			if (d2 < ps->r2_cut)
			{
				ps->pid1.push_back(pid[i]);
				ps->pid2.push_back(pid[j]);
			}
		}
}

void Cell::Add_Pairs(Cell* c, Pair_Set* ps)
{
	for (int i = 0; i < pid.size(); i++)
		for (int j = 0; j < c->pid.size(); j++)
		{
			C2DVector dr = ps->sceneset->scene[ps->step].particle[pid[i]].r - ps->sceneset->scene[ps->step].particle[c->pid[j]].r;
			dr.x -= 2*ps->sceneset->L_min*((int) (dr.x / ps->sceneset->L_min));
			dr.y -= 2*ps->sceneset->L_min*((int) (dr.y / ps->sceneset->L_min));
			Real d2 = dr.Square();
			if (d2 < ps->r2_cut)
			{
				ps->pid1.push_back(pid[i]);
				ps->pid2.push_back(c->pid[j]);
			}
		}
}

#endif
