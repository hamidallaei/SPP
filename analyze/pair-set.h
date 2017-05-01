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
	Real r_cut, r2_cut, r_min, r2_min, r_max, r2_max;
	static SceneSet* sceneset;

	Pair_Set();
	~Pair_Set();

// Finding particles closer than r_cut
	bool Find_Close_Particle(Real r_cut, int t);
// Finding particles between r_min and r_max.
	bool Find_Particle(Real r_min, Real r_max, int t);
// finding mean square distance of pairs t frames in advance
	Real Find_Mean_Square_Distance(int tau);
// finding lyapunov exponent in short time interval
	Real Find_Short_Lyapunov_Exponent(int tau);
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
	return (Find_Particle(0, r_cut, t));
}

// Find pairs of particles that their distance d is r_min < d < r_max
bool Pair_Set::Find_Particle(Real input_r_min, Real input_r_max, int t)
{
	if (sceneset == NULL)
	{
		cout << "Error: sceneset is not set to point any varialbe (NULL)" << endl;
		return (false);
	}
	if (t > sceneset->Nf)
	{
		cout << "Error: The specifide time for finding pairs " << t << " is larger than the whole number of snapshots" << endl;
		return (false);
	}
	step = t;

	r_min = input_r_min;
	r2_min = r_min*r_min;
	r_max = input_r_max;
	r2_max = r_min*r_max;

	SavingVector L(sceneset->L_min);
	int grid_dim_x, grid_dim_y;
	if (L.x > L.y)
	{
		grid_dim_x = (int) (2*sceneset->L_min.x / r_max);
		grid_dim_x--;

		grid_dim_x = min(1000,grid_dim_x);
		grid_dim_y = (int) round(L.y * grid_dim_x / L.x);
	}
	else
	{
		grid_dim_y = (int) (2*sceneset->L_min.y / r_max);
		grid_dim_y--;

		grid_dim_y = min(1000,grid_dim_y);
		grid_dim_x = (int) round(L.x * grid_dim_y / L.y);
	}

	Cell** c = new Cell*[grid_dim_x];
	for (int i = 0; i < grid_dim_x; i++)
		c[i] = new Cell[grid_dim_y];
	
	for (int i = 0; i < sceneset->scene[step].Ns; i++)
	{
		int x = (int) floor((sceneset->scene[step].sparticle[i].r.x + sceneset->L_min.x)*grid_dim_x / (2*sceneset->L_min.x));
		int y = (int) floor((sceneset->scene[step].sparticle[i].r.y + sceneset->L_min.y)*grid_dim_y / (2*sceneset->L_min.y));
		c[x][y].Add(i);
	}

	for (int x = 0; x < grid_dim_x; x++)
		for (int y = 0; y < grid_dim_y; y++)
		{
			c[x][y].Add_Pairs_Self(this);
			c[x][y].Add_Pairs(&c[(x+1)%grid_dim_x][y], this);
			c[x][y].Add_Pairs(&c[(x+1)%grid_dim_x][(y+1)%grid_dim_y], this);
			c[x][y].Add_Pairs(&c[(x+1)%grid_dim_x][(y-1+grid_dim_y)%grid_dim_y], this);
			c[x][y].Add_Pairs(&c[x][(y+1)%grid_dim_y], this);
		}

	for (int i = 0; i < grid_dim_x; i++)
		delete [] c[i];

	delete [] c;

	return (true);
}

Real Pair_Set::Find_Mean_Square_Distance(int tau)
{
	if ((tau+step) > sceneset->Nf)
	{
		cout << "Error: The specifide time for finding distances " << tau+step << " is longer than the whole number of snapshots" << endl;
		return (-1);
	}
	Real sum_d2 = 0;
	for (int n = 0; n < pid1.size(); n++)
	{
		SavingVector dr = sceneset->scene[tau+step].sparticle[pid1[n]].r - sceneset->scene[tau+step].sparticle[pid2[n]].r;
		dr.x -= 2*sceneset->L.x*((int) (dr.x / sceneset->L_min.x));
		dr.y -= 2*sceneset->L.y*((int) (dr.y / sceneset->L_min.y));
		Real d2 = dr.Square();
		sum_d2 += d2;
	}
	sum_d2 /= pid1.size();
	return (sum_d2);
}

Real Pair_Set::Find_Short_Lyapunov_Exponent(int tau)
{
	if ((tau+step) > sceneset->Nf)
	{
		cout << "Error: The specifide time for finding distances " << tau+step << " is longer than the whole number of snapshots" << endl;
		return (-1);
	}
	Real sum_lambda = 0;
	for (int n = 0; n < pid1.size(); n++)
	{
		SavingVector dr_tau = sceneset->scene[tau+step].sparticle[pid1[n]].r - sceneset->scene[tau+step].sparticle[pid2[n]].r;
		dr_tau.x -= 2*sceneset->L.x*((int) (dr_tau.x / sceneset->L_min.y));
		dr_tau.y -= 2*sceneset->L.y*((int) (dr_tau.y / sceneset->L_min.y));
		SavingVector dr0 = sceneset->scene[step].sparticle[pid1[n]].r - sceneset->scene[step].sparticle[pid2[n]].r;
		dr0.x -= 2*sceneset->L.x*((int) (dr0.x / sceneset->L_min.x));
		dr0.y -= 2*sceneset->L.y*((int) (dr0.y / sceneset->L_min.y));
		Real lambda = dr_tau.Square() / dr0.Square();
		lambda = 0.5 * (log(lambda) / tau);
		sum_lambda += lambda;
	}
	sum_lambda /= pid1.size();
	return (sum_lambda);
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
			SavingVector dr = ps->sceneset->scene[ps->step].sparticle[pid[i]].r - ps->sceneset->scene[ps->step].sparticle[pid[j]].r;
			dr.x -= 2*ps->sceneset->L.x*((int) (dr.x / ps->sceneset->L_min.y));
			dr.y -= 2*ps->sceneset->L.y*((int) (dr.y / ps->sceneset->L_min.y));
			Real d2 = dr.Square();
			if (d2 < ps->r2_max && d2 > ps->r2_min)
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
			SavingVector dr = ps->sceneset->scene[ps->step].sparticle[pid[i]].r - ps->sceneset->scene[ps->step].sparticle[c->pid[j]].r;
			dr.x -= 2*ps->sceneset->L_min.x*((int) (dr.x / ps->sceneset->L_min.x));
			dr.y -= 2*ps->sceneset->L_min.y*((int) (dr.y / ps->sceneset->L_min.y));
			Real d2 = dr.Square();
			if (d2 < ps->r2_cut)
			{
				ps->pid1.push_back(pid[i]);
				ps->pid2.push_back(c->pid[j]);
			}
		}
}

#endif
