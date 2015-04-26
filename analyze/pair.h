#ifndef _PAIRSET_
#define _PAIRSET_

#include "../shared/c2dvector.h"
#include <boost/tuple/tuple.hpp>
#include <vector>

class Pair_Set{
public:
	vector<int> pid1;
	vector<int> pid2;
	int step;
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

bool Pair_Set::Find_Close_Particle(Real r_cut, int t)
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
	Real r2_cut = r_cut*r_cut;
	for (int i = 0; i < Scene::number_of_particles; i++)
		for (int j = i+1; j < Scene::number_of_particles; j++)
		{
			C2DVector dr = sceneset->scene[step].particle[i].r - sceneset->scene[step].particle[j].r;
			dr.x -= 2*sceneset::L*((int) (dr.x / sceneset::L));
			dr.y -= 2*sceneset::L*((int) (dr.y / sceneset::L));
			Real d2 = dr.Square();
			if (d2 < r2_cut)
			{
				pid1.push_back(i);
				pid2.push_back(j);
			}
		}
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
		dr.x -= 2*sceneset::L*((int) (dr.x / sceneset::L));
		dr.y -= 2*sceneset::L*((int) (dr.y / sceneset::L));
		Real d2 = dr.Square();
		sum_d2 += d2;
	}
	sum_d2 /= pid1.size();
	return (sum_d2);
}

#endif
