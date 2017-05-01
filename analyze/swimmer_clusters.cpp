#include<iostream>
#include<cstdlib>
#include<cmath>
#include<vector>
#include <algorithm>    // sort

#include"analyze.h"

using namespace std;

/* 
This code reads position data of the swimmers' particles from -r-v.bin files and finds the number of clusters that these particles make inside the vesicle, according to both their angular positions (wrt the vesicle's center of mass) and their velocity directions.
*/


struct Angle{
	float theta;
	float theta_position;
	int id;
	Angle();
	Angle(Angle&);
	void operator = (const Angle &angle)
	{
		theta = angle.theta;
		id = angle.id;
	}
};

Angle::Angle()
{
}

Angle::Angle(Angle& angle)
{
	id = angle.id;
	theta = angle.theta;
}

int compare (const void * a, const void * b)
{
	float diff = ((Angle*) a)->theta - ((Angle*) b)->theta;
	int sgn = (diff > 0) ? 1 : ((diff < 0) ? -1 : 0);
	return ( sgn );
}

int compare_position (const void * a, const void * b)
{
	float diff = ((Angle*) a)->theta_position - ((Angle*) b)->theta_position;
	int sgn = (diff > 0) ? 1 : ((diff < 0) ? -1 : 0);
	return ( sgn );
}


int main(int argc, char** argv)
{
	SavingVector::Init_Rand(321);
	srand(time(NULL));

	stringstream ss("");
	ss << "clusters_phi.dat";
	ofstream out_file;
	out_file.open(ss.str().c_str());

	for (int i = 1; i < argc; i++)
	{
		string name = argv[i];
		SceneSet* sceneset = new SceneSet(name);
		bool read_state = sceneset->Read();
		if (read_state)
		{
			SavingVector box_dim(sceneset->L);

			Angle angles[sceneset->scene[0].Ns];

//			cout << sceneset->scene.size() << endl;  //8193

			for (int j = 1000; j < (sceneset->scene.size()-5); j+=4) 
			/* I subtract 5 and jump every 4 steps to equalize dim of output file (curvature file) with the files extracted from quant*; e.g. omegas, speed, .... j=1000 corresponds to from_row=251 in quant* */
			{
//			int j = 35;

				//// sort particles according to their angular positions and compute number of jumps(clusters)
				SavingVector rcm;
				rcm.Null();

				//// membrane's center of mass
				for (int k = 0; k < sceneset->scene[j].Nm; k++)
					rcm += sceneset->scene[j].mparticle[k].r;
				rcm /= sceneset->scene[j].Nm;

				for (int k = 0; k < sceneset->scene[j].Ns; k++)
				{
					angles[k].id = k;
					angles[k].theta = sceneset->scene[j].sparticle[k].theta;
					SavingVector dr = (sceneset->scene[j].sparticle[k].r - rcm);
					angles[k].theta_position = atan2(dr.y,dr.x);
				}
				qsort (angles, sceneset->scene[j].Ns, sizeof(Angle), compare_position);

				int nc_position = 0;
				for (int k = 0; k < sceneset->scene[j].Ns; k++)
				{
					float dtheta = angles[(k+1) % sceneset->scene[j].Ns].theta_position - angles[k].theta_position;
					if (dtheta > M_PI)
						dtheta -= 2*M_PI;
					if (dtheta < -M_PI)
						dtheta += 2*M_PI;
					dtheta = fabs(dtheta);
//					cout << angles[k].theta_position << "\t" << dtheta << endl;
					if (dtheta > 0.4)
					{
//						cout << "dtheta > 0.4 :\t" << dtheta << endl;
						nc_position++;
//						cout << "nc_position = \t" << nc_position << endl;
					}
				}
				//// if it hasn't find a jump, set nc_position equal to 1, which refers to a one large cluster.
				if (nc_position == 0)
					nc_position =1;

				//// sort particles according to their directions and compute number of jumps(clusters)
				qsort (angles, sceneset->scene[j].Ns, sizeof(Angle), compare);

				int nc_direction = 0;
				for (int k = 0; k < sceneset->scene[j].Ns; k++)
				{
					float dtheta = angles[(k+1) % sceneset->scene[j].Ns].theta - angles[k].theta;
					if (dtheta > M_PI)
						dtheta -= 2*M_PI;
					if (dtheta < -M_PI)
						dtheta += 2*M_PI;
					dtheta = fabs(dtheta);
//					cout << angles[k].theta << "\t" << dtheta << endl;
					if (dtheta > 0.75)
					{
//						cout << "dtheta > 0.75 :\t" << dtheta << endl;
						nc_direction++;
//						cout << "nc_direction = \t" << nc_direction << endl;
					}
				}
				//// if it hasn't find a jump, set nc_direction equal to 1, which refers to a one large cluster.
				if (nc_direction == 0)
					nc_direction =1;

//				for (int k = 0; k < sceneset->scene[j].Ns; k++)
//				{
//					cout << angles[k].theta << endl;
//				}
				out_file << nc_direction << "\t" << nc_position << endl;
//				out_file << sceneset->scene[j].t << "\t" << nc_direction << "\t" << nc_position << "\t" << max(nc_direction,nc_position) << endl;
			}

			delete sceneset;
		}
		else
			cout << "Was not able to open file: " << name << endl;
	}

	return 0;
}

