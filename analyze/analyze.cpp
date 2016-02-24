#include<iostream>
#include<cstdlib>
#include<vector>

#include"analyze.h"

using namespace std;

int main(int argc, char** argv)
{
	C2DVector::Init_Rand(321);
	srand(time(NULL));

	for (int i = 1; i < argc; i++)
	{
		string name = argv[i];
		SceneSet* sceneset = new SceneSet(name);
		bool read_state = sceneset->Read();
		if (read_state)
		{
			sceneset->L -= 0.5-0.1;
			Real box_dim = sceneset->L;

			boost::replace_all(name, "-r-v.bin", "");

//			sceneset->Save_Theta_Deviation(120, 0, sceneset->scene.size(), "theta-stat.dat");
//			sceneset->Plot_Fields(21, 40, name);
//			sceneset->Plot_Averaged_Fields(128, name);
//			sceneset->Plot_Averaged_Fields(64, name);
//			sceneset->Plot_Averaged_Fields(32, name);
//			sceneset->Plot_Averaged_Fields(25, name);
//			sceneset->Plot_Averaged_Fields_Section(41, 40, name);
//			sceneset->Plot_Averaged_Fields_Section(41, 38, name);
//			sceneset->Plot_Averaged_Fields_Section(41, 20, name);
//			sceneset->Plot_Density_Contour(61, 0.1, name);


			size_t pos1,pos2;
			pos1 = name.find("-Lx");
			if (pos1 > 0 && pos1 < name.length())
				name.erase(name.begin() + pos1,name.end());
			pos1 = name.find("-2Lx");
			if (pos1 > 0 && pos1 < name.length())
				name.erase(name.begin() + pos1,name.end());
			boost::replace_all(name, "rho=", "");
			boost::replace_all(name, "-noise=", "\t");
			boost::replace_all(name, "-cooling", "\t");
			boost::replace_all(name, "-g=", "\t");
			boost::replace_all(name, "-v=", "\t");
			boost::replace_all(name, "-alpha=", "\t");

			double p,dp,sigma2,G;
			Compute_Order_Parameters(sceneset, p,dp, sigma2, G);
			cout << name << "\t" << p << "\t" << dp << "\t" << sigma2 << "\t" << G << endl;

//			Stat<double>	angular_momentum_data;
//			Compute_Angular_Momentum(sceneset, &angular_momentum_data);
//			cout << name << "\t" << sceneset->L << "\t" << (angular_momentum_data.mean) << "\t" << angular_momentum_data.error << endl;
//			angular_momentum_data.Reset();
//			cout << name << "\t" << Local_Cohesion(sceneset, 10) << endl;

//			cout << "# " << name << endl;
//			Time_AutoCorrelation(sceneset, 10);
//			Spatial_AutoCorrelation(sceneset, 50, 5);

//			int r = rand() % Scene::number_of_particles;
//			Trajectory(sceneset,r);
//			Angle_Time(sceneset, 12);
//			for (int j = 0; j < Scene::number_of_particles; j++)
//			Angular_Velocity_Time(sceneset, j);

//			Compute_Fluctuation(sceneset);

//			Radial_Density(sceneset, 200);

			// The variables: Mean_Squared_Distance_Growth(SceneSet* s, int frames, int number_of_points, int number_of_pair_sets, Real r_cut)
//			Mean_Squared_Distance_Growth(sceneset, 200, 200, 40, 0.01); // Mean_Squared_Distance_Growth(SceneSet* s, int frames, int number_of_points, int number_of_pair_sets, Real r_cut)
//			Mean_Squared_Displacement_Growth(sceneset, sceneset->scene.size(), 400);// void Mean_Squared_Displacement_Growth(SceneSet* s, int frames, int number_of_points)
//			Lyapunov_Exponent(sceneset, 900, 200, 40, 0.1,0.2);

//			Pair_Distribution(sceneset, 6,400);

			delete sceneset;
		}
		else
			cout << "Was not able to open file: " << name << endl;
	}
		

	return 0;
}

