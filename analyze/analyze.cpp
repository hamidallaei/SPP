#include<iostream>
#include<cstdlib>
#include<vector>

#include"analyze.h"

using namespace std;

int main(int argc, char** argv)
{
	SavingVector::Init_Rand(321);
	srand(time(NULL));
	for (int i = 1; i < argc; i++)
	{
		string name = argv[i];
		SceneSet* sceneset = new SceneSet(name);
		bool read_state = sceneset->Read();
		if (read_state)
		{
//			sceneset->L -= 0.5-0.1;
			SavingVector box_dim(sceneset->L);

			boost::replace_all(name, "-r-v.bin", "");
			stringstream ss("");
			ss << "quantities-" << name << ".dat";
			ofstream out_file;
			out_file.open(ss.str().c_str());

//			sceneset->Save_Theta_Deviation(120, 0, sceneset->Nf, "theta-stat.dat");
//			sceneset->Plot_Fields(21, 40, name);
//			sceneset->Plot_Averaged_Fields(128, name);
//			sceneset->Plot_Averaged_Fields(64, name);
//			sceneset->Plot_Averaged_Fields(32, name);
//			sceneset->Plot_Averaged_Fields(25, name);
//			sceneset->Plot_Averaged_Fields_Section(41, 40, name);
//			sceneset->Plot_Averaged_Fields_Section(41, 38, name);
//			sceneset->Plot_Averaged_Fields_Section(41, 20, name);
//			sceneset->Plot_Density_Contour(61, 0.1, name);

//			std::size_t pos = name.find("Dr=");
//			std::string str = name.substr (pos);
//			pos = str.find("-");
//			str = str.substr(3,pos-3);
//			stringstream ss;
//			ss.str("");
//			ss << "theta_stat-Dr-" << str;

//			double p_c = atof(argv[2]);
//			double dp = atof(argv[3]);
//			sceneset->Accumulate_Theta(16, 30, p_c, dp, ss.str().c_str());
//			sceneset->Accumulate_Theta(16, 30, 0.01, dp, ss.str().c_str());
//			sceneset->Accumulate_Theta(16, 30, 0.05, dp, ss.str().c_str());
//			sceneset->Accumulate_Theta(16, 30, 0.1, dp, ss.str().c_str());
//			sceneset->Accumulate_Theta(16, 30, 0.2, dp, ss.str().c_str());
//			sceneset->Accumulate_Theta(16, 30, 0.3, dp, ss.str().c_str());
//			sceneset->Accumulate_Theta(16, 30, 0.4, dp, ss.str().c_str());
//			sceneset->Accumulate_Theta(16, 30, 0.5, dp, ss.str().c_str());
//			sceneset->Accumulate_Theta(16, 30, 0.6, dp, ss.str().c_str());
//			sceneset->Accumulate_Theta(16, 30, 0.7, dp, ss.str().c_str());
//			sceneset->Accumulate_Theta(16, 30, 0.8, dp, ss.str().c_str());
//			sceneset->Accumulate_Theta(16, 30, 0.9, dp, ss.str().c_str());
//			sceneset->Accumulate_Theta(16, 30, 0.95, dp, ss.str().c_str());
//			sceneset->Accumulate_Theta(16, 30, 0.99, dp, ss.str().c_str());


//			size_t pos1,pos2;
//			pos1 = name.find("-Lx");
//			if (pos1 > 0 && pos1 < name.length())
//				name.erase(name.begin() + pos1,name.end());
//			pos1 = name.find("-2Lx");
//			if (pos1 > 0 && pos1 < name.length())
//				name.erase(name.begin() + pos1,name.end());\
//			pos1 = name.find("-2L");
//			if (pos1 > 0 && pos1 < name.length())
//				name.erase(name.begin() + pos1,name.end());
//			pos1 = name.find("-v=");
//			pos2 = name.find("-noise=");
//			if (pos1 > 0 && pos1 < name.length())
//				name.erase(name.begin() + pos1 + 3,name.begin() + pos2);

//			boost::replace_all(name, "rho=", "");
//			boost::replace_all(name, "-noise=", "\t");
//			boost::replace_all(name, "-cooling", "\t");
//			boost::replace_all(name, "-g=", "\t");
//			boost::replace_all(name, "-v=", "\t");
//			boost::replace_all(name, "-alpha=", "\t");

//			Polarization_AutoCorr(sceneset);
//			Polarization_Time(sceneset);

			cout << name << endl;
			Quantities_Time(sceneset, 0*sceneset->Nf/2, out_file);

//			Stat<double> p;
//			Compute_Polarization(sceneset,&p);

//			double p,dp,sigma2,G;
//			Compute_Order_Parameters(sceneset, p,dp, sigma2, G);
//			cout << name << "\t" << p << "\t" << dp << "\t" << sigma2 << "\t" << G << endl;

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
//			Mean_Squared_Displacement_Growth(sceneset, sceneset->Nf, 400);// void Mean_Squared_Displacement_Growth(SceneSet* s, int frames, int number_of_points)
//			Lyapunov_Exponent(sceneset, 900, 200, 40, 0.1,0.2);

//			Pair_Distribution(sceneset, 6,400);



			delete sceneset;
		}
		else
			cout << "Was not able to open file: " << name << endl;
	}
		

	return 0;
}

