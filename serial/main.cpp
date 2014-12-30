#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <ctime>
#include <algorithm>
#include "c2dvector.h"
#include "parameters.h"
#include "particle.h"
#include "cell.h"
#include "box.h"

inline void timing_information(clock_t start_time, int i_step, int total_step)
{
		clock_t current_time = clock();
		int lapsed_time = (current_time - start_time) / (CLOCKS_PER_SEC);
		int remaining_time = (lapsed_time*(total_step - i_step)) / (i_step + 1);
		cout << "\r" << round(100.0*i_step / total_step) << "% lapsed time: " << lapsed_time << " s		remaining time: " << remaining_time << " s" << flush;
}


inline void equilibrium(Box* box, int equilibrium_step, int saving_period, ofstream& out_file)
{
	clock_t start_time = clock();
	cout << "equilibrium:" << endl;
//	box->Init();
	for (int i = 0; i < equilibrium_step; i+=cell_update_period)
	{
		box->Multi_Step(cell_update_period);
		timing_information(start_time,i,equilibrium_step);
	}
	cout << "Finished" << endl;
}


inline void data_gathering(Box* box, int total_step, int saving_period, ofstream& out_file)
{
	clock_t start_time = clock();
	cout << "gathering data:" << endl;
	int saving_time = 0;
	for (int i = 0; i < total_step; i+=cell_update_period)
	{
		box->Multi_Step(cell_update_period);
		timing_information(start_time,i,total_step);

		if ((i / cell_update_period) % saving_period == 0)
		{
			out_file << box;

//			if ((i / saving_period) % 100 == 0)
//			{
//				stringstream address;
//				address.str("");
//				address << "data-" << saving_time << ".dat";
//				ofstream velocity_field_file(address.str().c_str());
//				box->Save_Fields(velocity_field_file);
//				velocity_field_file.close();
//			}

//			saving_time++;
		}
	}

//	stringstream address;
//	address.str("");
//	address << "data-" << box->info.str() << ".dat";
//	ofstream velocity_field_file(address.str().c_str());
//	velocity_field_file.close();

	cout << "Finished" << endl;
}


int main(int argc, char *argv[])
{
	C2DVector::Init_Rand(time(NULL));
//	C2DVector::Init_Rand(321);

	Box box;

	box.Init(atof(argv[1]), atof(argv[2]), atof(argv[3]), atof(argv[4]));

	cout << box.info.str() << endl;

	stringstream address;
	address.str("");
	address << box.info.str() << "-r-v.bin";
	ofstream out_file(address.str().c_str());

	equilibrium(&box, equilibrium_step, saving_period, out_file);
	data_gathering(&box, total_step, saving_period, out_file);
	out_file.close();
}

