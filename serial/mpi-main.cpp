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

#include "mpi.h"

inline void timing_information(clock_t start_time, int i_step, int total_step)
{
		clock_t current_time = clock();
		int lapsed_time = (current_time - start_time) / (CLOCKS_PER_SEC);
		int remaining_time = (lapsed_time*(total_step - i_step)) / (i_step + 1);
		cout << "\r" << round(100.0*i_step / total_step) << "% lapsed time: " << lapsed_time << " s		remaining time: " << remaining_time << " s" << flush;
}


inline Real equilibrium(Box* box, int equilibrium_step, int saving_period, ofstream& out_file)
{
	clock_t start_time, end_time;
	start_time = clock();

	cout << "equilibrium:" << endl;

	for (int i = 0; i < equilibrium_step; i+=cell_update_period)
	{
		box->Multi_Step(cell_update_period);
//		timing_information(start_time,i,equilibrium_step);
	}
	cout << "Finished" << endl;

	end_time = clock();

	Real t = (Real) (end_time - start_time) / CLOCKS_PER_SEC;
	return(t);
}


inline Real data_gathering(Box* box, int total_step, int saving_period, ofstream& out_file)
{
	clock_t start_time, end_time;
	start_time = clock();

	cout << "gathering data:" << endl;
	int saving_time = 0;
	for (int i = 0; i < total_step; i+=cell_update_period)
	{
		box->Multi_Step(cell_update_period);
//		timing_information(start_time,i,total_step);

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

	end_time = clock();

	Real t = (Real) (end_time - start_time) / CLOCKS_PER_SEC;
	return(t);
}

bool Chek_Seeds(long int seed, int thisnode, int totalnodes)
{
	long int s[totalnodes];
	bool b;
	int int_b;
	MPI_Status status;
	if (thisnode == 0)
	{
		s[0] = seed;
		for (int i = 1; i < totalnodes; i++)
			MPI_Recv(&s[i],1,MPI_LONG_INT,i,1,MPI_COMM_WORLD, &status);
	}
	else
		MPI_Send(&seed,1,MPI_LONG_INT,0,1,MPI_COMM_WORLD);

	MPI_Barrier(MPI_COMM_WORLD);

	if (thisnode == 0)
	{
		b = true;
		for (int i = 0; i < totalnodes; i++)
			for (int j = i+1; j < totalnodes; j++)
				b = b && (s[i] != s[j]);
		int_b = b;
		for (int i = 1; i < totalnodes; i++)
			MPI_Send(&int_b,1,MPI_INT,i,1,MPI_COMM_WORLD);
	}
	else
	{
		MPI_Recv(&int_b,1,MPI_INT,0,1,MPI_COMM_WORLD, &status);
		if (int_b == 1)
			b = true;
		else
			b = false;
	}

	MPI_Barrier(MPI_COMM_WORLD);

	return (b);
}

int main(int argc, char *argv[])
{
	int thisnode, totalnodes;
	MPI_Status status;
	MPI_Init(&argc, &argv);

	int number_of_realizations = argc - 4;
	Real input_rho = atof(argv[1]);
	Real input_g = atof(argv[2]);
	Real input_alpha = atof(argv[3]);

	Real noise_list[number_of_realizations];
	for (int i = 0; i < number_of_realizations; i++)
		noise_list[i] = atof(argv[i + 4]);

	MPI_Comm_size(MPI_COMM_WORLD, &totalnodes);
	MPI_Comm_rank(MPI_COMM_WORLD, &thisnode);

	long int seed = 0 + thisnode*112488;

	while (!Chek_Seeds(seed, thisnode, totalnodes))
	{
		seed = time(NULL) + thisnode*112488;
		MPI_Barrier(MPI_COMM_WORLD);		
	}


	C2DVector::Init_Rand(seed);

	Real t_eq,t_sim;

	Box box;

	for (int i = 0; i < number_of_realizations; i++)
	{
		if ((i % totalnodes) == thisnode)
		{
			box.Init(input_rho, input_g, input_alpha, noise_list[i]);
			cout << "From processor " << thisnode << " Box information is: " << box.info.str() << endl;
			stringstream address;
			address.str("");
			address << box.info.str() << "-r-v.bin";
			ofstream out_file(address.str().c_str());

			t_eq = equilibrium(&box, equilibrium_step, saving_period, out_file);
			cout << i << "	From node: " << thisnode << "\t" << "elapsed time of equilibrium is \t" << t_eq << endl;

			t_sim = data_gathering(&box, total_step, saving_period, out_file);
			cout << i << "	From node: " << thisnode << "\t" << "elapsed time of simulation is \t" << t_sim << endl;

			out_file.close();
		}

	}

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
}

