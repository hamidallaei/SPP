#include "c2dvector.h"
#include "parameters.h"
#include "particle.h"
#include "cell.h"
#include "box.h"

inline void timing_information(Node* node, clock_t start_time, int i_step, int total_step)
{
	if (node->node_id == 0)
	{
		clock_t current_time = clock();
		int lapsed_time = (current_time - start_time) / (CLOCKS_PER_SEC);
		int remaining_time = (lapsed_time*(total_step - i_step)) / (i_step + 1);
		cout << "\r" << round(100.0*i_step / total_step) << "% lapsed time: " << lapsed_time << " s		remaining time: " << remaining_time << " s" << flush;
	}
	MPI_Barrier(MPI_COMM_WORLD);
}


inline Real equilibrium(Box* box, int equilibrium_step, int saving_period, ofstream& out_file)
{
	clock_t start_time, end_time;
	start_time = clock();

	if (box->thisnode->node_id == 0)
		cout << "equilibrium:" << endl;

	for (int i = 0; i < equilibrium_step; i+=cell_update_period)
	{
		box->Multi_Step(cell_update_period);
		timing_information(box->thisnode,start_time,i,equilibrium_step);
	}

	if (box->thisnode->node_id == 0)
		cout << "Finished" << endl;

	end_time = clock();

	Real t = (Real) (end_time - start_time) / CLOCKS_PER_SEC;
	return(t);
}


inline Real data_gathering(Box* box, int total_step, int saving_period, ofstream& out_file)
{
	clock_t start_time, end_time;
	start_time = clock();

	#ifdef TRACK_PARTICLE
	if (!flag)
		flag = true;
	#endif

	if (box->thisnode->node_id == 0)
		cout << "gathering data:" << endl;
	int saving_time = 0;

	for (int i = 0; i < total_step; i+=cell_update_period)
	{
		box->Multi_Step(cell_update_period);
		timing_information(box->thisnode,start_time,i,total_step);
		if ((i / cell_update_period) % saving_period == 0)
			out_file << box;
	}

	if (box->thisnode->node_id == 0)
		cout << "Finished" << endl;

	end_time = clock();

	Real t = (Real) (end_time - start_time) / CLOCKS_PER_SEC;
	return(t);
}



int main(int argc, char *argv[])
{
	int this_node_id, total_nodes;
	MPI_Status status;
	MPI_Init(&argc, &argv);

	Real input_rho = atof(argv[1]);
	Real input_g = atof(argv[2]);
	Real input_alpha = atof(argv[3]);
	Real noise_amplitude = atof(argv[4]);

	Node thisnode;

//	thisnode.seed = 0 + thisnode.node_id*112488;
	thisnode.seed = seed;

//	while (!thisnode.Chek_Seeds())
//	{
//		thisnode.seed = time(NULL) + thisnode.node_id*112488;
//		MPI_Barrier(MPI_COMM_WORLD);
//	}

	C2DVector::Init_Rand(seed);

	Real t_eq,t_sim;

	Box box;
	box.Init(&thisnode, input_rho, input_g, input_alpha, noise_amplitude);

	cout << "From processor " << thisnode.node_id << " Box information is: " << box.info.str() << endl;

	ofstream out_file;
	if (thisnode.node_id == 0)
	{
			stringstream address;
			address.str("");
			address << box.info.str() << "-r-v.bin";
			out_file.open(address.str().c_str());
	}
	MPI_Barrier(MPI_COMM_WORLD);

	t_eq = equilibrium(&box, equilibrium_step, saving_period, out_file);
	MPI_Barrier(MPI_COMM_WORLD);
	t_sim = data_gathering(&box, total_step, saving_period, out_file);
	MPI_Barrier(MPI_COMM_WORLD);

//	if (thisnode.node_id == 0)
//		for (int i = 0; i < box.N; i++)
//			cout << i << "\t" << box.particle[i].r << endl;

	if (thisnode.node_id == 0)
		out_file.close();

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
}

