#include "../shared/parameters.h"
#include "../shared/c2dvector.h"
#include "../shared/particle.h"
#include "../shared/cell.h"
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


inline Real equilibrium(Box* box, long int equilibrium_step, int saving_period, ofstream& out_file)
{
	clock_t start_time, end_time;
	start_time = clock();

	if (box->thisnode->node_id == 0)
		cout << "equilibrium:" << endl;

	for (long int i = 0; i < equilibrium_step; i+=cell_update_period)
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


inline Real data_gathering(Box* box, long int total_step, int saving_period, ofstream& out_file)
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

	for (long int i = 0; i < total_step; i+=cell_update_period)
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

void Change_Noise(int argc, char *argv[], Node* thisnode)
{
	Real input_rho = atof(argv[1]);
	Real input_g = atof(argv[2]);
	Real input_alpha = atof(argv[3]);

	vector<Real> noise_list;
	for (int i = 4; i < argc; i++)
		noise_list.push_back(atof(argv[i]));

	Real t_eq,t_sim;

	Box box;
	box.Init(thisnode, input_rho, input_g, input_alpha, 0);

	ofstream out_file;

	for (int i = 0; i < noise_list.size(); i++)
	{
		Particle::noise_amplitude = noise_list[i] / sqrt(dt); // noise amplitude depends on the step (dt) because of ito calculation. If we have epsilon in our differential equation and we descritise it with time steps dt, the noise in each step that we add is epsilon times sqrt(dt) if we factorise it with a dt we have dt*(epsilon/sqrt(dt)).
		box.info.str("");
		box.info << "rho=" << box.density <<  "-g=" << Particle::g << "-alpha=" << Particle::alpha << "-noise=" << noise_list[i] << "-cooling";

		if (thisnode->node_id == 0)
		{
			stringstream address;
			address.str("");
			address << box.info.str() << "-r-v.bin";
			out_file.open(address.str().c_str());
		}

		if (thisnode->node_id == 0)
			cout << " Box information is: " << box.info.str() << endl;

		MPI_Barrier(MPI_COMM_WORLD);
		t_eq = equilibrium(&box, equilibrium_step, saving_period, out_file);
		MPI_Barrier(MPI_COMM_WORLD);

		if (thisnode->node_id == 0)
			cout << " Done in " << (t_eq / 60.0) << " minutes" << endl;

		t_sim = data_gathering(&box, total_step, saving_period, out_file);
		MPI_Barrier(MPI_COMM_WORLD);

		if (thisnode->node_id == 0)
		{
			cout << " Done in " << (t_sim / 60.0) << " minutes" << endl;
			out_file.close();
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
}

void Change_Alpha(int argc, char *argv[], Node* thisnode)
{
	Real input_rho = atof(argv[1]);
	Real input_g = atof(argv[2]);
	Real input_noise = atof(argv[3]);

	vector<Real> alpha_list;
	for (int i = 4; i < argc; i++)
		alpha_list.push_back(atof(argv[i]));

	Real t_eq,t_sim;

	Box box;
	box.Init(thisnode, input_rho, input_g, alpha_list[0], 0);

	ofstream out_file;

	Particle::noise_amplitude = input_noise / sqrt(dt); // noise amplitude depends on the step (dt) because of ito calculation. If we have epsilon in our differential equation and we descritise it with time steps dt, the noise in each step that we add is epsilon times sqrt(dt) if we factorise it with a dt we have dt*(epsilon/sqrt(dt)).

	for (int i = 0; i < alpha_list.size(); i++)
	{
		Particle::alpha = alpha_list[i];
		box.info.str("");
		box.info << "rho=" << box.density <<  "-g=" << Particle::g << "-alpha=" << Particle::alpha << "-noise=" << input_noise << "-cooling";

		if (thisnode->node_id == 0)
		{
			stringstream address;
			address.str("");
			address << box.info.str() << "-r-v.bin";
			out_file.open(address.str().c_str());
		}

		if (thisnode->node_id == 0)
		{
			cout << " Box information is: " << box.info.str() << endl;
			Triangle_Lattice_Formation(box.particle, box.N, 1);
		}

		MPI_Barrier(MPI_COMM_WORLD);

		// Master node will broadcast the particles information
		thisnode->Root_Bcast();
		// Any node update cells, knowing particles and their cell that they are inside.
		thisnode->Full_Update_Cells();

		MPI_Barrier(MPI_COMM_WORLD);

		t_eq = equilibrium(&box, equilibrium_step, saving_period, out_file);
		MPI_Barrier(MPI_COMM_WORLD);

		if (thisnode->node_id == 0)
			cout << " Done in " << (t_eq / 60.0) << " minutes" << endl;

		t_sim = data_gathering(&box, total_step, saving_period, out_file);
		MPI_Barrier(MPI_COMM_WORLD);

		if (thisnode->node_id == 0)
		{
			cout << " Done in " << (t_sim / 60.0) << " minutes" << endl;
			out_file.close();
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
}

int main(int argc, char *argv[])
{
	int this_node_id, total_nodes;
	MPI_Status status;
	MPI_Init(&argc, &argv);





	Node thisnode;

	#ifdef COMPARE
		thisnode.seed = seed;
	#else
		thisnode.seed = time(NULL) + thisnode.node_id*112488;
		while (!thisnode.Chek_Seeds())
		{
			thisnode.seed = time(NULL) + thisnode.node_id*112488;
			MPI_Barrier(MPI_COMM_WORLD);
		}
	#endif

	C2DVector::Init_Rand(thisnode.seed);

//	Change_Noise(argc, argv, &thisnode);
	Change_Alpha(argc, argv, &thisnode);

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();

}

