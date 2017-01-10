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
		cout << "\r" << round(100.0*i_step / total_step) << "% lapsed time: " << lapsed_time << " s		remaining time: " << remaining_time << " s" << "\t P=" << node->polarization << flush;
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


inline Real data_gathering(Box* box, long int total_step, int saving_period, ofstream& out_file, ofstream& polarization_file)
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
		if ((i / cell_update_period) % saving_period == 0)
		{
			out_file << box;
			timing_information(box->thisnode,start_time,i,total_step);
		}
		box->Save_Polarization(polarization_file);
		box->Multi_Step(cell_update_period);
	}
	if ((total_step / cell_update_period) % saving_period == 0)
			out_file << box;

	if (box->thisnode->node_id == 0)
		cout << "Finished" << endl;

	end_time = clock();

	Real t = (Real) (end_time - start_time) / CLOCKS_PER_SEC;
	return(t);
}





void Run(Box& box, int argc, char *argv[])
{
	if (argc < 4)
	{
		if (box.thisnode->node_id == 0)
			cout << "arguments are: \n" << "density,\tg,\talpha,\tDr" << endl;
		exit(0);
	}
	string test = argv[1];
	int input_file = 0;

	if (test.length() > 10)
		input_file = 1;

	Real input_rho = atof(argv[1+input_file]);
	Real input_g = atof(argv[2+input_file]);
	Real input_g_phi = atof(argv[3+input_file]);
	Real input_alpha = atof(argv[4+input_file]);
	Real input_Rphi = atof(argv[5+input_file]);
	Real input_D_phi = atof(argv[6+input_file]);
	Real input_Dr = atof(argv[7+input_file]);
	Real input_K = atof(argv[8+input_file]);
	Real input_omega = atof(argv[9+input_file]);
	Real input_vmin = atof(argv[10+input_file]);
	Real input_vmax = atof(argv[11+input_file]);


	Real t_eq,t_sim;

	if (input_file == 0)
	{
		box.Init(box.thisnode, input_rho);
		if (box.thisnode->node_id == 0)
		{
// Positioning the particles
				Random_Formation(box.particle, box.N, 0.5); // Positioning partilces Randomly, but distant from 
			for (int i = 0; i < box.N; i++)
				box.particle[i].phi = gsl_ran_flat(C2DVector::gsl_r,-M_PI,M_PI);
		}
		box.Sync();
	}
	else
	{
		bool b = box.Positioning_Particles(box.thisnode, test);
		if (box.density != input_rho)
			cout << "Confilict of density!" << endl;
		if (!b)
			exit(0);
	}
	Particle::g = input_g;
	Particle::g_phi = input_g_phi;
	Particle::alpha = input_alpha;
	Particle::Rphi = input_Rphi;
	Particle::D_phi = input_D_phi;
	Particle::Dr = input_Dr;
	Particle::K = input_K;
	Particle::omega = input_omega;
	Particle::vmin = input_vmin;
	Particle::vmax = input_vmax;

	Particle::Set_Variables();

	ofstream out_file, polarization_file;

	box.info.str("");
	box.info << "rho=" << box.density <<  "-g=" << Particle::g << "-alpha=" << Particle::alpha << "-vmin=" << Particle::vmin << "-vmax=" << Particle::vmax << "-Dr=" << Particle::Dr << "-K=" << Particle::K << "-2Lx=" << Lx2 << "-2Ly=" << Ly2;

	if (box.thisnode->node_id == 0)
	{
		stringstream address;
		address.str("");
		address << box.info.str() << "-r-v.bin";
		out_file.open(address.str().c_str());
		address.str("");
		address << "polarization-time-" << box.info.str() << ".dat";
		polarization_file.open(address.str().c_str());
	}

	if (box.thisnode->node_id == 0)
		cout << " Box information is: " << box.info.str() << endl;

	MPI_Barrier(MPI_COMM_WORLD);
	t_eq = equilibrium(&box, equilibrium_step, saving_period, out_file);
	MPI_Barrier(MPI_COMM_WORLD);

	if (box.thisnode->node_id == 0)
		cout << " Done in " << (t_eq / 60.0) << " minutes" << endl;

	t_sim = data_gathering(&box, total_step, saving_period, out_file, polarization_file);
	MPI_Barrier(MPI_COMM_WORLD);

	if (box.thisnode->node_id == 0)
	{
		cout << " Done in " << (t_sim / 60.0) << " minutes" << endl;
		out_file.close();
		polarization_file.close();
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
		thisnode.Init_Rand(seed);
	#else
		thisnode.Init_Rand();
	#endif

	Box box;
	box.thisnode = &thisnode;

	Run(box, argc, argv);

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
}

