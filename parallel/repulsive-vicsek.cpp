#include "../shared/parameters.h"
#include "../shared/c2dvector.h"
#include "../shared/particle.h"
#include "../shared/cell.h"
#include "box.h"

inline void timing_information(Box* box, clock_t start_time, Real box_initial_time)
{
	if (box->thisnode->node_id == 0)
	{
		clock_t current_time = clock();
		int lapsed_time = (current_time - start_time) / (CLOCKS_PER_SEC);
		int remaining_time = (lapsed_time*(sim_time - box->t)) / (box->t -  box_initial_time + 1);
		cout << "\r" << box->t << "\t" << round(100.0*box->t / sim_time) << "% lapsed time: " << lapsed_time << " s		remaining time: " << remaining_time << " s" << flush;
	}
	MPI_Barrier(MPI_COMM_WORLD);
}

bool does_file_exist(const char *fileName)
{
	std::ifstream infile(fileName);
	bool result = infile.good();
	infile.close();
	return result;
}

inline Real equilibrium(Box* box, long int equilibrium_step, int saving_period, ofstream& out_file)
{
	clock_t start_time, end_time;
	start_time = clock();
	Real box_initial_time = box->t;

	if (box->thisnode->node_id == 0)
		cout << "equilibrium:" << endl;

	for (long int i = 0; i < equilibrium_step; i+=cell_update_period)
	{
		box->Multi_Step(cell_update_period);
		timing_information(box,start_time, box_initial_time);
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
	Real box_initial_time = box->t;

	#ifdef TRACK_PARTICLE
	if (!flag)
		flag = true;
	#endif

	if (box->thisnode->node_id == 0)
		cout << "gathering data:" << endl;
	int saving_time = 0;

	if (box->t < dt)
			out_file << box;

	Real initial_dt = dt;
	Set_time_step(dt);

	int i = 0;
	while (box->t < 1.0)
	{
		box->Multi_Step(cell_update_period);
		i += cell_update_period;
		if ((i / cell_update_period) % saving_period == 0)
		{
			out_file << box;
			timing_information(box,start_time, box_initial_time);
		}
	}

	Set_time_step(initial_dt);

	while (box->t < sim_time)
	{
		box->Multi_Step(cell_update_period);
		i += cell_update_period;
		if ((i / cell_update_period) % saving_period == 0)
		{
			out_file << box;
			timing_information(box,start_time, box_initial_time);
		}
	}

	if (box->thisnode->node_id == 0)
		cout << "Finished" << endl;

	end_time = clock();

	Real t = (Real) (end_time - start_time) / CLOCKS_PER_SEC;
	return(t);
}

void Single_Run(Box& box, int argc, char *argv[])
{
	if (argc < 4)
	{
		if (box.thisnode->node_id == 0)
			cout << "arguments are: \n" << "packingfraction,\tg,\tDr" << endl;
		exit(0);
	}
	Real input_packing_fraction = atof(argv[2]);
	Real input_A_p =atof(argv[3]);
	Real input_g = atof(argv[4]); 

	Real input_rho = 4*input_packing_fraction / (M_PI*Particle::sigma_p*Particle::sigma_p);

	Real input_noise = atof(argv[5]);

	bool FROM_FILE = false;
	string name;
	if (argc == 7)
	{
		name = argv[6];
		FROM_FILE = true;
	}

	Real t_eq,t_sim;

	Particle::Set_nb(1);
	Particle::Set_F0(1);
	Particle::Set_sigma_p(1);
	Particle::Set_repulsion_radius(1.05);
	Particle::Set_alignment_radius(1.1);
	Particle::Set_A_p(input_A_p);
	Particle::Set_g(1.0);

	box.Init(box.thisnode, input_rho);
	box.packing_fraction = input_packing_fraction;

	Particle::Set_Dr(input_noise); // This will set the noise amplitude as well. Noise amplitude depends on the step (dt) because of ito calculation. If we have epsilon in our differential equation and we descritise it with time steps dt, the noise in each step that we add is epsilon times sqrt(dt) if we factorise it with a dt we have dt*(epsilon/sqrt(dt)).

	box.info.str("");
	box.info << "phi=" << box.packing_fraction <<  "-g=" << Particle::g << "-Ap=" << RepulsiveParticle::A_p  << "-noise=" << input_noise << "-dt_inv=" << 1.0/dt << "-L=" << Lx << "-N=" << box.Ns;

	ofstream out_file;

	stringstream address;
	address.str("");
	address << box.info.str() << "-r-v.bin";

	if (FROM_FILE)
	{
		box.Positioning_Particles(name);
		if (box.thisnode->node_id == 0)
			out_file.open(address.str().c_str(), ios::out | ios::app | ios::binary);
	}
	else
	{
		if (does_file_exist(address.str().c_str()))
		{
			box.Positioning_Particles(address.str());
			if (box.thisnode->node_id == 0)
				out_file.open(address.str().c_str(), ios::out | ios::app | ios::binary);
		}
		else
		{
			if (box.thisnode->node_id == 0)
			{
				Triangle_Lattice_Formation(box.particle,box.Ns,0);
				for (int i = 0; i < box.Ns; i++)
					box.particle[i].r_original = box.particle[i].r;
				out_file.open(address.str().c_str());
			}
			box.Sync();
		}
	}

	if (box.thisnode->node_id == 0)
		cout << " Box information is: " << box.info.str() << endl;

		t_sim = data_gathering(&box, total_step, saving_period, out_file);
		MPI_Barrier(MPI_COMM_WORLD);

		if (box.thisnode->node_id == 0)
		{
			cout << " Done in " << (t_sim / 60.0) << " minutes" << endl;
			out_file.close();
		}
	MPI_Barrier(MPI_COMM_WORLD);
}

void Change_Noise(Box& box, int argc, char *argv[])
{
	if (argc < 5)
	{
		if (box.thisnode == 0)
			cout << "arguments are: \n" << "density,\tg,\tepsilon" << endl;
		exit(0);
	}
	Real input_packing_fraction = atof(argv[2]);
	Real input_g = atof(argv[3]);

	Real input_rho = 4*input_packing_fraction / (M_PI*Particle::sigma_p*Particle::sigma_p);

	vector<Real> noise_list;
	for (int i = 4; i < argc; i++)
		noise_list.push_back(atof(argv[i]));

	Real t_eq,t_sim;


	Particle::Set_nb(1);
	Particle::Set_F0(1);
	Particle::Set_sigma_p(0.95);
	Particle::Set_repulsion_radius(1.0);
	Particle::Set_alignment_radius(1.5);
	Particle::Set_A_p(200.0);
	Particle::Set_g(input_g);

	box.Init(box.thisnode, input_rho);
	box.packing_fraction = input_packing_fraction;

	ofstream out_file;

	if (box.thisnode->node_id == 0)
	{
		Triangle_Lattice_Formation(box.particle,box.Ns,0);
		for (int i = 0; i < box.Ns; i++)
			box.particle[i].r_original = box.particle[i].r;
	}
	box.Sync();
	for (int i = 0; i < box.Ns; i++)
		box.particle[i].r_original = box.particle[i].r;

	for (int i = 0; i < noise_list.size(); i++)
	{
		Particle::Set_Dr(noise_list[i]); // This will set the noise amplitude as well. Noise amplitude depends on the step (dt) because of ito calculation. If we have epsilon in our differential equation and we descritise it with time steps dt, the noise in each step that we add is epsilon times sqrt(dt) if we factorise it with a dt we have dt*(epsilon/sqrt(dt)).

		box.info.str("");
		box.info << "phi=" << box.packing_fraction <<  "-g=" << Particle::g << "-Ap=" << RepulsiveParticle::A_p  << "-noise=" << noise_list[i] << "-dt_inv=" << 1.0/dt << "-L=" << Lx << "-N=" << box.Ns;

		if (box.thisnode->node_id == 0)
		{
			stringstream address;
			address.str("");
			address << box.info.str() << "-r-v.bin";
			out_file.open(address.str().c_str());
		}

		if (box.thisnode->node_id == 0)
			cout << " Box information is: " << box.info.str() << endl;

//		MPI_Barrier(MPI_COMM_WORLD);
//		t_eq = equilibrium(&box, equilibrium_step, saving_period, out_file);
//		MPI_Barrier(MPI_COMM_WORLD);

//		if (box.thisnode->node_id == 0)
//			cout << " Done in " << (t_eq / 60.0) << " minutes" << endl;

		t_sim = data_gathering(&box, total_step, saving_period, out_file);
		MPI_Barrier(MPI_COMM_WORLD);

		if (box.thisnode->node_id == 0)
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

	Real input_Lx = atof(argv[1]);
	Real input_packing_fraction = atof(argv[2]);
	Real input_g = atof(argv[3]);

	Box box(input_Lx, input_Lx, input_packing_fraction);

	Node thisnode;

	#ifdef COMPARE
		thisnode.Init_Rand(seed);
	#else
		thisnode.Init_Rand();
	#endif

	box.thisnode = &thisnode;

	//Change_Noise(box, argc, argv);
	Single_Run(box, argc, argv);

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
}

