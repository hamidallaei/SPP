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

// Intialize the box from a file, this includes reading particles information, updating cells and sending information to all nodes. 
bool Init_Box_From_File(Box& box, const string input_name)
{
	Real input_density;
	Real input_g;
	Real input_alpha;
	Real input_v;
	Real input_noise;
	Real input_K;
	Real input_Lx;
	Real input_Ly;

	string name = input_name;

	boost::replace_all(name, "-r-v.bin", "");
	boost::replace_all(name, "rho=", "");
	boost::replace_all(name, "-g=", "\t");
	boost::replace_all(name, "-alpha=", "\t");
	boost::replace_all(name, "-v=", "\t");
	boost::replace_all(name, "-noise=", "\t");
	boost::replace_all(name, "-Dr=", "\t");
	boost::replace_all(name, "-K=", "\t");
	boost::replace_all(name, "-2Lx=", "\t");
	boost::replace_all(name, "-2Ly=", "\t");

	stringstream ss_name(name);
	ss_name >> input_density;
	ss_name >> input_g;
	ss_name >> input_alpha;
	ss_name >> input_v;
	ss_name >> input_noise;
	ss_name >> input_K;
	ss_name >> input_Lx;
	ss_name >> input_Ly;
	input_Lx /= 2.0;
	input_Ly /= 2.0;
	if (input_Lx != Lx_int)
	{
		cout << "The specified box size " << input_Lx << " is not the same as the size in binary file which is " << Lx_int << " please recompile the code with the right Lx_int in parameters.h file." << endl;
		return false;
	}

	if (input_Ly != Ly_int)
	{
		cout << "The specified box size " << input_Ly << " is not the same as the size in binary file which is " << Ly_int << " please recompile the code with the right Lx_int in parameters.h file." << endl;
		return false;
	}

	box.density = input_density;
	box.N = (int) round(Lx2*Ly2*input_density);

	Particle::g = input_g;
	Particle::alpha = input_alpha;
	Particle::speed = input_v;
	Particle::Dr = input_noise;
	Particle::K = input_K;

	box.Positioning_Particles(box.thisnode, input_name);

	return (true);
}



void Change_Noise(Box& box, int argc, char *argv[])
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
	Real input_alpha = atof(argv[3+input_file]);

	vector<Real> noise_list;
	for (int i = 4+input_file; i < argc; i++)
		noise_list.push_back(atof(argv[i]));

	Real t_eq,t_sim;

	Particle::noise_amplitude = 0;
	Particle::Dr = noise_list[0];
	Particle::K = K;
	Particle::Kamp = sqrt(2*Particle::K*dt);

	if (input_file == 0)
	{
		box.Init(box.thisnode, input_rho);
		if (box.thisnode->node_id == 0)
		{
// Positioning the particles
			Dc = input_g * input_rho * (1-input_alpha) / 2.0;
			if (Particle::Dr > Dc)
				Random_Formation(box.particle, box.N, 0.5); // Positioning partilces Randomly, but distant from walls (the last argument is the distance from walls)
			else
				Polar_Formation(box.particle,box.N);
//			Triangle_Lattice_Formation(box.particle, box.N, 1);
//			Polar_Formation(box.particle,box.N);
//			Random_Formation_Circle(box.particle, box.N, Lx-1); // Positioning partilces Randomly, but distant from walls
//			Single_Vortex_Formation(box.particle, box.N);
//			Four_Vortex_Formation(box.particle, box.N);
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
	Particle::alpha = input_alpha;

	ofstream out_file, polarization_file;

	for (int i = 0; i < noise_list.size(); i++)
	{
		Particle::Dr = noise_list[i];
		Particle::noise_amplitude = sqrt(2*noise_list[i]) / sqrt(dt); // noise amplitude depends on the step (dt) because of ito calculation. If we have epsilon in our differential equation and we descritise it with time steps dt, the noise in each step that we add is epsilon times sqrt(dt) if we factorise it with a dt we have dt*(epsilon/sqrt(dt)).
		box.info.str("");
		box.info << "rho=" << box.density <<  "-g=" << Particle::g << "-alpha=" << Particle::alpha << "-v=" << Particle::speed << "-Dr=" << noise_list[i] << "-K=" << Particle::K << "-2Lx=" << Lx2 << "-2Ly=" << Ly2;

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
	}
	MPI_Barrier(MPI_COMM_WORLD);
}



void Init_Nodes(Node& thisnode)
{
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
	MPI_Barrier(MPI_COMM_WORLD);
}

int main(int argc, char *argv[])
{
	int this_node_id, total_nodes;
	MPI_Status status;
	MPI_Init(&argc, &argv);

	Node thisnode;
	Init_Nodes(thisnode);

	Box box;
	box.thisnode = &thisnode;
	Particle::speed = speed;

	Change_Noise(box, argc, argv);

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
}

