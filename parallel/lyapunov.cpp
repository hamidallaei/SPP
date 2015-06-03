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
//		timing_information(box->thisnode,start_time,i,equilibrium_step);
	}

	if (box->thisnode->node_id == 0)
		cout << "Finished" << endl;

	end_time = clock();

	Real t = (Real) (end_time - start_time) / CLOCKS_PER_SEC;
	return(t);
}

inline Real Lyapunov_Exponent(Box* box, int number_of_repeat, Real interval, Real duration, ofstream& outfile)
{
	clock_t start_time, end_time;
	start_time = clock();

	#ifdef TRACK_PARTICLE
	if (!flag)
		flag = true;
	#endif

	if (box->thisnode->node_id == 0)
		cout << "Finding deviiations data:" << endl;

	outfile << "0" << "\t" << "1" << "\t" << "0" << endl;
	for (Real t = interval; t < duration; t+=interval)
	{
		if (box->thisnode->node_id == 0)
			cout << "tau is: " << t << endl;
		int i = (int) round(t/dt);
		Real error;
		Real mrc = box->Mean_Rlative_Change(i, number_of_repeat, error);
		outfile << i*dt << "\t" << mrc << "\t" << error << endl;
//		timing_information(box->thisnode,start_time,i,total_step);
		MPI_Barrier(MPI_COMM_WORLD);
	}
	
	if (box->thisnode->node_id == 0)
		cout << "Finished" << endl;

	outfile.close();
	
	end_time = clock();

	Real t = (Real) (end_time - start_time) / CLOCKS_PER_SEC;
	return(t);
}

void Run(int argc, char *argv[], Node* thisnode)
{
	Real input_rho = atof(argv[1]);
	Real input_mu_plus = atof(argv[2]);
	Real input_mu_minus = atof(argv[3]);
	Real input_D = atof(argv[4]);

	Real t_eq,t_sim;

	Box box;
	box.Init(thisnode, input_rho);

	MarkusParticle::mu_plus = input_mu_plus;
	MarkusParticle::mu_minus = input_mu_minus;

	Particle::D_phi = input_D;
	Particle::noise_amplitude = sqrt(2*input_D) / sqrt(dt); // noise amplitude depends on the step (dt) because of ito calculation. If we have epsilon in our differential equation and we descritise it with time steps dt, the noise in each step that we add is epsilon times sqrt(dt) if we factorise it with a dt we have dt*(epsilon/sqrt(dt)).
	box.info.str("");
	box.info << "rho=" << box.density <<  "-mu+=" << Particle::mu_plus << "-mu-=" << Particle::mu_minus << "-Dphi=" << Particle::D_phi << "-L=" << Lx;

	ofstream out_file,traj_file;
	if (thisnode->node_id == 0)
	{
		stringstream address;
		address.str("");
		address << box.info.str() << "-r-v.bin";
		traj_file.open(address.str().c_str());
		address.str("");
		address << box.info.str() << "-deviation.dat";
		out_file.open(address.str().c_str());
	}

	if (thisnode->node_id == 0)
		cout << " Box information is: " << box.info.str() << endl;

	MPI_Barrier(MPI_COMM_WORLD);
	t_eq = equilibrium(&box, equilibrium_step, saving_period, out_file);
	MPI_Barrier(MPI_COMM_WORLD);

	if (thisnode->node_id == 0)
		cout << " Done in " << floor(t_eq / 60.0) << " minutes and " << t_eq - 60*floor(t_eq / 60.0) << " s" << endl;

	MPI_Barrier(MPI_COMM_WORLD);
	t_sim = Lyapunov_Exponent(&box, 4, 0.5,20.0,out_file);
	MPI_Barrier(MPI_COMM_WORLD);

	traj_file << &box;
	if (thisnode->node_id == 0)
	{
		cout << " Done in " << floor(t_sim / 60.0) << " minutes and " << t_sim - 60*floor(t_sim / 60.0) << " s" << endl;
		out_file.close();
		traj_file.close();
	}

	MPI_Barrier(MPI_COMM_WORLD);
}

bool Run_From_File(int argc, char *argv[], Node* thisnode)
{
	Real t_eq,t_sim;

	Box box;
	if (!box.Init(thisnode, argv[1]))
		return false;

	Particle::noise_amplitude = sqrt(2*Particle::D_phi) / sqrt(dt); // noise amplitude depends on the step (dt) because of ito calculation. If we have epsilon in our differential equation and we descritise it with time steps dt, the noise in each step that we add is epsilon times sqrt(dt) if we factorise it with a dt we have dt*(epsilon/sqrt(dt)).
	box.info.str("");
	box.info << "rho=" << box.density <<  "-mu+=" << Particle::mu_plus << "-mu-=" << Particle::mu_minus << "-Dphi=" << Particle::D_phi << "-L=" << Lx;

	ofstream out_file,traj_file;
	if (thisnode->node_id == 0)
	{
		stringstream address;
		address.str("");
		address << box.info.str() << "-r-v.bin";
		traj_file.open(address.str().c_str());
		address.str("");
		address << box.info.str() << "-deviation.dat";
		out_file.open(address.str().c_str());
	}

	if (thisnode->node_id == 0)
		cout << " Box information is: " << box.info.str() << endl;

	MPI_Barrier(MPI_COMM_WORLD);
	t_sim = Lyapunov_Exponent(&box, 4, 0.5,20.0, out_file);
	MPI_Barrier(MPI_COMM_WORLD);

	if (thisnode->node_id == 0)
	{
		cout << " Done in " << floor(t_sim / 60.0) << " minutes and " << t_sim - 60*floor(t_sim / 60.0) << " s" << endl;
		out_file.close();
		traj_file.close();
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	return true;
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

	if (argc > 2)
		Run(argc, argv, &thisnode);
	else
		Run_From_File(argc, argv, &thisnode);

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();

}

