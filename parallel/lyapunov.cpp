#include "../shared/parameters.h"
#include "../shared/c2dvector.h"
#include "../shared/particle.h"
#include "../shared/cell.h"
#include "../shared/vector-set.h"
#include "box.h"

inline void timing_information(Node* node, clock_t start_time, int i_step, int total_step)
{
	if (node->node_id == 0)
	{
		clock_t current_time = clock();
		int lapsed_time = (current_time - start_time) / (CLOCKS_PER_SEC);
		int remaining_time = (lapsed_time*(total_step - i_step)) / (i_step + 1);
		if (i_step % 1000 == 0)
			cout << "\r" << round(100.0*i_step / total_step) << "% lapsed time: " << lapsed_time << " s		remaining time: " << remaining_time << " s" << flush;
	}
	MPI_Barrier(MPI_COMM_WORLD);
}

inline Real equilibrium(Box* box, long int equilibrium_step, int saving_period)
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

void Evolution(Box& box, VectorSet& us, int tau)
{
	VectorSet vs0(us);
	VectorSet vs(us.direction_num, box.N);
	
	static State_Hyper_Vector gamma(box.N);
	static State_Hyper_Vector gamma_prime(box.N);
	static State_Hyper_Vector gamma_0(box.N);

	box.Save(gamma_0);
	MPI_Barrier(MPI_COMM_WORLD);

	box.Multi_Step(tau, 20);

	box.Save(gamma);

	for (int i = 0; i < us.direction_num; i++)
	{
		box.Load(gamma_0);
		MPI_Barrier(MPI_COMM_WORLD);
		box.Add_Deviation(vs0.v[i]);
		box.Multi_Step(tau, 20);
		box.Save(gamma_prime);
		vs.v[i] = (gamma_prime - gamma);
	}

	vs.Renormalize(us);
	for (int i = 0; i < us.direction_num; i++)
	{
		us.v[i].growth = 0;
//		for (int j = 0; j < number_of_directions; j++)
//			us.v[i].growth += (vs.v[j] * us.v[i]) / (vs0.v[j] * us.v[i]);
		us.v[i].growth = (vs.v[i] * us.v[i]) / (vs0.v[i] * us.v[i]);
//		us.v[i].growth = vs.v[i].Magnitude() / vs0.v[i].Magnitude();
	}

	us.Scale();
	
	box.Load(gamma);
}

void Init_Deviations(Box& box, VectorSet& us, int tau, int iteration)
{
	us.Rand();
	MPI_Barrier(MPI_COMM_WORLD);
	
	for (int i = 0; i < iteration; i++)
		Evolution(box, us, tau);

	if (box.thisnode->node_id == 0)
	for (int i = 0; i < us.direction_num; i++)
		cout << us.v[i].particle[0].r / us.amplitude << "\t" << us.v[i].particle[0].theta / us.amplitude << endl;
}

// void Evolution_Tracking(Box& box, VectorSet& us, int tau, int interval, ofstream& outfile)
// {
// 	VectorSet vs0(us);
// 	VectorSet vs(us.N, N);
// 
// 	static State_Hyper_Vector gamma(N);
// 	static State_Hyper_Vector gamma_prime(N);
// 	static State_Hyper_Vector gamma_0(N);
// 
// 	for (int t = interval; t < tau; t += interval)
// 	{
// 		box->Save(gamma_0);
// 		MPI_Barrier(MPI_COMM_WORLD);
// 
// 		Multi_Step(interval, 20);
// 
// 		box->Save(gamma);
// 
// 		
// 		for (int i = 0; i < number_of_directions; i++)
// 		{
// 			box->Load(gamma_0);
// 			MPI_Barrier(MPI_COMM_WORLD);
// 			box->Add_Deviation(vs0.v[i]);
// 			Multi_Step(interval, 20);
// 			box->Save(gamma_prime);
// 			vs.v[i] = (gamma_prime - gamma);
// 		}
// 
// 		vs.Renormalize(us);
// 		vs0 = us;
// 		for (int i = 0; i < number_of_directions; i++)
// 		{
// 			us.v[i].growth = 0;
// 			us.v[i].growth = (vs.v[i] * us.v[i]) / (vs0.v[i] * us.v[i]);
// 		}
// 	}
// 
// 	box->Load(gamma);
// }

inline Real Lyapunov_Exponent(Box& box, int spectrum_levels, Real interval, Real duration, ofstream& outfile)
{
	clock_t start_time, end_time;
	start_time = clock();

	VectorSet::amplitude = 1e-4;

	#ifdef TRACK_PARTICLE
	if (!flag)
		flag = true;
	#endif

	outfile << 0;
	for (int i = 0; i < spectrum_levels; i++)
		outfile << "\t" << 1;
	outfile << endl;

	VectorSet us(spectrum_levels, box.N);

	Init_Deviations(box, us, round(1.0/dt), 1000);
	MPI_Barrier(MPI_COMM_WORLD);
	
	for (Real t = interval; t < duration; t+=interval)
	{
		if (box.thisnode->node_id == 0)
			cout << "tau is: " << t << endl;
		int i = (int) round(t/dt);
		Evolution(box, us, i);

		outfile << i*dt << "\t";
		for (int j = 0; j < spectrum_levels; j++)
			outfile << us.v[j].growth << "\t";
		outfile << endl;
//		timing_information(box.thisnode,start_time,i,total_step);
		MPI_Barrier(MPI_COMM_WORLD);
	}

	if (box.thisnode->node_id == 0)
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
		address << "deviation-" << box.info.str() << ".dat";
		out_file.open(address.str().c_str());
	}

	if (thisnode->node_id == 0)
		cout << " Box information is: " << box.info.str() << endl;

	MPI_Barrier(MPI_COMM_WORLD);
	t_eq = equilibrium(&box, equilibrium_step, saving_period);
	MPI_Barrier(MPI_COMM_WORLD);

	if (thisnode->node_id == 0)
		cout << " Done in " << floor(t_eq / 60.0) << " minutes and " << t_eq - 60*floor(t_eq / 60.0) << " s" << endl;

	MPI_Barrier(MPI_COMM_WORLD);
 	t_sim = Lyapunov_Exponent(box, 3, 1,10, out_file);
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

	ofstream out_file;
	if (thisnode->node_id == 0)
	{
		stringstream address;
		address.str("");
		address << "deviation-" << box.info.str() << ".dat";
		out_file.open(address.str().c_str());
	}

	if (thisnode->node_id == 0)
		cout << " Box information is: " << box.info.str() << endl;

	MPI_Barrier(MPI_COMM_WORLD);
	t_eq = equilibrium(&box, equilibrium_step, saving_period);
	MPI_Barrier(MPI_COMM_WORLD);

	if (thisnode->node_id == 0)
		cout << " Done in " << floor(t_eq / 60.0) << " minutes and " << t_eq - 60*floor(t_eq / 60.0) << " s" << endl;

	MPI_Barrier(MPI_COMM_WORLD);
	t_sim = Lyapunov_Exponent(box, 10, 0.01,10, out_file);
	MPI_Barrier(MPI_COMM_WORLD);

	if (thisnode->node_id == 0)
	{
		cout << " Done in " << floor(t_sim / 60.0) << " minutes and " << t_sim - 60*floor(t_sim / 60.0) << " s" << endl;
		out_file.close();
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	return true;
}

void Init_Nodes(Node& thisnode, int input_seed = seed)
{
	#ifdef COMPARE
		thisnode.seed = input_seed;
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

