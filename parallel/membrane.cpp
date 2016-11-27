#include "../shared/parameters.h"
#include "../shared/c2dvector.h"
#include "../shared/particle.h"
#include "../shared/cell.h"
#include "beadbox.h"

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

inline Real data_gathering(Box* box, long int total_step, int saving_period, ofstream& out_file, ofstream& variables_file)
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

	box->Save_Membrane_Position();
	for (long int i = 0; i < total_step; i+=cell_update_period)
	{
		if ((i / cell_update_period) % saving_period == 0)
		{
			out_file << box;
			timing_information(box->thisnode,start_time,i,total_step);
		}
		box->Save_All_Variables(variables_file);
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
	if (argc < 3)
	{
		if (box.thisnode->node_id == 0)
			cout << "arguments are: \n" << "Number of active particles,\tNumber of membrane particles" << endl;
		exit(0);
	}
	string test = argv[1];
	int input_file = 0;

	if (test.length() > 10)
		input_file = 1;

	int input_chain_length = atoi(argv[1+input_file]);
	int input_Nm = atoi(argv[2+input_file]);
	Real input_packing_fraction = atof(argv[3+input_file]);
	Real input_chiral_radius = atof(argv[4+input_file]);

	box.packing_fraction = input_packing_fraction;
	int input_Ns = (int) round(input_packing_fraction*input_Nm*input_Nm / (input_chain_length*M_PI*M_PI));
//	box.packing_fraction = box.Ns*input_chain_length*M_PI*M_PI / (box.Nm*box.Nm);

	Real t_eq,t_sim;

	Particle::Dr = 0.1;
	Particle::noise_amplitude = sqrt(2*Particle::Dr*dt);
	Particle::R0 = input_chiral_radius;

//	Particle::lambda = 0.1; // tumbling rate
//	Particle::t_tumble = 0.1/Particle::lambda; // tumbling duration
//	Particle::torque_tumble = 0.2; // torque strength of a tumble




// The following must be before box.init
	for (int i = 0; i < input_Nm; i++)
		box.particle[i].Set_Parameters(1,0.0);
	for (int i = input_Nm; i < input_Ns+input_Nm; i++)
		box.particle[i].Set_Parameters(input_chain_length,1.0);

	box.Init(box.thisnode, input_Ns, input_Nm);


	if (input_file == 0)
	{
		if (box.thisnode->node_id == 0)
		{
			Ring_Membrane(box.particle, box.Nm);
			Confined_In_Ring_Membrane(box.particle, box.Ns, box.Nm);
//			Triangle_Lattice_Formation(&box.particle[box.Nm], box.Ns, 2*Particle::sigma_p);
		}
		box.Sync();
	}

	box.info.str("");
	box.info << "Nm=" << box.Nm;
	box.info << "-phi=" << box.packing_fraction;
	box.info << "-R0=" << Particle::R0;

	ofstream out_file;
	ofstream variables_file;

	if (box.thisnode->node_id == 0)
	{
		stringstream address;
		address.str("");
		address << "r-v-" << box.info.str() << ".bin";
		out_file.open(address.str().c_str());
		address.str("");
		address << "variables-" << box.info.str() << ".dat";
		variables_file.open(address.str().c_str());
	}

	if (box.thisnode->node_id == 0)
		cout << " Box information is: " << box.info.str() << endl;

		t_sim = data_gathering(&box, total_step, saving_period, out_file, variables_file);
		MPI_Barrier(MPI_COMM_WORLD);

		if (box.thisnode->node_id == 0)
		{
			cout << " Done in " << (t_sim / 60.0) << " minutes" << endl;
			out_file.close();
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

	Run(box, argc, argv);

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
}

