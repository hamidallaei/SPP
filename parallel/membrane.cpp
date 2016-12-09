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

inline Real data_gathering(Box* box, long int total_step, int trajectory_saving_period, int quantities_saving_period, ofstream& out_file, ofstream& variables_file)
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
		if ((i / cell_update_period) % trajectory_saving_period == 0)
		{
			out_file << box;
			timing_information(box->thisnode,start_time,i,total_step);
		}
		box->Multi_Step(cell_update_period);
		if ((i / cell_update_period) % quantities_saving_period == 0)
			box->Save_All_Variables(variables_file);
	}
	if ((total_step / cell_update_period) % trajectory_saving_period == 0)
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

	if (box.thisnode->node_id == 0)
	{
		cout << "Arguments are:\t" << endl;
		cout << "1-membrane elasticity" << endl;
		cout << "2-membrane radius" << endl;
		cout << "3-chirality radius" << endl;
		cout << "4-packing fraction" << endl;
	}

	int input_chain_length = 2;
	Real input_membrane_elasticity = atof(argv[2+input_file]);
	Real input_membrane_radius = atof(argv[3+input_file]);
	Real input_chiral_radius = atof(argv[4+input_file]);
	Real input_packing_fraction = atof(argv[5+input_file]);

//	int input_Nm = (int) round(2*M_PI*input_membrane_radius / Particle::sigma_p);
//	int input_Ns = (int) round(input_packing_fraction*input_Nm*input_Nm / (input_chain_length*M_PI*M_PI));
	int input_Nm = (int) round(M_PI/asin(0.5*Particle::sigma_p/input_membrane_radius));
	int input_Ns = (int) round(input_packing_fraction/( input_chain_length*sin(M_PI/input_Nm)*sin(M_PI/input_Nm) ));

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
	box.membrane_elasticity = input_membrane_elasticity;
	box.membrane_radius = input_membrane_radius;
	box.packing_fraction = input_packing_fraction;

	if (input_file == 0)
	{
		if (box.thisnode->node_id == 0)
		{
			Ring_Membrane(box.particle, box.Nm);
			Confined_In_Ring_Membrane(box.particle, box.Ns, box.Nm);
		}
		box.Sync();
	}

	box.info.str("");
	box.info << "k="<< box.membrane_elasticity;
	box.info << "-R="<< box.membrane_radius;
	box.info << "-r=" << Particle::R0;
	box.info << "-phi=" << box.packing_fraction;
	box.info << "-ABP";

	ofstream out_file;
	ofstream variables_file;

	if (box.thisnode->node_id == 0)
	{
		stringstream address;
		address.str("");
		address << "r-v-" << box.info.str() << ".bin";
		out_file.open(address.str().c_str());
		address.str("");
		address << "quantities-" << box.info.str() << ".dat";
		variables_file.open(address.str().c_str());
	}

	if (box.thisnode->node_id == 0)
		cout << " Box information is: " << box.info.str() << endl;

		int quantities_saving_period = ( (int) round(1/dt) ) / cell_update_period;
		t_sim = data_gathering(&box, total_step, saving_period, quantities_saving_period, out_file, variables_file);
		MPI_Barrier(MPI_COMM_WORLD);

		if (box.thisnode->node_id == 0)
		{
			cout << " Done in " << (t_sim / 60.0) << " minutes" << endl;
			out_file.close();
		}

	MPI_Barrier(MPI_COMM_WORLD);
}



void Init_Nodes(Node& thisnode, const int& seed)
{
//	thisnode.seed = time(NULL) + thisnode.node_id*112488;
//	while (!thisnode.Chek_Seeds())
//	{
//		thisnode.seed = time(NULL) + thisnode.node_id*112488;
//		MPI_Barrier(MPI_COMM_WORLD);
//	}
	thisnode.seed = seed +  thisnode.node_id*112488;
	C2DVector::Init_Rand(thisnode.seed);
	MPI_Barrier(MPI_COMM_WORLD);
}

int main(int argc, char *argv[])
{
	int this_node_id, total_nodes;
	MPI_Status status;
	MPI_Init(&argc, &argv);

	Node thisnode;
	input_seed = atoi(argv[1]);

	Init_Nodes(thisnode, input_seed);

	Box box;
	box.thisnode = &thisnode;

	Run(box, argc, argv);

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
}

