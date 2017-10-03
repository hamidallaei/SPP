#include "../shared/parameters.h"
#include "../shared/c2dvector.h"
#include "../shared/particle.h"
#include "../shared/cell.h"
#include "beadbox.h"

int input_seed;

inline void timing_information(Box* box, clock_t start_time, Real box_initial_time)
{
		clock_t current_time = clock();
		int lapsed_time = (current_time - start_time) / (CLOCKS_PER_SEC);
		int remaining_time = (lapsed_time*(sim_time - box->t)) / (box->t -  box_initial_time + 1);
		cout << "\r" << round(100.0*box->t / sim_time) << "% lapsed time: " << lapsed_time << " s		remaining time: " << remaining_time << " s" << flush;
}

bool does_file_exist(const char *fileName)
{
	std::ifstream infile(fileName);
	bool result = infile.good();
	infile.close();
	return result;
}

inline Real data_gathering(Box* box, long int total_step, int trajectory_saving_period, int quantities_saving_period, ofstream& out_file, ofstream& variables_file)
{
	clock_t start_time, end_time;
	start_time = clock();
	Real box_initial_time = box->t;

	#ifdef TRACK_PARTICLE
	if (!flag)
		flag = true;
	#endif

	cout << "gathering data:" << endl;
	int saving_time = 0;

	if (box->t < dt)
	{
		box->Save_Particles_Positions();
		out_file << box;
		timing_information(box, start_time, box_initial_time);
	}

	int i = 0;
	while (box->t < sim_time)
	{
		box->Multi_Step(cell_update_period);
		i += cell_update_period;
		if ((i / cell_update_period) % trajectory_saving_period == 0)
		{
			out_file << box;
			timing_information(box, start_time, box_initial_time);
		}
		if ((i / cell_update_period) % quantities_saving_period == 0)
			box->Save_All_Variables(variables_file);
	}
	if ((total_step / cell_update_period) % trajectory_saving_period == 0)
		out_file << box;

	cout << "Finished" << endl;

	end_time = clock();

	Real t = (Real) (end_time - start_time) / CLOCKS_PER_SEC;
	return(t);
}

Real input_membrane_elasticity;
Real input_membrane_radius;
Real input_chiral_radius;
Real input_packing_fraction;
Real input_cover_fraction;
int input_file = 0;
int input_chain_length = 3;

void Read_Arguments(int argc, char *argv[])
{
	if (argc < 5)
	{
			cout << "arguments are: \n" << "Number of active particles,\tNumber of membrane particles" << endl;
		exit(0);
	}
	string test = argv[1];
	int input_file = 0;

	if (test.length() > 10)
		input_file = 1;

	cout << "Arguments are:\t" << endl;
	cout << "1-membrane elasticity" << endl;
	cout << "2-membrane radius" << endl;
	cout << "3-chirality radius" << endl;
	cout << "4-cover fraction" << endl;

	input_membrane_elasticity = atof(argv[1+input_file]);
	input_membrane_radius = atof(argv[2+input_file]);
	input_chiral_radius = atof(argv[3+input_file]);
	input_cover_fraction = atof(argv[4+input_file]);

}

void Run(Box& box, int argc, char *argv[])
{
	Real membrane_l_eq = Particle::sigma_m + 13*Particle::A_p*(Particle::repulsion_radius_m - Particle::sigma_m)/(input_membrane_elasticity*Particle::sigma_m + 13*Particle::A_p);
	membrane_l_eq = Particle::repulsion_radius_m;
	int input_Nm = (int) round(M_PI/asin(0.5*membrane_l_eq/input_membrane_radius));
//	int input_Ns = (int) round(input_packing_fraction/( input_chain_length*sin(M_PI/input_Nm)*sin(M_PI/input_Nm) ));
	int input_Ns = (int) round(M_PI * input_cover_fraction * Particle::sigma_p * (1 / sin(M_PI/input_Nm/Particle::repulsion_radius_m)  - 2*Particle::sigma_p ));

	bool FROM_FILE = false;
	string name;
	if (argc == 7)
	{
		name = argv[6];
		FROM_FILE = true;
	}

	Real t_eq,t_sim;

//	Particle::lambda = 0.1; // tumbling rate
//	Particle::t_tumble = 0.1/Particle::lambda; // tumbling duration
//	Particle::torque_tumble = 0.2; // torque strength of a tumble


	Particle::Set_Dr(1.0/8);
	Particle::Set_separation(1.0/(input_chain_length-1));

// The following must be before box.init
	for (int i = 0; i < input_Nm; i++)
	{
		box.particle[i].Set_nb(1);
		box.particle[i].Set_F0(0.0);
		box.particle[i].Set_R0(input_chiral_radius);
	}
	for (int i = input_Nm; i < input_Ns+input_Nm; i++)
	{
		box.particle[i].Set_nb(input_chain_length);
		box.particle[i].Set_F0(1.0);
		box.particle[i].Set_R0(input_chiral_radius);
//		box.particle[i].Set_Parameters(input_chain_length,1.0);
	}

	box.Init(input_Ns, input_Nm);
	box.membrane_elasticity = input_membrane_elasticity;
	box.membrane_radius = input_membrane_radius;
	box.packing_fraction = input_packing_fraction;
	box.cover_fraction = input_cover_fraction;

	box.info.str("");
	box.info << "k="<< box.membrane_elasticity;
	box.info << "-R="<< box.membrane_radius;
	box.info << "-r=" << Particle::R0;
	box.info << "-psi=" << box.cover_fraction;
	box.info << "-Dr=" << Particle::Dr;
	box.info << "-seed=" << input_seed;
	box.info << "-ABP";

	ofstream out_file;
	ofstream variables_file;

	stringstream traj_address;
	traj_address.str("");
	traj_address << "r-v-" << box.info.str() << ".bin";


	stringstream quant_address;
	quant_address.str("");
	quant_address << "quantities-" << box.info.str() << ".dat";


	if (FROM_FILE)
	{
		box.Positioning_Particles(traj_address.str());
		out_file.open(traj_address.str().c_str(), ios::out | ios::app | ios::binary);
		variables_file.open(quant_address.str().c_str(), ios::out | ios::app | ios::binary);
	}
	else
	{
		if (does_file_exist(traj_address.str().c_str()))
		{
			box.Positioning_Particles(traj_address.str());
			out_file.open(traj_address.str().c_str(), ios::out | ios::app | ios::binary);
			variables_file.open(quant_address.str().c_str(), ios::out | ios::app | ios::binary);
		}
		else
		{
			Ring_Membrane(box.particle, Particle::repulsion_radius_m, box.Nm);
			Confined_In_Ring_Membrane(box.particle, Particle::repulsion_radius_m, box.Ns, box.Nm);
			for (int i = 0; i < box.Ns; i++)
				box.particle[i].r_original = box.particle[i].r;
			out_file.open(traj_address.str().c_str());
			variables_file.open(quant_address.str().c_str());
		}
	}

	box.Sync();

	cout << " Box information is: " << box.info.str() << endl;

	cout << box.particle[184].m_parallel << "\t" << box.particle[184].m_perpendicular << "\t" << box.particle[184].k_perpendicular << endl;
	cout << box.particle[184].torque0 << "\t" << box.particle[184].F0 << "\t" << box.particle[184].R0 << endl;

	int quantities_saving_period = ( (int) round(1.0/dt) ) / cell_update_period;
	t_sim = data_gathering(&box, total_step, saving_period, quantities_saving_period, out_file, variables_file);

	cout << " Done in " << (t_sim / 60.0) << " minutes" << endl;
	out_file.close();
}

int main(int argc, char *argv[])
{
	input_seed = atoi(argv[5]);
	C2DVector::Init_Rand(input_seed);
	Read_Arguments(argc, argv);
	Box box(input_membrane_radius);
	Run(box, argc, argv);
}

