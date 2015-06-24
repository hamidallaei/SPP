#include "../shared/parameters.h"
#include "../shared/c2dvector.h"
#include "../shared/particle.h"
#include "../shared/cell.h"
#include "../shared/set-up.h"
#include "box.h"

inline void timing_information(clock_t start_time, int i_step, int total_step)
{
		clock_t current_time = clock();
		int lapsed_time = (current_time - start_time) / (CLOCKS_PER_SEC);
		int remaining_time = (lapsed_time*(total_step - i_step)) / (i_step + 1);
		cout << "\r" << round(100.0*i_step / total_step) << "% lapsed time: " << lapsed_time << " s		remaining time: " << remaining_time << " s" << flush;
}


inline void equilibrium(Box* box, int equilibrium_step, int saving_period, ofstream& out_file)
{
	clock_t start_time = clock();
	cout << "equilibrium:" << endl;
//	box->Init();
	for (int i = 0; i < equilibrium_step; i+=cell_update_period)
	{
		box->Multi_Step(cell_update_period);
		timing_information(start_time,i,equilibrium_step);
	}
	cout << "Finished" << endl;
}


inline void data_gathering(Box* box, int total_step, int saving_period, ofstream& out_file)
{
	clock_t start_time = clock();

	#ifdef TRACK_PARTICLE
		if (!flag)
			flag = true;
	#endif

	cout << "gathering data:" << endl;
	int saving_time = 0;
	for (int i = 0; i < total_step; i+=cell_update_period)
	{
		box->Multi_Step(cell_update_period);
		timing_information(start_time,i,total_step);

		if ((i / cell_update_period) % saving_period == 0)
			out_file << box;
	}

	cout << "Finished" << endl;
}

// Initialize the wall positions and numbers.
void Init_Topology(Box& box)
{
	box.geometry.Reset();
	#ifndef PERIODIC_BOUNDARY_CONDITION
		box.geometry.Add_Wall(-Lx, -Ly, -Lx, Ly);
		box.geometry.Add_Wall(-Lx, Ly, Lx, Ly);
		box.geometry.Add_Wall(Lx, Ly, Lx, -Ly);
		box.geometry.Add_Wall(Lx, -Ly, -Lx, -Ly);
	#endif
}

void Init(Box* box, Real input_density, Real input_mu_plus, Real input_mu_minus, Real input_Dphi)
{
	box->density = input_density;
	box->N = (int) round(Lx2*Ly2*box->density);
	cout << "number_of_particles = " << box->N << endl;

	Particle::mu_plus = input_mu_plus;
	Particle::mu_minus = input_mu_minus;
	Particle::D_phi = input_Dphi;
	Particle::noise_amplitude = sqrt(2*Particle::D_phi) / sqrt(dt);

// Positioning the particles
//	Polar_Formation(box->particle,box->N);
//	Triangle_Lattice_Formation(box->particle, box->N, 1);
	Random_Formation(box->particle, box->N, 0); // Positioning partilces Randomly, but distant from walls (the last argument is the distance from walls)
//	Random_Formation_Circle(box->particle, box->N, Lx-1); // Positioning partilces Randomly, but distant from walls
//	Single_Vortex_Formation(box->particle, box->N);
//	Four_Vortex_Formation(box->particle, box->N);

	#ifndef PERIODIC_BOUNDARY_CONDITION
		box->geometry.Add_Wall(Lx, Ly, Lx, -Ly);
		box->geometry.Add_Wall(Lx, -Ly, -Lx, -Ly);
		box->geometry.Add_Wall(-Lx, -Ly, -Lx, Ly);
		box->geometry.Add_Wall(-Lx, Ly, Lx, Ly);
	#endif

	box->Update_Cells();

	box->info.str("");
	box->info << "rho=" << box->density <<  "-mu+=" << Particle::mu_plus << "-mu-=" << Particle::mu_minus << "-Dphi=" << Particle::D_phi << "-L=" << Lx;
}

// Intialize the box from a file, this includes reading particles information, updating cells and sending information to all nodes. Unfortunately this is only for Markus partiles
bool Init(Box& box, const string input_name)
{
	#ifdef TRACK_PARTICLE
		track_p = &particle[track];
	#endif

	Real input_kapa;
	Real input_mu_plus;
	Real input_mu_minus;
	Real input_Dphi;
	Real input_L;

	string name = input_name;
	stringstream address(name);
	ifstream is;
	is.open(address.str().c_str());
	if (!is.is_open())
		return false;

	boost::replace_all(name, "-r-v.bin", "");
	boost::replace_all(name, "rho=", "");
	boost::replace_all(name, "-k=", "");
	boost::replace_all(name, "-mu+=", "\t");
	boost::replace_all(name, "-mu-=", "\t");
	boost::replace_all(name, "-Dphi=", "\t");
	boost::replace_all(name, "-L=", "\t");

	stringstream ss_name(name);
	ss_name >> box.density;
	ss_name >> input_kapa;
	ss_name >> input_mu_plus;
	ss_name >> input_mu_minus;
	ss_name >> input_Dphi;
	ss_name >> input_L;
	if (input_L != Lx_int)
	{
		cout << "The specified box size " << input_L << " is not the same as the size in binary file which is " << Lx_int << " please recompile the code with the right Lx_int in parameters.h file." << endl;
		return false;
	}

	Particle::kapa = input_kapa;
	Particle::mu_plus = input_mu_plus;
	Particle::mu_minus = input_mu_minus;
	Particle::D_phi = input_Dphi;

	is.read((char*) &box.N, sizeof(int) / sizeof(char));
	if (box.N < 0 || box.N > 1000000)
		return (false);

	for (int i = 0; i < box.N; i++)
	{
		is >> box.particle[i].r;
		is >> box.particle[i].v;
	}

	is.close();

	Init_Topology(box); // Adding walls

	cout << "number_of_particles = " << box.N << endl; // Printing number of particles.

	box.Update_Cells();

	box.info.str("");
	box.info << "rho=" << box.density <<  "-mu+=" << Particle::mu_plus << "-mu-=" << Particle::mu_minus << "-Dphi=" << Particle::D_phi << "-L=" << Lx;
	return (true);
}

int main(int argc, char *argv[])
{
	#ifdef COMPARE
		C2DVector::Init_Rand(seed);
	#else
		C2DVector::Init_Rand(time(NULL));
	#endif

	Box box;

	Init(&box, atof(argv[1]), atof(argv[2]), atof(argv[3]), atof(argv[4]));

	cout << box.info.str() << endl;

	stringstream address;
	address.str("");
	address << box.info.str() << "-r-v.bin";
	ofstream out_file(address.str().c_str());

	equilibrium(&box, equilibrium_step, saving_period, out_file);
	data_gathering(&box, total_step, saving_period, out_file);
	out_file.close();
}

