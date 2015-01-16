#include "../shared/parameters.h"
#include "../shared/c2dvector.h"
#include "../shared/particle.h"
#include "../shared/cell.h"
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

void Init(Box* box, Real input_density, Real g, Real alpha, Real noise_amplitude)
{

	box->wall_num = 4;
	box->density = input_density;
	box->N = (int) round(Lx2*Ly2*box->density);
	cout << "number_of_particles = " << box->N << endl;

	Particle::noise_amplitude = noise_amplitude / sqrt(dt);
	ContinuousParticle::g = g;
	ContinuousParticle::alpha = alpha;

	Triangle_Lattice_Formation(box->particle, box->N);
//	Single_Vortex_Formation(particle, N);
//	Four_Vortex_Formation(particle, N);

	#ifndef PERIODIC_BOUNDARY_CONDITION
		box->wall[0].Init(-Lx,-Ly,-Lx, Ly);
		box->wall[1].Init(-Lx, Ly, Lx, Ly);
		box->wall[2].Init( Lx, Ly, Lx,-Ly);
		box->wall[3].Init( Lx,-Ly,-Lx,-Ly);
	#endif

	box->Update_Cells();

	box->info.str("");
	box->info << "rho=" << box->density <<  "-g=" << Particle::g << "-alpha=" << Particle::alpha << "-noise=" << noise_amplitude;
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

