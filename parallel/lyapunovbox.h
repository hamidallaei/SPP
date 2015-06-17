#ifndef __LyapunovBox__
#define __LyapunovBox__

#include "../shared/parameters.h"
#include "../shared/c2dvector.h"
#include "../shared/particle.h"
#include "../shared/cell.h"
#include "../shared/vector-set.h"
#include "box.h"

class LyapunovBox: public Box{
public:
	VectorSet us,vs,vs0; // The us (unit set) is the unit vector showing direction of the largest lyapunov exponents.
	vector<Real> t,tau;
	vector<GrowthRatio> ratio;
	vector<State_Hyper_Vector> gamma;

	ofstream outfile;
	
	void Init_Deviation(int direction_num);
	void Init_Time(const int, const int);
	
	void Add_Deviation(const State_Hyper_Vector&); // Add argument state vector as a deviation to the state of the box
	void Evolution();
	Real Lyapunov_Exponent(const int, const int, const int, const int, const int); // Finding the largest lyapunov exponent
};



void LyapunovBox::Init_Deviation(int direction_num)
{
	us.direction_num = direction_num;
	us.particle_num = N;
	us.amplitude = 1e-6;
	us.Init();
	us.Rand();
	vs.Init();
	vs0.Init();
	GrowthRatio::direction_num = direction_num;
}


// Add argument state vector as a deviation to the state of the box
void LyapunovBox::Add_Deviation(const State_Hyper_Vector& dsv)
{
	thisnode->Root_Gather();
	if (N != dsv.N)
	{
		cout << "Error: Number of particles in state vectors differ from box" << endl;
		exit(0);
	}
	for (int i = 0; i < N; i++)
	{
		particle[i].r += dsv.particle[i].r;
		particle[i].r.Periodic_Transform();
		particle[i].theta += dsv.particle[i].theta;
		particle[i].v.x = cos(particle[i].theta);
		particle[i].v.y = sin(particle[i].theta);
	}
	thisnode->Root_Bcast();
	thisnode->Full_Update_Cells();
	#ifdef verlet_list
		thisnode->Update_Neighbor_List();
	#endif
}


void LyapunovBox::Init_Time(const int interval, const int durution)
{
	t.clear();
	tau.clear();
	ratio.clear();
	for (Real i = interval; i < durution; i+=interval)
	{
		t.push_back((int) round(i/dt));
		GrowthRatio temp_ratio;
		ratio.push_back(temp_ratio);
	}
	tau.push_back(t[0]);
	for (int i = 0; i < t.size() - 1; i++)
	{
		tau.push_back(t[i+1]-t[i]);
	}
}


void LyapunovBox::Evolution()
{
	vs0 = us;
	vs0.Scale();

	static State_Hyper_Vector gamma_0(N);
	static State_Hyper_Vector gamma_prime(N);
	static State_Hyper_Vector temp_gamma(N);

	gamma.clear();
	Save(gamma_0);

	for (int i = 0; i < tau.size(); i++)
	{
		Multi_Step(tau[i], 20);
		Save(temp_gamma);
		gamma.push_back(temp_gamma);
	}

	for (int i = 0; i < us.direction_num; i++)
	{
		Load(gamma_0);
		Add_Deviation(vs0.v[i]);
		for (int j = 0; j < tau.size(); j++)
		{
			Multi_Step(tau[j], 20);
			Save(gamma_prime);
			vs.v[i] = (gamma_prime - gamma[j]);
			Real temp = (vs.v[i] * us.v[i]) / (us.amplitude);
			ratio[j].r[i] = temp;
			ratio[j].r2[i] = temp*temp;
		}
	}

	vs.Renormalize(us);
	Load(gamma.back());
}

// Finding the largest lyapunov exponent
Real LyapunovBox::Lyapunov_Exponent(const int eq_interval, const int eq_duration, const int interval, const int duration, const int direction_num)
{
	clock_t start_time, end_time;
	start_time = clock();
	
	Init_Deviation(direction_num);
	Init_Time(eq_interval, eq_duration);
	Evolution();

	end_time = clock();
	if (thisnode->node_id == 0)
		cout << "Finded unit orthonormal vectors in: " << (end_time - start_time) / CLOCKS_PER_SEC << " s" << endl;
	start_time = end_time;
	
	Init_Time(interval, duration);
	Evolution();
	outfile << 0;
	for (int i = 0; i < direction_num; i++)
		outfile << "\t" << 1;
	outfile << endl;
	for (int i = 0; i < t.size(); i++)
	{
		outfile << dt*t[i];
		for (int j = 0; j < direction_num; j++)
			outfile << "\t" << ratio[i].r[j];
		outfile << endl;
	}

	end_time = clock();
	Real running_time = (end_time - start_time) / CLOCKS_PER_SEC;
	if (thisnode->node_id == 0)
		cout << "Finded unit orthonormal vectors in: " << (end_time - start_time) / CLOCKS_PER_SEC << " s" << endl;
	return (running_time);
}

#endif