#ifndef __LyapunovBox__
#define __LyapunovBox__

#include "../shared/parameters.h"
#include "../shared/c2dvector.h"
#include "../shared/particle.h"
#include "../shared/cell.h"
#include "../shared/vector-set.h"
#include "../serial/box.h"
#include "node.h"
#include "mpi.h"

class LyapunovBox: public Box{
public:
	int thisnode, totalnode;
	int track_id;
	VectorSet us,vs,dvs; // The us (unit set) is the unit vector showing direction of the largest lyapunov exponents.
	vector<Real> t,tau;
	vector<GrowthRatio> ratio;

	ofstream outfile, trajfile;

	Node* node_head; // Node is a class that has information about the node_id and its boundaries, neighbores and etc.

	LyapunovBox();

	void Track_Particle(vector<Particle>&);

	void Send_State_Hyper_Vector(const State_Hyper_Vector& shv, int dest, int tag);
	void Recv_State_Hyper_Vector(const State_Hyper_Vector& shv, int source, int tag);
	void Root_Bcast_State_Hyper_Vector(const State_Hyper_Vector& shv);
	void Root_Gather_Vector_Set(const VectorSet& v);
	void Root_Bcast_Vector_Set(const VectorSet& v);

	void Init_Deviation(int direction_num);
	void Init_Time(const Real, const Real);

	void Evolution_Reorthonormalize(bool save);
	Real Recursive_Orthonormalize(Real interval, int iteration_num, bool save);
	Real Find_Growth(bool save);

	Real Lyapunov_Exponent(const Real, const Real, const Real, const Real, const int); // Finding the largest lyapunov exponent
	Real Lyapunov_Exponent(const Real, const Real, const Real, const Real, const Real, const Real, const int); // Finding the largest lyapunov exponent
	Real Recursive_Lyapunov_Exponent(const Real, const int, const Real, const int, const int);
	Real Simple_Lyapunov_Exponent(const Real eq_interval, const Real eq_duration, const Real interval, const Real duration, const int direction_num);

	Real Polarization();
};

LyapunovBox::LyapunovBox() : Box()
{
	MPI_Comm_rank(MPI_COMM_WORLD, &thisnode);
	MPI_Comm_size(MPI_COMM_WORLD, &totalnode);
}

void LyapunovBox::Track_Particle(vector<Particle>& trajectory)
{
	trajectory.push_back(particle[track_id]);
}

void LyapunovBox::Send_State_Hyper_Vector(const State_Hyper_Vector& shv, int dest, int tag)
{
// Allocation
	double* data_buffer = new double[3*N];

// Assigning the buffer arrays
	int counter = 0; // count particles. We can not use a shift very simply becasue we need to write both index_buffer and data_buffer.
	for (int i = 0; i < N; i++)
	{
		data_buffer[3*i] = shv.particle[i].r.x;
		data_buffer[3*i+1] = shv.particle[i].r.y;
		data_buffer[3*i+2] = shv.particle[i].theta;
		counter++;
	}

	MPI_Send(data_buffer, 3*N, MPI_DOUBLE, dest, 0,MPI_COMM_WORLD);
	MPI_Send(shv.gsl_r->state, shv.gsl_r->type->size, MPI_BYTE, dest, 0,MPI_COMM_WORLD);

	delete [] data_buffer;
}

void LyapunovBox::Recv_State_Hyper_Vector(const State_Hyper_Vector& shv, int source, int tag)
{
// Allocation
	double* data_buffer = new double[3*N];

	MPI_Status status;
	MPI_Recv(data_buffer, 3*N, MPI_DOUBLE, source, 0,MPI_COMM_WORLD, &status);

// Assigning the buffer arrays
	int counter = 0; // count particles. We can not use a shift very simply becasue we need to write both index_buffer and data_buffer.
	for (int i = 0; i < N; i++)
	{
		shv.particle[i].r.x = data_buffer[3*i];
		shv.particle[i].r.y = data_buffer[3*i+1];
		shv.particle[i].theta = data_buffer[3*i+2];
		counter++;
	}

	MPI_Recv(shv.gsl_r->state, shv.gsl_r->type->size, MPI_BYTE, source, 0,MPI_COMM_WORLD, &status);

	delete [] data_buffer;
}

void LyapunovBox::Root_Bcast_State_Hyper_Vector(const State_Hyper_Vector& shv)
{
// Allocation
	double* data_buffer = new double[3*N];

// Assigning the buffer arrays
	if (thisnode == 0)
		for (int i = 0; i < N; i++)
		{
			data_buffer[3*i] = shv.particle[i].r.x;
			data_buffer[3*i+1] = shv.particle[i].r.y;
			data_buffer[3*i+2] = shv.particle[i].theta;
		}

	MPI_Barrier(MPI_COMM_WORLD);

	MPI_Bcast(data_buffer, 3*N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(shv.gsl_r->state, shv.gsl_r->type->size, MPI_BYTE, 0,MPI_COMM_WORLD);

	if (thisnode != 0)
	{
		for (int i = 0; i < N; i++)
		{
			shv.particle[i].r.x = data_buffer[3*i];
			shv.particle[i].r.y = data_buffer[3*i+1];
			shv.particle[i].theta = data_buffer[3*i+2];
		}
	}

	delete [] data_buffer;

	MPI_Barrier(MPI_COMM_WORLD);
}


void LyapunovBox::Root_Gather_Vector_Set(const VectorSet& v)
{
	for (int i = 1 ; i < v.direction_num; i++)
	{
		if (thisnode != 0)
		{
			if (i % totalnode == thisnode)
				Send_State_Hyper_Vector(v.v[i],0,thisnode);
		}
		else
		{
			if (i % totalnode != 0)
				Recv_State_Hyper_Vector(v.v[i],i % totalnode,i % totalnode);
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);
}

void LyapunovBox::Root_Bcast_Vector_Set(const VectorSet& v)
{
	for (int i = 1; i < v.direction_num; i++)
	{
		if (thisnode != 0)
		{
			if (i % totalnode == thisnode)
				Recv_State_Hyper_Vector(v.v[i], 0, i % totalnode);
		}
		else
		{
			if (i % totalnode != 0)
				Send_State_Hyper_Vector(v.v[i], i % totalnode, i % totalnode);
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);
}


void LyapunovBox::Init_Deviation(int direction_num)
{
	us.direction_num = direction_num;
	us.particle_num = N;
	us.amplitude = 1e-10;
	us.Init();
	us.Rand();
// 	for (int i = 1; i < direction_num; i++)
// 	{
// 		us.v[i].Null();
// 		if ((i-1) / 2 < us.v[i].N)
// 		{
// 			if (i % 2 == 1)
// 				us.v[i].particle[(i-1) / 2].r.x = 1;
// 			else
// 				us.v[i].particle[(i-1) / 2].r.y = 1;
// 		}
// 		else
// 			us.v[i].particle[i - 1 - 2*us.v[i].N].theta = 1;
// 	}
// 	if (thisnode == 0)
// 		cout << us << endl;
	vs.Init();
	dvs.Init();
	GrowthRatio::direction_num = direction_num;
}

void LyapunovBox::Init_Time(const Real interval, const Real durution)
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

void LyapunovBox::Evolution_Reorthonormalize(bool save = false)
{
	if (thisnode == 0)
	{
		Save(vs.v[0]);
		dvs = us;
		dvs.Scale();
		if (save)
		{
			outfile << 0;
			for (int i = 0; i < us.direction_num; i++)
				outfile << "\t" << 1;
			outfile << endl;
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);
	Root_Bcast_Vector_Set(dvs);

	Real tt = 0;
	for (int j = 0; j < tau.size(); j++)
	{
		tt += tau[j];

		if (thisnode == 0 && (dt*tt > 100))
		{
			tt = 0;
			cout << "System is in time " << dt*t[j] << endl;
//			if (save)
//				cout << us << endl;
		}

		Root_Bcast_State_Hyper_Vector(vs.v[0]);
		Root_Bcast_Vector_Set(dvs);

		for (int i = 0; i < us.direction_num; i++)
		{
			if (i % totalnode == thisnode)
			{
				vs.v[i] = vs.v[0] + dvs.v[i];
				Load(vs.v[i]);
				Multi_Step(tau[j], 20);
				Save(vs.v[i]);
			}
		}

		Root_Gather_Vector_Set(vs);

		if (thisnode == 0)
		{
			for (int i = 1; i < us.direction_num; i++)
			{
				dvs.v[i] = (vs.v[i] - vs.v[0]);
			}
//			dvs.v[1].particle[0].theta = 1e-12;
//			dvs.v[2].particle[0].theta = 1e-9;
//			cout << "unperturbed\t" << setprecision(20) << sin(vs.v[0].particle[0].theta) << endl;
//			cout << "perturbed1\t" << setprecision(20) << sin((vs.v[0].particle[0].theta + dvs.v[1].particle[0].theta)) << endl;
//			cout << "perturbed2\t" << setprecision(20) << sin((vs.v[0] + dvs.v[2]).particle[0].theta) << endl;
			dvs.Renormalize(us);
			for (int i = 1; i < us.direction_num; i++)
			{
//				dvs.v[i] = us.v[i]*dvs.v[i].Magnitude();
				dvs.v[i] = us.v[i]*(dvs.v[i] * us.v[i]);
				Real temp_ratio = (dvs.v[i].Magnitude()) / (us.amplitude);
				ratio[j].r[i] = temp_ratio;
				ratio[j].r2[i] = temp_ratio*temp_ratio;
			}
			if (save)
			{
				outfile << dt*t[j];
				for (int i = 1; i < us.direction_num; i++)
					outfile << "\t" << ratio[j].r[i];
				outfile << endl;
			}
		}
	}
	if (thisnode == 0)
		Load(vs.v[0]);

	MPI_Barrier(MPI_COMM_WORLD);
}

Real LyapunovBox::Recursive_Orthonormalize(Real interval, int iteration_num, bool save = false)
{
	Real lambda;
	State_Hyper_Vector gamma0(N);

// Saving the initial state and Making deviations
	if (thisnode == 0)
	{
		Save(gamma0);
		vs.v[0] = gamma0;
		dvs = us;
		dvs.Scale();
		if (save)
		{
			outfile << 0;
			for (int i = 0; i < us.direction_num; i++)
				outfile << "\t" << 1;
			outfile << endl;
		}
	}

// Broad casting the deviations
	MPI_Barrier(MPI_COMM_WORLD);
	Root_Bcast_State_Hyper_Vector(gamma0);

	Real tt = 0;

	int tau = (int) round(interval / dt); // tau is number of steps in the interval (integer)
	for (int j = 0; j < iteration_num; j++)
	{
		tt += dt*tau;
		if (thisnode == 0 && (tt > 100))
		{
			tt = 0;
			cout << "System is in iteration " << j << endl;
//			if (save)
//				cout << us << endl;
		}

		vs.v[0] = gamma0;
// Broad casting the deviations
		MPI_Barrier(MPI_COMM_WORLD);
		Root_Bcast_Vector_Set(dvs);

		for (int i = 0; i < us.direction_num; i++)
		{
			if (i % totalnode == thisnode)
			{
				vs.v[i] = gamma0 + dvs.v[i];
				Load(vs.v[i]);
				Multi_Step(tau, 20);
				Save(vs.v[i]);
			}
		}

		MPI_Barrier(MPI_COMM_WORLD);
		Root_Gather_Vector_Set(vs);

		if (thisnode == 0)
		{
			for (int i = 1; i < us.direction_num; i++)
			{
				dvs.v[i] = (vs.v[i] - vs.v[0]);
			}

			dvs.Renormalize(us);
			Real temp_ratio[us.direction_num];
			for (int i = 1; i < us.direction_num; i++)
			{
//				dvs.v[i] = us.v[i]*dvs.v[i].Magnitude();
				dvs.v[i] = us.v[i]*(dvs.v[i] * us.v[i]);
				temp_ratio[i] = (dvs.v[i].Magnitude()) / (us.amplitude);
			}
			if (save)
			{
				outfile << dt*j*tau;
				for (int i = 1; i < us.direction_num; i++)
					outfile << "\t" << temp_ratio[i];
				outfile << endl;
			}
			lambda = log(temp_ratio[1]) / (dt*j*tau);
		}
	}

	if (thisnode == 0)
		Load(gamma0);

	MPI_Barrier(MPI_COMM_WORLD);

	return (lambda);
}

// Finding growth without orthonormalization
Real LyapunovBox::Find_Growth(bool save = false)
{
	Real lambda;
	State_Hyper_Vector gamma0(N);

// Saving the initial state and Making deviations
	if (thisnode == 0)
	{
		Save(gamma0);
		vs.v[0] = gamma0;
		dvs = us;
		dvs.Scale();
		if (save)
		{
			outfile << 0;
			for (int i = 0; i < us.direction_num; i++)
				outfile << "\t" << 1;
			outfile << endl;
		}
	}

// Broad casting the deviations
	MPI_Barrier(MPI_COMM_WORLD);
	Root_Bcast_State_Hyper_Vector(gamma0);
	Root_Bcast_Vector_Set(dvs);
	MPI_Barrier(MPI_COMM_WORLD);

	Real tt = 0;
	vs.v[0] = gamma0;

	for (int i = 0; i < us.direction_num; i++)
	{
		if (i % totalnode == thisnode)
			vs.v[i] = gamma0 + dvs.v[i];
	}

// Going over directions
	for (int j = 0; j < tau.size(); j++)
	{
		tt += tau[j]*dt;

		for (int i = 0; i < us.direction_num; i++)
		{
			if (i % totalnode == thisnode)
			{
				Load(vs.v[i]);
				Multi_Step(tau[j], 20);
				Save(vs.v[i]);
			}
		}

		MPI_Barrier(MPI_COMM_WORLD);
		Root_Gather_Vector_Set(vs);

		if (thisnode == 0)
		{
			for (int i = 1; i < us.direction_num; i++)
				dvs.v[i] = (vs.v[i] - vs.v[0]);

			Real temp_ratio[us.direction_num];
			for (int i = 1; i < us.direction_num; i++)
			{
				temp_ratio[i] = (dvs.v[i].Magnitude()) / (us.amplitude);
			}
			if (save)
			{
				outfile << tt;
				for (int i = 1; i < us.direction_num; i++)
					outfile << "\t" << temp_ratio[i];
				outfile << endl;
			}
			lambda = log(temp_ratio[1]) / (tau[j]);

			gamma0 = vs.v[0];
		}
		MPI_Barrier(MPI_COMM_WORLD);
		Root_Bcast_State_Hyper_Vector(gamma0);
	}
	MPI_Barrier(MPI_COMM_WORLD);

	return (lambda);
}

// Finding the largest lyapunov exponent
Real LyapunovBox::Lyapunov_Exponent(const Real eq_interval, const Real eq_duration, const Real interval, const Real duration, const int direction_num)
{
	clock_t start_time, end_time;
	start_time = clock();
	
	Init_Deviation(direction_num);
	Init_Time(eq_interval, eq_duration);
	Evolution_Reorthonormalize(false);

	end_time = clock();
	if (thisnode == 0)
		cout << "Finded unit orthonormal vectors in: " << (end_time - start_time) / CLOCKS_PER_SEC << " s" << endl;
	start_time = end_time;
	
	Init_Time(interval, duration);
	if (thisnode == 0)
		cout << "Initialized time set for lyapunov computation " << endl;
	Evolution_Reorthonormalize(true);

	end_time = clock();
	Real running_time = (end_time - start_time) / CLOCKS_PER_SEC;
	if (thisnode == 0)
		cout << "Finded unit orthonormal vectors in: " << (end_time - start_time) / CLOCKS_PER_SEC << " s" << endl;
	return (running_time);
}

Real LyapunovBox::Lyapunov_Exponent(const Real eq_interval0, const Real eq_duration0, const Real eq_interval1, const Real eq_duration1, const Real interval, const Real duration, const int direction_num)
{
	clock_t start_time, end_time;
	start_time = clock();

	Init_Deviation(direction_num);
	Init_Time(eq_interval0, eq_duration0);
	Evolution_Reorthonormalize(false);

	end_time = clock();
	if (thisnode == 0)
	{
		cout << "Finded unit orthonormal vectors in: " << (end_time - start_time) / CLOCKS_PER_SEC << " s" << endl;
	}
	start_time = end_time;
	
	Init_Time(eq_interval1, eq_duration1);
	Evolution_Reorthonormalize(false);

	end_time = clock();
	if (thisnode == 0)
	{
		cout << "Finded unit orthonormal vectors more precise in: " << (end_time - start_time) / CLOCKS_PER_SEC << " s" << endl;
	}
	start_time = end_time;	

	Init_Time(interval, duration);
	if (thisnode == 0)
		cout << "Initialized time set for lyapunov computation " << endl;
	Evolution_Reorthonormalize(true);

	end_time = clock();
	Real running_time = (end_time - start_time) / CLOCKS_PER_SEC;
	if (thisnode == 0)
		cout << "Finded unit orthonormal vectors in: " << (end_time - start_time) / CLOCKS_PER_SEC << " s" << endl;

	return (running_time);
}

Real LyapunovBox::Recursive_Lyapunov_Exponent(const Real interval0, const int iteration0, const Real interval, const int iteration, const int direction_num)
{
	clock_t start_time, end_time;
	start_time = clock();

	Init_Deviation(direction_num);
	Recursive_Orthonormalize(interval0, iteration0, false);

	end_time = clock();
	if (thisnode == 0)
	{
		cout << "Finded unit orthonormal vectors in: " << (end_time - start_time) / CLOCKS_PER_SEC << " s" << endl;
	}
	start_time = end_time;

	Real lambda = Recursive_Orthonormalize(interval, iteration, false);

	end_time = clock();
	Real running_time = (end_time - start_time) / CLOCKS_PER_SEC;
	if (thisnode == 0)
		cout << "Finded unit orthonormal vectors in: " << (end_time - start_time) / CLOCKS_PER_SEC << " s" << endl;

	return (lambda);
}


// Lyapunov exponent without rotation of deviations. This procedure must have a diffusion limit for large time
Real LyapunovBox::Simple_Lyapunov_Exponent(const Real eq_interval, const Real eq_duration, const Real interval, const Real duration, const int direction_num)
{
	clock_t start_time, end_time;
	start_time = clock();
	
	Init_Deviation(direction_num);
	Init_Time(eq_interval, eq_duration);
	Evolution_Reorthonormalize(false);

	end_time = clock();
	if (thisnode == 0)
		cout << "Finded unit orthonormal vectors in: " << (end_time - start_time) / CLOCKS_PER_SEC << " s" << endl;
	start_time = end_time;
	
	Init_Time(interval, duration);
	if (thisnode == 0)
		cout << "Initialized time set for lyapunov computation " << endl;
	Find_Growth(true);

	end_time = clock();
	Real running_time = (end_time - start_time) / CLOCKS_PER_SEC;
	if (thisnode == 0)
		cout << "Finded unit orthonormal vectors in: " << (end_time - start_time) / CLOCKS_PER_SEC << " s" << endl;
	return (running_time);
}

Real LyapunovBox::Polarization()
{
	Real phi;
	C2DVector p;
	p.x = p.y = 0;
	if (thisnode == 0)
	{
		for (int i = 0; i < N; i++)
			p += particle[i].v;
		p /= N;
		phi = sqrt(p.Square());
	}
	return (phi);
}

#endif
