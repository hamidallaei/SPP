#ifndef _BOX_
#define _BOX_

#include <iostream>
#include <sstream>
#include <fstream>
#include "c2dvector.h"
#include "parameters.h"
#include "particle.h"
#include "cell.h"
#include "wall.h"


class Box{
public:
	int N, wall_num;
	Particle particle[max_N];
	Wall wall[8];
	Cell cell[divisor_x][divisor_y];

	Real density, volume_fraction;
	stringstream info;

	Box();
	void Init(Real input_volume_fraction, Real g = 2, Real alpha=0.5, Real noise_amplitude = 0.1);
	void Square_Lattice_Formation();
	void Triangle_Lattice_Formation();
	void Single_Vortex_Formation();
	void Four_Vortex_Formation();
	void Clump_Formation(int size);
	void Update_Cells();
	void Interact();
	void Move();
	void One_Step();
	void Multi_Step(int steps);
	void Translate(C2DVector d);
	void Center();

	void Make_Traj(Real scale, std::ofstream& data_file);

	friend std::ostream& operator<<(std::ostream& os, Box* box);
	friend std::istream& operator>>(std::istream& is, Box* box);
};

Box::Box()
{
	N = 0;
	density = 0;
}

void Box::Init(Real input_density, Real g, Real alpha, Real noise_amplitude)
{
	Cell::particle = particle;

	wall_num = 4;
	density = input_density;
	N = (int) round(L2*L2*density);
	cout << "number_of_particles = " << N << endl;

	Particle::noise_amplitude = noise_amplitude / sqrt(dt);
	ContinuousParticle::g = g;
	ContinuousParticle::alpha = alpha;

	Triangle_Lattice_Formation();
//	Single_Vortex_Formation();
//	Four_Vortex_Formation();

	for (int i = 0; i < divisor_x; i++)
		for (int j = 0; j < divisor_y; j++)
			cell[i][j].Init((Real) L*(2*i-divisor_x + 0.5)/divisor_x, (Real) L*(2*j-divisor_y + 0.5)/divisor_y);

	Update_Cells();

	#ifndef PERIODIC_BOUNDARY_CONDITION
		wall[0].Init(-L,-L,-L,L);
		wall[1].Init(-L,L,L,L);
		wall[2].Init(L,L,L,-L);
		wall[3].Init(L,-L,-L,-L);
	#endif

	info.str("");
	info << "rho=" << density <<  "-g=" << Particle::g << "-alpha=" << Particle::alpha << "-noise=" << noise_amplitude;
}

void Box::Single_Vortex_Formation()
{
	C2DVector v_cm,v;
	v_cm.Null();

	Triangle_Lattice_Formation();

	for (int i = 0; i < N; i++)
		{
			v.x = -(particle[i].r.y);
			v.y = (particle[i].r.x);
			v = v / sqrt(v.Square());
			particle[i].Init(particle[i].r,v);
			v_cm += (particle[i].v / N);
		}

	for (int i = 0; i < N; i++)
		particle[i].v -= v_cm;
}

void Box::Four_Vortex_Formation()
{
	C2DVector v_cm,r,v;
	int Nx = (int) sqrt(N) + 1;
	for (int i = 0; i < Nx/2; i++)
		for (int j = 0; j < Nx/2; j++)
		{
			r.x = -L + 2*i*(L/Nx) + L/Nx;
			r.y = -L + 2*j*(L/Nx) + L/Nx;
			v.x = -(2*(r.y/L) + 1);
			v.y = (2*(r.x/L) + 1);
			particle[i + j*Nx].Init(r,v);
			v_cm += (particle[i + j*Nx].v / N);
		}
	for (int i = Nx/2; i < Nx; i++)
		for (int j = 0; j < Nx/2; j++)
		{
			r.x = -L + 2*i*(L/Nx) + L/Nx;
			r.y = -L + 2*j*(L/Nx) + L/Nx;
			v.x = (2*(r.y/L) + 1);
			v.y = -(2*(r.x/L) - 1);
			particle[i + j*Nx].Init(r,v);
			v_cm += (particle[i + j*Nx].v / N);
		}
	for (int i = 0; i < Nx/2; i++)
		for (int j = Nx/2; j < Nx; j++)
		{
			r.x = -L + 2*i*(L/Nx) + L/Nx;
			r.y = -L + 2*j*(L/Nx) + L/Nx;
			particle[i + j*Nx].Init(r);
			v.x = (2*(r.y/L) - 1);
			v.y = -(2*(r.x/L) + 1);
			particle[i + j*Nx].Init(r,v);
			v_cm += (particle[i + j*Nx].v / N);
		}
	for (int i = Nx/2; i < Nx; i++)
		for (int j = Nx/2; j < Nx; j++)
		{
			r.x = -L + 2*i*(L/Nx) + L/Nx;
			r.y = -L + 2*j*(L/Nx) + L/Nx;
			particle[i + j*Nx].Init(r);
			v.x = -(2*(r.y/L) - 1);
			v.y = (2*(r.x/L) - 1);
			particle[i + j*Nx].Init(r,v);
			v_cm += (particle[i + j*Nx].v / N);
		}
	for (int i = 0; i < N; i++)
		particle[i].v -= v_cm;
}

void Box::Square_Lattice_Formation()
{
	C2DVector v_cm,r;
	int Nx = (int) sqrt(N) + 1;
	for (int i = 0; i < N; i++)
	{
		r.x = ((2.0*(i % Nx))/(Nx) - 1)*L;
		r.y = ((2.0*(i / Nx))/(Nx) - 1)*L;
		r = r*0.95;

		particle[i].Init(r);
		v_cm += (particle[i].v / N);
	}
	for (int i = 0; i < N; i++)
		particle[i].v -= v_cm;
}

void Box::Triangle_Lattice_Formation()
{
	C2DVector v_cm,r,basis_1, basis_2;
	basis_1.x = 1;
	basis_1.y = 0;
	basis_2.x = 0.5;
	basis_2.y = sqrt(3)/2;

	int Nx = (int) sqrt(N*sqrt(3)/2) + 1;
	int Ny = N / Nx + 1;

	basis_1 = basis_1*((L2 - 2*sigma) / Nx);
	basis_2 = basis_2*((L2 - 2*sigma) / Nx);

//	basis_1 = basis_1*(sigma*1.15);
//	basis_2 = basis_2*(sigma*1.15);
//	int Nx = (int) (L2 / basis_1.x);




	for (int i = 0; i < N; i++)
	{
		r = basis_1*(i % Nx) + basis_2*(i / Nx) - basis_1*(i / (2*Nx));
		r.x -= (L - sigma);
		r.y -= (L - sigma);
		if ((r.y > L) || (r.y < -L) || (r.x > L) || (r.x < -L))
			cout << r/L << endl;
//		r.Periodic_Transform();
		particle[i].Init(r);
		v_cm += (particle[i].v / N);
	}
	for (int i = 0; i < N; i++)
		particle[i].v -= v_cm;
}

void Box::Clump_Formation(int size)
{
	C2DVector v_cm, r, basis_1, basis_2;
	N = size;

	basis_1.x = 1;
	basis_1.y = 0;
	basis_2.x = 0.5;
	basis_2.y = sqrt(3)/2;

	int Nx = (int) sqrt(N*sqrt(3)/2) + 1;
	int Ny = N / Nx + 1;

	basis_1 *= 0.9;
	basis_2 *= 0.9;

	Real theta = gsl_ran_flat (C2DVector::gsl_r, -PI, PI);
	v_cm.x = cos(theta);
	v_cm.y = sin(theta);

	for (int i = 0; i < N; i++)
	{
		r = basis_1*(i % Nx) + basis_2*(i / Nx) - basis_1*(i / (2*Nx));
		particle[i].Init(r);
		particle[i].v = v_cm;
		particle[i].theta = theta;
	}
}

void Box::Update_Cells()
{
	for (int x = 0; x < divisor_x; x++)
		for (int y = 0; y < divisor_y; y++)
			cell[x][y].Delete();

	#pragma omp parallel for
	for (int i = 0; i < N; i++)
	{
		int x,y;
		x = (int) (particle[i].r.x + L)*divisor_x / L2;
		y = (int) (particle[i].r.y + L)*divisor_y / L2;

		#ifdef DEBUG
		if ((x >= divisor_x) || (x < 0) || (y >= divisor_y) || (y < 0))
		{
			cout << "\n Particle number " << i << " is Out of the box" << endl << flush;
			cout << "Particle Position is " << particle[i].r << endl;
			exit(0);
		}
		#endif

		cell[x][y].Add(i);
	}
}


void Box::Interact()
{
	#ifdef PERIODIC_BOUNDARY_CONDITION
	for (int x = 0; x < divisor_x; x++)
		for (int y = 0; y < divisor_y; y++)
		{
			cell[x][y].Self_Interact();
			cell[x][y].Interact(&cell[(x+1)%diviso_x][y]);
			cell[x][y].Interact(&cell[x][(y+1)%divisor_x]);
			cell[x][y].Interact(&cell[(x+1)%divisor_x][(y+1)%divisor_y]);
			cell[x][y].Interact(&cell[(x+1)%divisor_x][(y-1+divisor_y)%divisor_y]);
		}
	#else
	for (int x = 0; x < divisor_x-1; x++)
		for (int y = 1; y < divisor_y-1; y++)
		{
			cell[x][y].Self_Interact();
			cell[x][y].Interact(&cell[(x+1)%divisor_y][y]);
			cell[x][y].Interact(&cell[x][(y+1)%divisor_y]);
			cell[x][y].Interact(&cell[(x+1)%divisor_y][(y+1)%divisor_y]);
			cell[x][y].Interact(&cell[(x+1)%divisor_y][(y-1+divisor_y)%divisor_y]);
		}
	for (int y = 1; y < divisor_y-1; y++)
	{
		cell[divisor_x-1][y].Self_Interact();
		cell[divisor_x-1][y].Interact(&cell[divisor_x-1][y+1]);
	}
	for (int x = 0; x < divisor_x-1; x++)
	{
		cell[x][0].Self_Interact();
		cell[x][0].Interact(&cell[x+1][0]);
		cell[x][0].Interact(&cell[x][1]);
		cell[x][0].Interact(&cell[x+1][1]);

		cell[x][divisor_x-1].Self_Interact();
		cell[x][divisor_x-1].Interact(&cell[x+1][divisor_y-1]);
		cell[x][divisor_x-1].Interact(&cell[x+1][divisor_y-2]);
	}

	cell[divisor_x-1][divisor_y-1].Self_Interact();
	cell[divisor_x-1][0].Self_Interact();
	cell[divisor_x-1][0].Interact(&cell[divisor_y-1][1]);

	for (int i = 0; i < wall_num; i++)
		for (int j = 0; j < N; j++)
			wall[i].Interact(&particle[j]);
	#endif

//	for (int i = 0; i < N; i++)
//		for (int j = 0; j < 4; j++)
//			particle[i].Interact(&wall[j]);
}


void Box::Move()
{
	for (int i = 0; i < N; i++)
		particle[i].Move();
}

void Box::One_Step()
{
	Interact();
	Move();
}

void Box::Multi_Step(int steps)
{
	for (int i = 0; i < steps; i++)
	{
		Interact();
		Move();
	}
	Update_Cells();
}

void Box::Translate(C2DVector d)
{
	for (int i = 0; i < N; i++)
	{
		particle[i].r += d;
		particle[i].r.Periodic_Transform();
	}
	Update_Cells();
}

void Box::Center()
{
	C2DVector cm;
	for (int i = 0; i < N; i++)
		cm -= particle[i].r;
	cm = cm / N;
	Translate(cm);
}



void Box::Make_Traj(Real scale, ofstream& data_file)
{
	data_file << N << endl;
	data_file << "something" << endl;
	for (int i = 0; i < N; i++)
		data_file << "H	" << particle[i].r * scale << "\t" << 0.0 << endl;
}

std::ostream& operator<<(std::ostream& os, Box* box)
{
	os.write((char*) &box->N, sizeof(box->N) / sizeof(char));
	for (int i = 0; i < box->N; i++)
	{
		box->particle[i].r.write(os);
		C2DVector v;
		v.x = cos(box->particle[i].theta);
		v.y = sin(box->particle[i].theta);
		v.write(os);
	}
}

std::istream& operator>>(std::istream& is, Box* box)
{
	is.read((char*) &box->N, sizeof(int) / sizeof(char));
	for (int i = 0; i < box->N; i++)
	{
		is >> box->particle[i].r;
		is >> box->particle[i].v;
	}
}

#endif
