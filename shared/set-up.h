#ifndef _SETUP_
#define _SETUP_

// This will set up the box for simulation. That means distributing particles in the box and their formation like random velocity or a clump or a vortex. This will also defines geometry that the particles are interacting with. 

#include "../shared/parameters.h"
#include "../shared/c2dvector.h"
#include "../shared/particle.h"
#include "../shared/geometry.h"


void Square_Lattice_Formation(Particle* particle, int N); // Positioning partilces in a square lattice
void Triangle_Lattice_Formation(Particle* particle, int N, double sigma); // Positioning partilces in a triangular lattice. This is denser.
void Random_Formation(Particle* particle, int N); // Positioning partilces Randomly
void Random_Formation(Particle* particle, int N, double sigma); // Positioning partilces Randomly, but distant from walls
void Random_Formation_Circle(Particle* particle, int N, double r); // Positioning partilces Randomly in a circle with radius r
void Polar_Formation(Particle* particle, int N); // Polar state
void Single_Vortex_Formation(Particle* particle, int N); // Vortex initial condition.
void Four_Vortex_Formation(Particle* particle, int N); // Four vortex inside the box. left top, left bot, right top and right bot.
void Clump_Formation(Particle* particle, int N, int size); // Positioning particles in a clump that is moving in some direction.
void Square_Ring_Formation(Particle* particle, int N); // Positioning particles in a square shape ring(the center is empty)
void Star_Trap_Initialization(Geometry* geometry, int N_hands, Real half_delta); // Defining geometry of a star shaped trap


void Single_Vortex_Formation(Particle* particle, int N)
{
	C2DVector v_cm,v;
	v_cm.Null();

	Triangle_Lattice_Formation(particle, N, 1);

	for (int i = 0; i < N; i++)
		{
			v.x = -(particle[i].r.y);
			v.y = (particle[i].r.x);
			v = v / sqrt(v.Square());
			particle[i].Init(particle[i].r,v);
			v_cm += (particle[i].v / N);
		}
}

void Four_Vortex_Formation(Particle* particle, int N)
{
	C2DVector v_cm,r,v;
	int Nx = (int) sqrt(Lx*N/Ly) + 1;
	int Ny = (int) sqrt(Ly*N/Lx) + 1;
	for (int i = 0; i < Nx/2; i++)
		for (int j = 0; j < Ny/2; j++)
		{
			r.x = -Lx + 2*i*(Lx/Nx) + Lx/Nx;
			r.y = -Ly + 2*j*(Ly/Ny) + Ly/Ny;
			v.x = -(2*(r.y/Lx) + 1);
			v.y = (2*(r.x/Ly) + 1);
			particle[i + j*Nx].Init(r,v);
			v_cm += (particle[i + j*Nx].v / N);
		}
	for (int i = Nx/2; i < Nx; i++)
		for (int j = 0; j < Nx/2; j++)
		{
			r.x = -Lx + 2*i*(Lx/Nx) + Lx/Nx;
			r.y = -Ly + 2*j*(Ly/Ny) + Ly/Ny;
			v.x = (2*(r.y/Lx) + 1);
			v.y = -(2*(r.x/Ly) - 1);
			particle[i + j*Nx].Init(r,v);
			v_cm += (particle[i + j*Nx].v / N);
		}
	for (int i = 0; i < Nx/2; i++)
		for (int j = Nx/2; j < Nx; j++)
		{
			r.x = -Lx + 2*i*(Lx/Nx) + Lx/Nx;
			r.y = -Ly + 2*j*(Ly/Ny) + Ly/Ny;
			particle[i + j*Nx].Init(r);
			v.x = (2*(r.y/Ly) - 1);
			v.y = -(2*(r.x/Lx) + 1);
			particle[i + j*Nx].Init(r,v);
			v_cm += (particle[i + j*Nx].v / N);
		}
	for (int i = Nx/2; i < Nx; i++)
		for (int j = Nx/2; j < Nx; j++)
		{
			r.x = -Lx + 2*i*(Lx/Nx) + Lx/Nx;
			r.y = -Ly + 2*j*(Ly/Ny) + Ly/Ny;
			particle[i + j*Nx].Init(r);
			v.x = -(2*(r.y/Ly) - 1);
			v.y = (2*(r.x/Lx) - 1);
			particle[i + j*Nx].Init(r,v);
			v_cm += (particle[i + j*Nx].v / N);
		}
}

void Square_Lattice_Formation(Particle* particle, int N)
{
	C2DVector v_cm,r;
	int Nx = (int) sqrt(Lx*N/Ly) + 1;
	int Ny = (int) sqrt(Ly*N/Lx) + 1;
	for (int i = 0; i < N; i++)
	{
		r.x = ((2.0*(i / Ny))/(Nx) - 1)*Lx;
		r.y = ((2.0*(i % Ny))/(Ny) - 1)*Ly;
		r = r*0.95;

		particle[i].Init(r);
		v_cm += (particle[i].v / N);
	}
}

void Triangle_Lattice_Formation(Particle* particle, int N, double sigma)
{
	C2DVector v_cm,r,basis_1, basis_2;
	basis_1.x = 1;
	basis_1.y = 0;
	basis_2.x = 0.5;
	basis_2.y = sqrt(3)/2;

	int Nx = (int) sqrt((Lx*N*sqrt(3))/(2*Ly)) + 1;
	int Ny = N / Nx + 1;

	basis_1 = basis_1*((Lx2 - 2*sigma) / Nx);
	basis_2 = basis_2*((Lx2 - 2*sigma) / Nx);

//	basis_1 = basis_1*(sigma*1.15);
//	basis_2 = basis_2*(sigma*1.15);
//	int Nx = (int) (L2 / basis_1.x);

	for (int i = 0; i < N; i++)
	{
		r = basis_1*(i % Nx) + basis_2*(i / Nx) - basis_1*(i / (2*Nx));
		r.x -= (Lx - sigma);
		r.y -= (Ly - sigma);
		if ((r.y > Ly) || (r.y < -Ly) || (r.x > Lx) || (r.x < -Lx))
			cout << r.x/Lx << "\t" << r.y/Ly << endl;
//		r.Periodic_Transform();
		particle[i].Init(r);
		v_cm += (particle[i].v / N);
	}
}

void Random_Formation(Particle* particle, int N)
{
	Random_Formation(particle, N, 0);
}

void Random_Formation(Particle* particle, int N, double sigma)
{
	C2DVector r;
	for (int i = 0; i < N; i++)
	{
		r.Rand(Lx-sigma, Ly-sigma);
		if ((r.y > Ly) || (r.y < -Ly) || (r.x > Lx) || (r.x < -Lx))
			cout << r.x/Lx << "\t" << r.y/Ly << endl;
		r.Periodic_Transform();
		particle[i].Init(r);
	}
}

void Random_Formation_Circle(Particle* particle, int N, double radius)
{
	C2DVector r;
	for (int i = 0; i < N; i++)
	{
		Real a = gsl_ran_flat(r.gsl_r, 0, (radius*radius)/2);
		a = sqrt(2*a);
		Real b = gsl_ran_flat(r.gsl_r, -M_PI, M_PI);
		r.x = a*cos(b);
		r.y = a*sin(b);
		r.Periodic_Transform();
		particle[i].Init(r);
	}
}

// Polar state
void Polar_Formation(Particle* particle, int N)
{
//	Triangle_Lattice_Formation(particle, N, 0);
	Random_Formation(particle,N);

	for (int i = 0; i < N; i++)
	{
		particle[i].theta = 0;
		particle[i].v.x = 1;
		particle[i].v.y = 0;
	}
}

void Clump_Formation(Particle* particle, int N, int size)
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


void Square_Ring_Formation(Particle* particle, int N)
{
	C2DVector r;

	for (int i = 0; i < N; i++)
	{
		// particles appear in a box [-L/2.,L/2.]
		particle[i].Init(Lx/2.);
		// moving particles to a the outer shell 
		if (fabs(particle[i].r.x) > Lx/4. || fabs(particle[i].r.y) > Ly/4.)
			particle[i].r *= 2.;
		else if (particle[i].r.x > 0. && particle[i].r.y > 0.)
		{
			particle[i].r.x += Lx/2.;
			particle[i].r.y += Ly/2.;
		}
		else if (particle[i].r.x > 0. && particle[i].r.y < 0.)
		{
			particle[i].r.x += Lx/2.;
			particle[i].r.y -= Ly/2.;
		}
		else if (particle[i].r.x < 0. && particle[i].r.y > 0.)
		{
			particle[i].r.x -= Lx/2.;
			particle[i].r.y += Ly/2.;
		}
		else if (particle[i].r.x < 0. && particle[i].r.y < 0.)
		{
			particle[i].r.x -= Lx/2.;
			particle[i].r.y -= Ly/2.;
		}
	}
}

int r_big = 15;
int r_small = 6;

void Star_Trap_Initialization(Geometry* geometry, int N_hands, Real half_delta)
{
	// Star Trap (points)
	C2DVector end_point[100];
	Real theta, alpha, beta;
	theta = 0.; 
	end_point[0].x = r_big*cos(theta);
	end_point[0].y = r_big*sin(theta);

	alpha = theta + (1./2.)*(2.*PI/N_hands) - half_delta;
	beta  = theta + (1./2.)*(2.*PI/N_hands) + half_delta;

	end_point[1].x = r_small*cos(alpha);
	end_point[1].y = r_small*sin(alpha);

	end_point[2].x = r_small*cos(beta);
	end_point[2].y = r_small*sin(beta);

	for (int i = 3; i < 3*N_hands; i++)
		end_point[i] = end_point[i-3].Rotate(2.*PI/N_hands);

	// Star Trap (walls)
	for (int i = 0; i < (3*(N_hands-1)+1); i+=3)
	{
		geometry->Add_Wall(end_point[i], end_point[i+1]);
		geometry->Add_Wall(end_point[i+2], end_point[(i+3)%(3*N_hands)]);
	}
}


#endif
