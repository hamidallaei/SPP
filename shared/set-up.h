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

	basis_1 = basis_1*((Lx2 - sigma) / Nx);
	basis_2 = basis_2*((Lx2 - sigma) / Nx);

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
		r.Periodic_Transform();
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


void Star_Trap_Initialization(Geometry* geometry, int N_hands, Real half_delta)
{
	int r_big = 15;
	int r_small = 6;

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


void Ring_Membrane(Particle* particle, const Real bead_diameter, int N_m)
{/*
	Put all of the membrane beads on a circular ring
	N_m = Number of membrane beads
*/
	C2DVector r;
/*	Real membrane_radius = 0.5*N_m*bead_diameter/M_PI;*/
	Real membrane_radius = 0.5*bead_diameter/sin(M_PI/N_m);
	Real theta =  0;

	cout << "Membrane Radius:\t" << membrane_radius << "\tBox dim:\t" << 2*Lx << endl;

	for (int i = 0; i < N_m; i++)
	{
		r.x = membrane_radius*cos(theta);
		r.y = membrane_radius*sin(theta);
		particle[i].r = r;

		theta += 2*M_PI /N_m;
	}
}


void Confined_In_Ring_Membrane(Particle* particle, const Real bead_diameter, int N_s, int N_m)
{/*
	Put all of the active beads inside the circular ring
	we divide active beads to some parts and then put each part on a ring 
		adjacent to each other and the membrane.
	N_s = Number of swimmers
	N_m = Number of membrane beads
	bead_diameter = Particle::sigma_p
*/
	C2DVector r;
/*	Real membrane_radius = 0.5*N_m*Particle::sigma_p / M_PI;*/
	Real membrane_radius = 0.5*bead_diameter/sin(M_PI/N_m);
	int chain_length = particle[N_m].nb; // number of beads in each swimmer(chain)

	int i = 0;
	int n_1 = 0;
	int ring_index = 1;
	Real ring_theta = 0;
	while (n_1 < N_s)
	{ /* ring_radius is radius of each ring of active particles
		n_ring is number of active particles on that ring
		ring_index renotes the index of each ring
		*/
		Real ring_radius = membrane_radius - ring_index * chain_length * (bead_diameter - 0.1);// - bead_diameter/2;
		Real n_ring = (int) floor(2 * M_PI * ring_radius / (bead_diameter - 0.1));
		for (int i = n_1; i < (n_1 + n_ring) && i < N_s; i++)
		{
			r.x = (ring_radius+(chain_length-1)*Particle::sigma_p/2)*cos(ring_theta);
			r.y = (ring_radius+(chain_length-1)*Particle::sigma_p/2)*sin(ring_theta);
			particle[N_m+i].r = r;
			particle[N_m+i].theta = ring_theta + M_PI*(i%2);
			particle[N_m+i].v.x = cos(particle[N_m+i].theta);
			particle[N_m+i].v.y = sin(particle[N_m+i].theta);
			if (n_ring < (N_s - n_1))
				ring_theta += 2 * M_PI /n_ring;
			else
				ring_theta += 2 * M_PI / (N_s - n_1);
		}
		n_1 += n_ring;
		ring_index += 1;
		ring_theta = 0;
	}
}


void Square_Membrane(Particle* particle, const Real bead_diameter, int N_l)
{/*
	Put all of the membrane beads on a square N_l * N_l
	N_l = Number of membrane beads on each side of the square
	bead_diameter = Particle::sigma_p
*/
	C2DVector r;

	for (int i = 0; i < N_l; i++)
	{
		r.x = (-N_l + 1) * bead_diameter /2 + i * bead_diameter;
		r.y = (-N_l + 1) * bead_diameter /2;
		particle[i].r = r;
	}
	for (int i = N_l; i < 2*N_l-1; i++)
	{
		r.x = (N_l - 1) * bead_diameter /2;
		r.y = (-N_l + 1) * bead_diameter /2 + (i - N_l + 1)* bead_diameter;
		particle[i].r = r;
	}
	for (int i = 2*N_l -1; i < 3*N_l-2; i++)
	{
		r.x = (N_l - 1) * bead_diameter /2 - (i - 2*N_l + 2) * bead_diameter;
		r.y = (N_l - 1) * bead_diameter /2;
		particle[i].r = r;
	}
	for (int i = 3*N_l-2; i < 4*N_l-4; i++)
	{
		r.x = (-N_l + 1) * bead_diameter /2;
		r.y = (N_l - 1) * bead_diameter /2 - (i - 3*N_l + 3) * bead_diameter;
		particle[i].r = r;
	}
}

void Confined_In_Square_Membrane(Particle* particle, const Real bead_diameter,int N_s, int N_l)
{/*
	Put all of the active beads inside the square membrane (rows)
	N_s = Number of swimmers
	N_l = Number of membrane beads on each side of the square
	bead_diameter = Particle::sigma_p
*/
	int chain_length = particle[4*N_l-4].nb; //number of beads in each swimmer(chain)
	if (N_s > (N_l-2) * (int) ((N_l-2)/chain_length) )
	{
		cout << "ERROR: Number of swimmers are exceeding the membrane capacity!!!" << endl;
		exit(0);
	}

	C2DVector r;

	int n_x = (int) ceil (sqrt(N_s*chain_length)); //number of columns
	int n_y = (int) ceil (N_s/n_x);			  //number of rows
	double b = (N_l-2)/n_x;					  //aspect ratio of unit cell
	for (int j = 0; j < n_y+1; j+=1) 
	{
		for (int i = 0; i < n_x+1; i++)
		{
				r.x = 0.5*(-N_l + 2 + b) * bead_diameter + i * b * bead_diameter;
				r.y = 0.5*(-N_l + 2 + chain_length) * bead_diameter + j * b * chain_length * bead_diameter;
				particle[4*N_l-4 + j * n_x + i].r = r;
				particle[4*N_l-4 + j * n_x + i].theta = M_PI*(i%2) + M_PI*(j%2) + M_PI/2;
				particle[4*N_l-4 + j * n_x + i].v.x = cos(particle[4*N_l-4+i].theta);
				particle[4*N_l-4 + j * n_x + i].v.y = sin(particle[4*N_l-4+i].theta);
		}
	}
}


#endif
