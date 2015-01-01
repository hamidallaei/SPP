#ifndef _BOX_
#define _BOX_

#include "../shared/parameters.h"
#include "../shared/c2dvector.h"
#include "../shared/particle.h"
#include "../shared/cell.h"
#include "../shared/wall.h"
#include "node.h"


class Box{
public:
	int N, wall_num; // N is the number of particles and wallnum is the number of walls in the system.
	Particle particle[max_N]; // Array of particles that we are going to simulate.
	Wall wall[8]; // Array of walls in our system.

	Real density, volume_fraction;
	stringstream info; // information stream that contains the simulation information, like noise, density and etc. this will be used for the saving name of the system.

	Node* thisnode; // Node is a class that has information about the node_id and its boundaries, neighbores and etc.

	Box();
	void Init_Topology(); // Initialize the wall positions and numbers.
	void Init(Node* input_node,Real input_volume_fraction, Real g = 2, Real alpha=0.5, Real noise_amplitude = 0.1); // Intialize the box, positioning particles, giving them velocities, updating cells and sending information to all nodes.
	void Square_Lattice_Formation(); // Positioning partilces in a square lattice
	void Triangle_Lattice_Formation(); // Positioning partilces in a triangular lattice. This is denser.
	void Single_Vortex_Formation(); // Vortex initial condition.
	void Four_Vortex_Formation(); // Four vortex inside the box. left top, left bot, right top and right bot.
	void Clump_Formation(int size); // Positioning particles in a clump that is moving in some direction.
	void Interact(); // Here the intractio of particles are computed that is the applied tourque to each particle.
	void Move(); // Move all particles of this node.
	void One_Step(); // One full step, composed of interaction computation and move.
	void Multi_Step(int steps); // Several steps befor a cell upgrade.
	void Translate(C2DVector d); // Translate position of all particles with vector d

	friend std::ostream& operator<<(std::ostream& os, Box* box); // Save
	friend std::istream& operator>>(std::istream& is, Box* box); // Input
};

Box::Box()
{
	N = 0;
	density = 0;
	wall_num = 0;
}

// Initialize the wall positions and numbers.
void Box::Init_Topology()
{
	thisnode->Get_Box_Info(N,particle);
	thisnode->Init_Topology();
	#ifndef PERIODIC_BOUNDARY_CONDITION
		wall_num = 4;
		wall[0].Init(-Lx,-Ly,-Lx,Ly);
		wall[1].Init(-Lx,Ly,Lx,Ly);
		wall[2].Init(Lx,Ly,Lx,-Ly);
		wall[3].Init(Lx,-Ly,-Lx,-Ly);
	#endif
}

// Intialize the box, positioning particles, giving them velocities, updating cells and sending information to all nodes.
void Box::Init(Node* input_node, Real input_density, Real g, Real alpha, Real noise_amplitude)
{
	#ifdef TRACK_PARTICLE
	track_p = &particle[track];
	#endif

	thisnode = input_node;
	density = input_density;
	N = (int) round(L2*L2*density);

	Particle::noise_amplitude = noise_amplitude / sqrt(dt); // noise amplitude depends on the step (dt) because of ito calculation. If we have epsilon in our differential equation and we descritise it with time steps dt, the noise in each step that we add is epsilon times sqrt(dt) if we factorise it with a dt we have dt*(epsilon/sqrt(dt)).
	ContinuousParticle::g = g;
	ContinuousParticle::alpha = alpha;

	Init_Topology(); // Adding walls

	if (thisnode->node_id == 0)
	{
		cout << "number_of_particles = " << N << endl; // Printing number of particles.
// Positioning the particles
		Triangle_Lattice_Formation();
	//	Single_Vortex_Formation();
	//	Four_Vortex_Formation();
	}
	MPI_Barrier(MPI_COMM_WORLD);

// Master node will broadcast the particles information
	thisnode->Root_Bcast();
// Any node update cells, knowing particles and their cell that they are inside.
	thisnode->Full_Update_Cells();

// Buliding up info stream. In next versions we will take this part out of box, making our libraries more abstract for any simulation of SPP.
	info.str("");
	info << "rho=" << density <<  "-g=" << Particle::g << "-alpha=" << Particle::alpha << "-noise=" << noise_amplitude;
	MPI_Barrier(MPI_COMM_WORLD);
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
		if ((r.y > Ly) || (r.y < -Ly) || (r.x > Lx) || (r.x < -Lx))
			cout << r/L << endl;
//		r.Periodic_Transform();
		particle[i].Init(r);
		v_cm += (particle[i].v / N);
	}
	for (int i = 0; i < N; i++)
		particle[i].v -= v_cm; // There are issues about this mass center motion removal. Angle of partilces must be changed not their velocities.
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



// Here the intractio of particles are computed that is the applied tourque to each particle.
void Box::Interact()
{
// each node sends its information of boundary cells to the correspounding node.
	#ifdef PERIODIC_BOUNDARY_CONDITION
	for (int i = 0; i < thisnode->boundary.size(); i++)
	{
		thisnode->boundary[i].Send_Data();
		thisnode->boundary[(i+4)%8].Receive_Data();
	}
	#else
		for (int i = 0; i < thisnode->boundary.size(); i++)
			thisnode->boundary[i].Send_Data(); // Send information of i'th boundary of thisnode to the neighboring node that shares this boundary.
		for (int i = 0; i < thisnode->boundary.size(); i++)
			thisnode->boundary[i].Receive_Data(); // Receive information of i'th boundary of neighboring node.
// Sum up interaction of the walls with the particles of thisnode (in the absence of periodic boundary condition)
		for (int i = 0; i < wall_num; i++)
			for (int x = thisnode->head_cell_idx; x < thisnode->tail_cell_idx; x++)
				for (int y = thisnode->head_cell_idy; y < thisnode->tail_cell_idy; y++)
					for (int j = 0; j < thisnode->cell[x][y].pid.size(); j++)
						wall[i].Interact(&particle[thisnode->cell[x][y].pid[j]]);
	#endif

	thisnode->Self_Interact(); // Sum up interaction of particles within thisnode
	thisnode->Boundary_Interact(); // Sum up interaction of particles in the neighboring nodes.
}

// Move all particles of this node.
void Box::Move()
{
	thisnode->Move();
}

// One full step, composed of interaction computation and move.
void Box::One_Step()
{
	Interact();
	Move();
	MPI_Barrier(MPI_COMM_WORLD);
}

// Several steps befor a cell upgrade.
void Box::Multi_Step(int steps)
{
	for (int i = 0; i < steps; i++)
	{
		Interact();
		Move();
		MPI_Barrier(MPI_COMM_WORLD); // Barier guranty that the move step of all particles is done. Therefor in interact function we are using updated particles.
	}
	thisnode->Quick_Update_Cells();
}

// Translate position of all particles with vector d
void Box::Translate(C2DVector d)
{
	thisnode->Root_Gather();
	for (int i = 0; i < N; i++)
	{
		particle[i].r += d;
		particle[i].r.Periodic_Transform();
	}
	thisnode->Root_Bcast();
	thisnode->Full_Update_Cells();
}

// Saving the particle information (position and velocities) to a standard output stream (probably a file). This must be called by only the root.
std::ostream& operator<<(std::ostream& os, Box* box)
{
	box->thisnode->Root_Gather();
	box->thisnode->Root_Bcast();

	if (box->thisnode->node_id == 0)
	{
		os.write((char*) &box->N, sizeof(box->N) / sizeof(char));
		for (int i = 0; i < box->N; i++)
		{
			box->particle[i].r.write(os);
			C2DVector v;
			v.x = cos(box->particle[i].theta);
			v.y = sin(box->particle[i].theta);
			v.write(os);

//			cout << i << "\t" << setprecision(100) << box->particle[i].theta << endl;
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
}

// Reading the particle information (position and velocities) from a standard input stream (probably a file).
std::istream& operator>>(std::istream& is, Box* box)
{
	if (box->thisnode->node_id == 0)
	{
		is.read((char*) &box->N, sizeof(int) / sizeof(char));
		for (int i = 0; i < box->N; i++)
		{
			is >> box->particle[i].r;
			is >> box->particle[i].v;
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
}

#endif
