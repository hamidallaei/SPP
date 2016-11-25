#ifndef _BOX_
#define _BOX_

#include "../shared/parameters.h"
#include "../shared/c2dvector.h"
#include "../shared/particle.h"
#include "../shared/cell.h"
#include "../shared/wall.h"
#include "../shared/set-up.h"
#include "../shared/state-hyper-vector.h"
#include "node.h"

#include <boost/algorithm/string.hpp>

class Box{
public:
	int N, Nb, Np, Nw; // N is the number of particles, Nb is number of beads, Np is number of active particles and Nw is the number of 
	Particle* particle; // Array of particles that we are going to simulate.

	Real t;
	Real packing_fraction;

	stringstream info; // information stream that contains the simulation information, like noise, density and etc. this will be used for the saving name of the system.

	Node* thisnode; // Node is a class that has information about the node_id and its boundaries, neighbores and etc.

	Box();
	// I believe that it is better to move these init functions to main files
	void Init_Topology(); // Initialize the wall positions and numbers.
	void Init(Node* input_node, int, int); // Intialize the box, positioning particles, giving them velocities, updating cells and sending information to all nodes.
	bool Positioning_Particles(Node* input_node, const string name); // Intialize the box from a file, this includes reading particles information, updating cells and sending information to all nodes.
	void Sync();

	void Interact_Wall_Beads(); // Here the wall particles interact via a spring
	void Interact(); // Here the intractio of particles are computed that is the applied tourque to each particle.
	void Move(); // Move all particles of this node.
	void One_Step(); // One full step, composed of interaction computation and move.
	void Multi_Step(int steps); // Several steps befor a cell upgrade.
	void Multi_Step(int steps, int interval); // Several steps with a cell upgrade call after each interval.
	void Translate(C2DVector d); // Translate position of all particles with vector d
	void Save_Polarization(std::ostream& os); // Save polarization of the particles inside the box

	friend std::ostream& operator<<(std::ostream& os, Box* box); // Save
	friend std::istream& operator>>(std::istream& is, Box* box); // Input
};

Box::Box()
{
	N = 0;
	particle = new Particle[max_N];
}

// Initialize the wall positions and numbers.
void Box::Init_Topology()
{
	thisnode->Get_Box_Info(N,particle);
	thisnode->Init_Topology();
}

// Sync positions of the particles with other nodes
void Box::Sync()
{
	MPI_Barrier(MPI_COMM_WORLD);
	thisnode->Root_Bcast();
// Any node update cells, knowing particles and their cell that they are inside.
	thisnode->Full_Update_Cells();

	#ifdef verlet_list
	thisnode->Update_Neighbor_List();
	#endif

	MPI_Barrier(MPI_COMM_WORLD);
}


// Intialize the box, positioning particles, giving them velocities, updating cells and sending information to all nodes.
void Box::Init(Node* input_node, const int input_Np, const int input_Nw)
{
	#ifdef TRACK_PARTICLE
		track_p = &particle[track];
	#endif

	thisnode = input_node;

	Np = input_Np;
	Nw = input_Nw;
	N = Np + Nw;

	Nb = Nw;
	for (int i = Nw; i < N; i++)
	{
			Nb += particle[i].nb;
	}

	packing_fraction=0.5;

	Init_Topology(); // Adding walls

// Positioning the particles at first time. Note that, the positions can be tunned in the main file as well
	if (thisnode->node_id == 0)
	{
		cout << "number_of_particles = " << N << endl; // Printing number of particles.
		cout << "number_of_beads = " << Nb << endl; // Printing number of particles.
//		Triangle_Lattice_Formation(particle, N, 1);
	}
	Sync();

// Buliding up info stream. In next versions we will take this part out of box, making our libraries more abstract for any simulation of SPP.
	info.str("");
	MPI_Barrier(MPI_COMM_WORLD);
}


// Intialize the box from a file, this includes reading particles information, updating cells and sending information to all nodes. 
bool Box::Positioning_Particles(Node* input_node, const string input_name)
{
	#ifdef TRACK_PARTICLE
		track_p = &particle[track];
	#endif

	thisnode = input_node;

	string name = input_name;
	stringstream address(name);
	ifstream is;
	is.open(address.str().c_str());
	if (!is.is_open())
		return false;

	if (thisnode->node_id == 0)
	{
		is.read((char*) &N, sizeof(int) / sizeof(char));
		if (N < 0 || N > 1000000)
			return (false);
		while (!is.eof())
		{
			for (int i = 0; i < N; i++)
			{
				is >> particle[i].r;
				is >> particle[i].v;
			}
			is.read((char*) &N, sizeof(int) / sizeof(char));
		}
		for (int i = 0; i < N; i++)
			particle[i].theta = atan2(particle[i].v.y,particle[i].v.x);
		cout << "number_of_particles = " << N << endl; // Printing number of particles.
	}
	else
	{
		is.read((char*) &N, sizeof(int) / sizeof(char));
		if (N < 0 || N > 1000000)
			return (false);
	}

	is.close();

	Init_Topology(); // Adding walls

	Sync();

// Buliding up info stream. In next versions we will take this part out of box, making our libraries more abstract for any simulation of SPP.
	info.str("");
	MPI_Barrier(MPI_COMM_WORLD);

	return (true);
}


void Box::Interact_Wall_Beads()
{
	for (int i = 0; i < Nw; i++)
	{
		C2DVector dr = particle[i].r - particle[(i+1)%Nw].r;
		dr.Periodic_Transform();
		Real d = sqrt(dr.Square());
		C2DVector f = Spring(dr, d, Particle::sigma_p, 500);
		particle[i].f += f;
		particle[(i+1)%Nw].f -= f;
	}
}

// Here the intractio of particles are computed that is the applied tourque to each particle.
void Box::Interact()
{
	thisnode->Send_Receive_Data();
	MPI_Barrier(MPI_COMM_WORLD);

	#ifdef verlet_list
// with verlet list:
	thisnode->Neighbor_List_Interact();
	#else
// without verlet list:
	thisnode->Self_Interact(); // Sum up interaction of particles within thisnode
	thisnode->Boundary_Interact(); // Sum up interaction of particles in the neighboring nodes.
	#endif
}

// Move all particles of this node.
void Box::Move()
{
	thisnode->Move();
}

// One full step, composed of interaction computation and move.
void Box::One_Step()
{
	Interact_Wall_Beads();
	Interact();
	Move();
	t += dt;
	MPI_Barrier(MPI_COMM_WORLD);
}

// Several steps befor a cell upgrade.
void Box::Multi_Step(int steps)
{
	for (int i = 0; i < steps; i++)
	{
		Interact_Wall_Beads();
		Interact();
		Move();
		MPI_Barrier(MPI_COMM_WORLD); // Barier guranty that the move step of all particles is done. Therefor in interact function we are using updated particles.
	}
	t += dt*steps;
//	cout << "Updating cells" << endl;
	thisnode->Quick_Update_Cells();
//	cout << "Updating finished" << endl;
	#ifdef verlet_list
	thisnode->Update_Neighbor_List();
	#endif
}

// Several steps with a cell upgrade call after each interval. Warning, I see no cell update function call! I have to fix it!
void Box::Multi_Step(int steps, int interval)
{
	for (int i = 0; i < steps/interval; i++)
		Multi_Step(interval);
	Multi_Step(steps % interval);
	t += interval;
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
	#ifdef verlet_list
	thisnode->Update_Neighbor_List();
	#endif
}


// Saving the particle information (position and velocities) to a standard output stream (probably a file). This must be called by only the root.
std::ostream& operator<<(std::ostream& os, Box* box)
{
	box->thisnode->Root_Gather();
	box->thisnode->Root_Bcast();

	if (box->thisnode->node_id == 0)
	{
		os.write((char*) &box->Nb, sizeof(box->Nb) / sizeof(char));
		for (int i = 0; i < box->N; i++)
		{
			box->particle[i].Write(os);
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
