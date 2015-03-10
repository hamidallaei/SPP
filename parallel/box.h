#ifndef _BOX_
#define _BOX_

#include "../shared/parameters.h"
#include "../shared/c2dvector.h"
#include "../shared/particle.h"
#include "../shared/cell.h"
#include "../shared/wall.h"
#include "../shared/set-up.h"
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
		wall[0].Init(-Lx,-Ly,-Lx, Ly);
		wall[1].Init(-Lx, Ly, Lx, Ly);
		wall[2].Init( Lx, Ly, Lx,-Ly);
		wall[3].Init( Lx,-Ly,-Lx,-Ly);
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
	N = (int) round(Lx2*Ly2*density);

	Particle::noise_amplitude = noise_amplitude / sqrt(dt); // noise amplitude depends on the step (dt) because of ito calculation. If we have epsilon in our differential equation and we descritise it with time steps dt, the noise in each step that we add is epsilon times sqrt(dt) if we factorise it with a dt we have dt*(epsilon/sqrt(dt)).
	ContinuousParticle::g = g;
	ContinuousParticle::alpha = alpha;

	Init_Topology(); // Adding walls

	if (thisnode->node_id == 0)
	{
		cout << "number_of_particles = " << N << endl; // Printing number of particles.
// Positioning the particles
//		Triangle_Lattice_Formation(particle, N, 1);
		Single_Vortex_Formation(particle, N);
	//	Four_Vortex_Formation(particle, N);
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

// Here the intractio of particles are computed that is the applied tourque to each particle.
void Box::Interact()
{
// each node sends its information of boundary cells to the correspounding node.
	for (int i = 0; i < thisnode->boundary.size(); i++)
	{
		if (thisnode->boundary[i].is_active)
			thisnode->boundary[i].Send_Data(); // Send information of i'th boundary of thisnode to the neighboring node that shares this boundary (if there is any).
		if (thisnode->boundary[(i+4)%8].is_active)
			thisnode->boundary[(i+4)%8].Receive_Data(); // Receive information of i'th boundary of neighboring node (if the neighbor exitst).
	}

	thisnode->Self_Interact(); // Sum up interaction of particles within thisnode
	thisnode->Boundary_Interact(); // Sum up interaction of particles in the neighboring nodes.

	#ifndef PERIODIC_BOUNDARY_CONDITION
// Sum up interaction of the walls with the particles of thisnode (in the absence of periodic boundary condition)
		for (int i = 0; i < wall_num; i++)
			for (int x = thisnode->head_cell_idx; x < thisnode->tail_cell_idx; x++)
				for (int y = thisnode->head_cell_idy; y < thisnode->tail_cell_idy; y++)
					for (int j = 0; j < thisnode->cell[x][y].pid.size(); j++)
						wall[i].Interact(&particle[thisnode->cell[x][y].pid[j]]);
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
