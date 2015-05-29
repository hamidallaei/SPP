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
	void Init(Node* input_node, Real input_density); // Intialize the box, positioning particles, giving them velocities, updating cells and sending information to all nodes.

	void Load(const State_Hyper_Vector&);
	void Save(State_Hyper_Vector&);
	
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
void Box::Init(Node* input_node, Real input_density)
{
	#ifdef TRACK_PARTICLE
	track_p = &particle[track];
	#endif

	thisnode = input_node;
	density = input_density;
	N = (int) round(Lx2*Ly2*density);

	Init_Topology(); // Adding walls

	if (thisnode->node_id == 0)
	{
		cout << "number_of_particles = " << N << endl; // Printing number of particles.
// Positioning the particles
//		Triangle_Lattice_Formation(particle, N, 1);
//		Random_Formation(particle, N, 0); // Positioning partilces Randomly, but distant from walls
		Random_Formation_Circle(particle, N, Lx-1); // Positioning partilces Randomly, but distant from walls
//		Single_Vortex_Formation(particle, N);
	//	Four_Vortex_Formation(particle, N);
	}
	MPI_Barrier(MPI_COMM_WORLD);

// Master node will broadcast the particles information
	thisnode->Root_Bcast();
// Any node update cells, knowing particles and their cell that they are inside.
	thisnode->Full_Update_Cells();

	#ifdef verlet_list
	thisnode->Update_Neighbor_List();
	#endif

// Buliding up info stream. In next versions we will take this part out of box, making our libraries more abstract for any simulation of SPP.
	info.str("");
	MPI_Barrier(MPI_COMM_WORLD);
}

// Loading a state to the box.
void Box::Load(const State_Hyper_Vector& sv)
{
	if (N != sv.N)
	{
		cout << "Error: Number of particles in state vectors differ from box" << endl;
		exit();
	}
	for (int i = 0; i < N; i++)
	{
		particle[i].r = sv.particle[i].r;
		particle[i].theta = sv.particle[i].theta;
		particle[i].v.x = cos(particle[i].theta);
		particle[i].v.y = sin(particle[i].theta);
	}
	sv.Set_C2DVector_Rand_Generator();
}

// Saving state of the box.
void Box::Save(State_Hyper_Vector& sv)
{
	if (N != sv.N)
	{
		cout << "Error: Number of particles in state vectors differ from box" << endl;
		exit();
	}
	for (int i = 0; i < N; i++)
	{
		sv.particle[i].r = particle[i].r;
		sv.particle[i].theta = particle[i].theta;
	}
	sv.Get_C2DVector_Rand_Generator();
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

	#ifndef PERIODIC_BOUNDARY_CONDITION
		#ifdef CIRCULAR_BOX
// Sum up interaction of the circualr boundary with the particles of thisnode (in the absence of periodic boundary condition)
		for (int x = thisnode->head_cell_idx; x < thisnode->tail_cell_idx; x++)
			for (int y = thisnode->head_cell_idy; y < thisnode->tail_cell_idy; y++)
				for (int j = 0; j < thisnode->cell[x][y].pid.size(); j++)
				{
					int i = thisnode->cell[x][y].pid[j];
					Real r = sqrt(particle[i].r.Square());
					if (r > (Lx-1))
						particle[i].torque += 40*Particle::g*(particle[i].v.y*particle[i].r.x - particle[i].v.x*particle[i].r.y) / (2*M_PI*r*(Lx-r));
				}
		#else
// Sum up interaction of the walls with the particles of thisnode (in the absence of periodic boundary condition)
		for (int i = 0; i < wall_num; i++)
			for (int x = thisnode->head_cell_idx; x < thisnode->tail_cell_idx; x++)
				for (int y = thisnode->head_cell_idy; y < thisnode->tail_cell_idy; y++)
					for (int j = 0; j < thisnode->cell[x][y].pid.size(); j++)
						wall[i].Interact(&particle[thisnode->cell[x][y].pid[j]]);
		#endif
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
	#ifdef verlet_list
	thisnode->Update_Neighbor_List();
	#endif
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
