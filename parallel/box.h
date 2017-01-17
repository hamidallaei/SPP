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
	int Ns, Nm, wall_num; // N is the number of particles and wallnum is the number of walls in the system.
	Particle* particle; // Array of particles that we are going to simulate.
	Wall wall[8]; // Array of walls in our system.

	Real t;
	Real density, packing_fraction;
	Real polarization;
	stringstream info; // information stream that contains the simulation information, like noise, density and etc. this will be used for the saving name of the system.

	Node* thisnode; // Node is a class that has information about the node_id and its boundaries, neighbores and etc.

	Box();
	// I believe that it is better to move these init functions to main files
	void Init_Topology(); // Initialize the wall positions and numbers.
	void Init(Node* input_node, Real input_density); // Intialize the box, positioning particles, giving them velocities, updating cells and sending information to all nodes.
	bool Positioning_Particles(Node* input_node, const string name); // Intialize the box from a file, this includes reading particles information, updating cells and sending information to all nodes.
	void Sync();

	void Load(const State_Hyper_Vector&); // Load new position and angles of particles and a gsl random generator from a state hyper vector
	void Save(State_Hyper_Vector&) const; // Save current position and angles of particles and a gsl random generator to a state hyper vector

	void Interact(); // Here the intractio of particles are computed that is the applied tourque to each particle.
	void Move(); // Move all particles of this node.
	void Move_Runge_Kutta2_1(); // Do the first step of second order Runge Kutta
	void Move_Runge_Kutta2_2(); // Do the second step of second order Runge Kutta

	void Move_Runge_Kutta4_1(); // Do the first step of forth order Runge Kutta
	void Move_Runge_Kutta4_2(); // Do the second step of forth order Runge Kutta
	void Move_Runge_Kutta4_3(); // Do the third step of forth order Runge Kutta
	void Move_Runge_Kutta4_4(); // Do the forth step of forth order Runge Kutta

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
	Ns = 0;
	Nm = 0;
	density = 0;
	wall_num = 0;
	particle = new Particle[max_N];
}

// Initialize the wall positions and numbers.
void Box::Init_Topology()
{
	thisnode->Get_Box_Info(Ns,particle);
	thisnode->Init_Topology();
	#ifndef PERIODIC_BOUNDARY_CONDITION
		wall_num = 4;
		wall[0].Init(-Lx,-Ly,-Lx, Ly);
		wall[1].Init(-Lx, Ly, Lx, Ly);
		wall[2].Init( Lx, Ly, Lx,-Ly);
		wall[3].Init( Lx,-Ly,-Lx,-Ly);
	#endif
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
void Box::Init(Node* input_node, Real input_density)
{
	#ifdef TRACK_PARTICLE
	track_p = &particle[track];
	#endif

	thisnode = input_node;

	density = input_density;
	Ns = (int) round(Lx2*Ly2*input_density);
	Nm = 0;

	Init_Topology(); // Adding walls

// Positioning the particles at first time. Note that, the positions can be tunned in the main file as well
	if (thisnode->node_id == 0)
	{
		cout << "number_of_particles = " << Ns << endl; // Printing number of particles.
//		Triangle_Lattice_Formation(particle, Ns, 1);
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
		is.read((char*) &Ns, sizeof(int) / sizeof(char));
		if (Ns < 0 || Ns > 1000000)
			return (false);
		while (!is.eof())
		{
			for (int i = 0; i < Ns; i++)
			{
				is >> particle[i].r;
				is >> particle[i].v;
			}
			is.read((char*) &Ns, sizeof(int) / sizeof(char));
		}
		for (int i = 0; i < Ns; i++)
			particle[i].theta = atan2(particle[i].v.y,particle[i].v.x);
		cout << "number_of_particles = " << Ns << endl; // Printing number of particles.
	}
	else
	{
		is.read((char*) &Ns, sizeof(int) / sizeof(char));
		if (Ns < 0 || Ns > 1000000)
			return (false);
	}

	is.close();

	density = Ns / (Lx2*Ly2);

	Init_Topology(); // Adding walls

	Sync();

// Buliding up info stream. In next versions we will take this part out of box, making our libraries more abstract for any simulation of SPP.
	info.str("");
	MPI_Barrier(MPI_COMM_WORLD);

	return (true);
}


// Loading a state to the box.
void Box::Load(const State_Hyper_Vector& sv)
{
	if (Ns != sv.N)
	{
		cout << "Error: Number of particles in state vectors differ from box" << endl;
		exit(0);
	}
	for (int i = 0; i < Ns; i++)
	{
		particle[i].r = sv.particle[i].r;
		particle[i].theta = sv.particle[i].theta;
		particle[i].v.x = cos(particle[i].theta);
		particle[i].v.y = sin(particle[i].theta);
	}
	sv.Set_C2DVector_Rand_Generator();
	MPI_Barrier(MPI_COMM_WORLD);
	thisnode->Root_Bcast();
	thisnode->Full_Update_Cells();
	#ifdef verlet_list
		thisnode->Update_Neighbor_List();
	#endif
}

// Saving state of the box.
void Box::Save(State_Hyper_Vector& sv) const
{
	if (Ns != sv.N)
	{
		cout << "Error: Number of particles in state vectors differ from box" << endl;
		exit(0);
	}
	thisnode->Root_Gather();
	thisnode->Root_Bcast();
	for (int i = 0; i < Ns; i++)
	{
		sv.particle[i].r = particle[i].r;
		sv.particle[i].theta = particle[i].theta;
	}
	sv.Get_C2DVector_Rand_Generator();
// We need to make sure that indexing of particles are the same to exactly recompute the same values. Therefor at a saving we update cells and neighore list therefore if we load the same sv and update cells and neighore list we will come to the same indexing
	MPI_Barrier(MPI_COMM_WORLD);
	thisnode->Full_Update_Cells();
	#ifdef verlet_list
		thisnode->Update_Neighbor_List();
	#endif
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

#ifdef RUNGE_KUTTA2
// Do the first step of second order Runge Kutta
void Box::Move_Runge_Kutta2_1()
{
	thisnode->Move_Runge_Kutta2_1();
}

// Do the second step of second order Runge Kutta
void Box::Move_Runge_Kutta2_2()
{
	thisnode->Move_Runge_Kutta2_2();
}
#endif

#ifdef RUNGE_KUTTA4
// Do the first step of forth order Runge Kutta
void Box::Move_Runge_Kutta4_1()
{
	thisnode->Move_Runge_Kutta4_1();
}

// Do the second step of forth Runge Kutta
void Box::Move_Runge_Kutta4_2()
{
	thisnode->Move_Runge_Kutta4_2();
}

// Do the third step of forth Runge Kutta
void Box::Move_Runge_Kutta4_3()
{
	thisnode->Move_Runge_Kutta4_3();
}

// Do the forth step of forth Runge Kutta
void Box::Move_Runge_Kutta4_4()
{
	thisnode->Move_Runge_Kutta4_4();
}
#endif

// One full step, composed of interaction computation and move.
void Box::One_Step()
{
	#ifdef RUNGE_KUTTA2
		Interact();
		Move_Runge_Kutta2_1();

		MPI_Barrier(MPI_COMM_WORLD);

		Interact();
		Move_Runge_Kutta2_2();
		MPI_Barrier(MPI_COMM_WORLD);
	#else
		#ifdef RUNGE_KUTTA4
			Interact();
			Move_Runge_Kutta4_1();
			MPI_Barrier(MPI_COMM_WORLD);

			Interact();
			Move_Runge_Kutta4_2();
			MPI_Barrier(MPI_COMM_WORLD);

			Interact();
			Move_Runge_Kutta4_3();
			MPI_Barrier(MPI_COMM_WORLD);

			Interact();
			Move_Runge_Kutta4_4();
			MPI_Barrier(MPI_COMM_WORLD);
		#else
			Interact();
			Move();
			MPI_Barrier(MPI_COMM_WORLD);
		#endif
	#endif

	t += dt;
	MPI_Barrier(MPI_COMM_WORLD);
}

// Several steps befor a cell upgrade.
void Box::Multi_Step(int steps)
{
	for (int i = 0; i < steps; i++)
	{
		#ifdef RUNGE_KUTTA2
			Interact();
			Move_Runge_Kutta2_1();

			MPI_Barrier(MPI_COMM_WORLD);

			Interact();
			Move_Runge_Kutta2_2();
			MPI_Barrier(MPI_COMM_WORLD);
		#else
			#ifdef RUNGE_KUTTA4
				Interact();
				Move_Runge_Kutta4_1();
				MPI_Barrier(MPI_COMM_WORLD);

				Interact();
				Move_Runge_Kutta4_2();
				MPI_Barrier(MPI_COMM_WORLD);

				Interact();
				Move_Runge_Kutta4_3();
				MPI_Barrier(MPI_COMM_WORLD);

				Interact();
				Move_Runge_Kutta4_4();
				MPI_Barrier(MPI_COMM_WORLD);
			#else
				Interact();
				Move();
				MPI_Barrier(MPI_COMM_WORLD); // Barier guranty that the move step of all particles is done. Therefor in interact function we are using updated particles.
			#endif
		#endif
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
	for (int i = 0; i < Ns; i++)
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

// Save polarization of the particles inside the box
void Box::Save_Polarization(std::ostream& os)
{
	thisnode->Compute_Polarization();
	polarization = thisnode->polarization;
	if (thisnode->node_id == 0)
		os << t << "\t" << polarization << endl;
}

// Saving the particle information (position and velocities) to a standard output stream (probably a file). This must be called by only the root.
std::ostream& operator<<(std::ostream& os, Box* box)
{
	box->thisnode->Root_Gather();
	box->thisnode->Root_Bcast();

	if (box->thisnode->node_id == 0)
	{
		os.write((char*) &(box->t), sizeof(Real) / sizeof(char));
		os.write((char*) &(Lx), sizeof(Real) / sizeof(char));
		os.write((char*) &(Ly), sizeof(Real) / sizeof(char));
		os.write((char*) &(box->particle[box->Nm].nb), sizeof(int) / sizeof(char));
		os.write((char*) &(box->Ns), sizeof(box->Ns) / sizeof(char));
		os.write((char*) &(box->Nm), sizeof(box->Nm) / sizeof(char));

		for (int i = box->Nm; i < box->Ns; i++)
		{
			box->particle[i].Write(os);
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
}

// Reading the particle information (position and velocities) from a standard input stream (probably a file).
std::istream& operator>>(std::istream& is, Box* box)
{
	if (box->thisnode->node_id == 0)
	{
		is.read((char*) &box->Ns, sizeof(int) / sizeof(char));
		for (int i = 0; i < box->Ns; i++)
		{
			is >> box->particle[i].r;
			is >> box->particle[i].v;
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
}

#endif
