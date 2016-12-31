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
	int N, Ns, Nm; // N = Ns+Nm is the number of particles, Nm is number of membrane particles, Ns is number of active particles
	Particle* particle; // Array of particles that we are going to simulate.

	Real t;
	Real packing_fraction;
	Real membrane_elasticity;
	Real membrane_radius;

	Real t_old; // old time, the previous time that coordinates are saved for further analysis
	C2DVector* rm_old; // old position of membrane beads
	C2DVector r_cm, v_cm;
	Real xx,yy,xy; // Qxx, Qyy, Qxy. Q = 1/N sum (\vec{r} - \vec{r}_cm) (\vec{r} - \vec{r}_cm)
	Real l; // Angular momentum per particle
	Real omega; // Angular velocity
	Real lambda1, lambda2; // Eigenvalues of Q tensor
	Real Rg; // Gyration radius
	Real Delta; // Asphericity

	stringstream info; // information stream that contains the simulation information, like noise, density and etc. this will be used for the saving name of the system.

	Node* thisnode; // Node is a class that has information about the node_id and its boundaries, neighbores and etc.

	Box();
	// I believe that it is better to move these init functions to main files
	void Init_Topology(); // Initialize the wall positions and numbers.
	void Init(Node* input_node, int, int); // Intialize the box, positioning particles, giving them velocities, updating cells and sending information to all nodes.
	bool Positioning_Particles(Node* input_node, const string name); // Intialize the box from a file, this includes reading particles information, updating cells and sending information to all nodes.
	void Sync();

	void Interact_Membrane_Beads(); // Here the wall particles interact via a spring
	void Interact(); // Here the intractio of particles are computed that is the applied tourque to each particle.
	void Move(); // Move all particles of this node.
	void Move_Runge_Kutta_1(); // Do the first step of Runge Kutta
	void Move_Runge_Kutta_2(); // Do the second step of Runge Kutta
	void One_Step(); // One full step, composed of interaction computation and move.
	void Multi_Step(int steps); // Several steps befor a cell upgrade.
	void Multi_Step(int steps, int interval); // Several steps with a cell upgrade call after each interval.
	void Translate(C2DVector d); // Translate position of all particles with vector d

	void RootGather(); // Gather all needed data for computing variables such as center of mass, angular momentum, angular velocity, and ...
	void Save_Membrane_Position(); // Save the position of membrane beads to compute their velocities in future.
	void Compute_All_Variables(); // Compute center of mass position and speed, I, ...
	void Save_All_Variables(std::ostream& os); // Save center of mass position and speed, I, ...

	friend std::ostream& operator<<(std::ostream& os, Box* box); // Save
	friend std::istream& operator>>(std::istream& is, Box* box); // Input
};

Box::Box()
{
	N = 0;
	particle = new Particle[max_N];
	rm_old = new C2DVector[5000];
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

	for (int i = 0; i < N; i++)
		particle[i].r_original = particle[i].r;

	MPI_Barrier(MPI_COMM_WORLD);
}


// Intialize the box, positioning particles, giving them velocities, updating cells and sending information to all nodes.
void Box::Init(Node* input_node, const int input_Ns, const int input_Nm)
{
	#ifdef TRACK_PARTICLE
		track_p = &particle[track];
	#endif

	thisnode = input_node;

	Ns = input_Ns;
	Nm = input_Nm;
	N = Ns + Nm;

	Init_Topology(); // Adding walls

// Positioning the particles at first time. Note that, the positions can be tunned in the main file as well
	if (thisnode->node_id == 0)
	{
		cout << "number_of_swimmer = " << Ns << endl; // Printing number of particles.
		cout << "number_of_membrane = " << Nm << endl; // Printing number of particles.
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


void Box::Interact_Membrane_Beads()
{
	for (int i = 0; i < Nm; i++)
	{
		C2DVector dr = particle[i].r - particle[(i+1)%Nm].r;
		dr.Periodic_Transform();
		Real d = sqrt(dr.Square());
		C2DVector f = Spring(dr, d, Particle::sigma_p, membrane_elasticity);
		particle[i].f += f;
		particle[(i+1)%Nm].f -= f;
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

// Do the first step of Runge Kutta
void Box::Move_Runge_Kutta_1()
{
	thisnode->Move_Runge_Kutta_1();
}

// Do the second step of Runge Kutta
void Box::Move_Runge_Kutta_2()
{
	thisnode->Move_Runge_Kutta_2();
}

// One full step, composed of interaction computation and move.
void Box::One_Step()
{
	#ifndef RUNGE_KUTTA
		Interact_Membrane_Beads();
		Interact();
		Move();
		MPI_Barrier(MPI_COMM_WORLD);
	#else
		Interact_Membrane_Beads();
		Interact();
		Move_Runge_Kutta_1();

		MPI_Barrier(MPI_COMM_WORLD);

		Interact_Membrane_Beads();
		Interact();
		Move_Runge_Kutta_2();
		MPI_Barrier(MPI_COMM_WORLD);
	#endif
	t += dt;
}

// Several steps befor a cell upgrade.
void Box::Multi_Step(int steps)
{
	for (int i = 0; i < steps; i++)
	{
		#ifndef RUNGE_KUTTA
			Interact_Membrane_Beads();
			Interact();
			Move();
			MPI_Barrier(MPI_COMM_WORLD); // Barier guranty that the move step of all particles is done. Therefor in interact function we are using updated particles.
		#else		
			Interact_Membrane_Beads();
			Interact();
			Move_Runge_Kutta_1();

			MPI_Barrier(MPI_COMM_WORLD);

			Interact_Membrane_Beads();
			Interact();
			Move_Runge_Kutta_2();

			MPI_Barrier(MPI_COMM_WORLD); // Barier guranty that the move step of all particles is done. Therefor in interact function we are using updated particles.
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

void Box::Save_Membrane_Position()
{
	t_old = t;
	for (int i = 0; i < Nm; i++)
		rm_old[i] = particle[i].r_original;
}

void Box::Compute_All_Variables() // Compute center of mass position and speed, I, ...
{
	r_cm.Null();
	v_cm.Null();
	xx = yy = xy = l = 0;
	int n_data = 8;

	for (int x = thisnode->head_cell_idx; x < thisnode->tail_cell_idx; x++)
		for (int y = thisnode->head_cell_idy; y < thisnode->tail_cell_idy; y++)
		{
			for (int i = 0; i < thisnode->cell[x][y].pid.size(); i++)
				if (thisnode->cell[x][y].pid[i] < Nm)
				{
					C2DVector r, v;
					r = particle[thisnode->cell[x][y].pid[i]].r_original;
					v = (r - rm_old[thisnode->cell[x][y].pid[i]]) / (t - t_old);
					r_cm += r;
					v_cm += v;
					l += r.x * v.y - r.y * v.x;
					xx += r.x*r.x;
					yy += r.y*r.y;
					xy += r.x*r.y;
				}
		}
	double buffer[n_data];
//	double buffer_sum[n_data] = {0};
	double buffer_sum[n_data];

	for (int i = 0; i < n_data; i++)
		buffer_sum[i] =  0;

	buffer[0] = r_cm.x;
	buffer[1] = r_cm.y;
	buffer[2] = v_cm.x;
	buffer[3] = v_cm.y;
	buffer[4] = xx;
	buffer[5] = yy;
	buffer[6] = xy;
	buffer[7] = l;

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Reduce(&buffer, &buffer_sum, n_data, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);

	r_cm.x = buffer_sum[0] / Nm;
	r_cm.y = buffer_sum[1] / Nm;
	v_cm.x = buffer_sum[2] / Nm;
	v_cm.y = buffer_sum[3] / Nm;
	xx = (buffer_sum[4] / Nm) - r_cm.x*r_cm.x;
	yy = (buffer_sum[5] / Nm) - r_cm.y*r_cm.y;
	xy = (buffer_sum[6] / Nm) - r_cm.x*r_cm.y;
	l = (buffer_sum[7] / Nm) - r_cm.x*v_cm.y + r_cm.y*v_cm.x;

	Real lambda1, lambda2;
	lambda1 = 0.5*(xx+yy) + 0.5*sqrt((xx-yy)*(xx-yy) + 4*xy*xy);
	lambda2 = 0.5*(xx+yy) - 0.5*sqrt((xx-yy)*(xx-yy) + 4*xy*xy);

	Rg = sqrt(xx+yy);
	Delta = (lambda1 - lambda2) / (lambda1 + lambda2);
	Delta = Delta*Delta;
	omega = l / (xx+yy);
	Save_Membrane_Position();
}

void Box::Save_All_Variables(std::ostream& os) // Save center of mass position and speed, I, ...
{
	Compute_All_Variables();
	static bool first_time = true;
	if (thisnode->node_id == 0)
	{
		if (first_time)
			os << "#\ttime\tx_cm\ty_cm\tvx_cm\tvy_cm\tangular momentum\tangular freq.\tGyration raduis\tAsphericity\tQxx\tQxy\tQyy" << endl;
		os << t << "\t" << r_cm << "\t" << v_cm << "\t" << l << "\t" << omega << "\t" << Rg << "\t" << Delta << "\t" << xx << "\t" << xy << "\t" << yy << endl;
		first_time = false;
	}
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

		for (int i = 0; i < box->Nm; i++)
		{
			box->particle[i].r_original.write(os);
		}

		for (int i = box->Nm; i < box->N; i++)
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
