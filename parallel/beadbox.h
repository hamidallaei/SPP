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
	Real* theta_old; // old angle of self-propelsion
	C2DVector* r_old; // old position of beads
	C2DVector mb_r_cm, mb_v_cm, sw_r_cm, sw_v_cm;
	Real mb_xx, mb_yy, mb_xy, sw_xx, sw_yy, sw_xy; // Qxx, Qyy, Qxy. Q = 1/N sum (\vec{r} - \vec{r}_cm) (\vec{r} - \vec{r}_cm)
	Real mb_l, sw_l; // Angular momentum per particle
	Real mb_omega, sw_omega, sw_spin; // Angular velocity
	Real mb_lambda1, mb_lambda2, sw_lambda1, sw_lambda2; // Eigenvalues of Q tensor
	Real mb_Rg, sw_Rg; // Gyration radius
	Real mb_Delta, sw_Delta; // Asphericity

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
	void Save_Particles_Positions(); // Save the position of membrane beads to compute their velocities in future.
	void Compute_All_Variables(); // Compute center of mass position and speed, I, ...
	void Save_All_Variables(std::ostream& os); // Save center of mass position and speed, I, ...

	friend std::ostream& operator<<(std::ostream& os, Box* box); // Save
	friend std::istream& operator>>(std::istream& is, Box* box); // Input
};

Box::Box()
{
	N = 0;
	particle = new Particle[max_N];
	r_old = new C2DVector[max_N];
	theta_old = new Real[max_N];
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

#ifdef RUNGE_KUTTA
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
#endif

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

void Box::Save_Particles_Positions()
{
	t_old = t;
	for (int i = 0; i < N; i++)
	{
		r_old[i] = particle[i].r_original;
		theta_old[i] = particle[i].theta;
	}
}

void Box::Compute_All_Variables() // Compute center of mass position and speed, I, ...
{
	mb_r_cm.Null();
	mb_v_cm.Null();
	sw_r_cm.Null();
	sw_v_cm.Null();
	mb_xx = mb_yy = mb_xy = mb_l = 0;
	sw_xx = sw_yy = sw_xy = sw_l = 0;
	int n_data = 17;

	for (int x = thisnode->head_cell_idx; x < thisnode->tail_cell_idx; x++)
		for (int y = thisnode->head_cell_idy; y < thisnode->tail_cell_idy; y++)
		{
			for (int i = 0; i < thisnode->cell[x][y].pid.size(); i++)
				if (thisnode->cell[x][y].pid[i] < Nm)
				{
					C2DVector r, v;
					r = particle[thisnode->cell[x][y].pid[i]].r_original;
					v = (r - r_old[thisnode->cell[x][y].pid[i]]) / (t - t_old);
					mb_r_cm += r;
					mb_v_cm += v;
					mb_l += r.x * v.y - r.y * v.x;
					mb_xx += r.x*r.x;
					mb_yy += r.y*r.y;
					mb_xy += r.x*r.y;
				}
				else
				{
					C2DVector r, v;
					Real spin, dtheta;
					r = particle[thisnode->cell[x][y].pid[i]].r_original;
					v = (r - r_old[thisnode->cell[x][y].pid[i]]) / (t - t_old);
					dtheta = particle[thisnode->cell[x][y].pid[i]].theta - theta_old[thisnode->cell[x][y].pid[i]];
					dtheta -= 2*M_PI*((int) floor(dtheta / (2*M_PI)));
					dtheta -= 2*M_PI*((int) floor(dtheta / (M_PI)));
					spin = (dtheta) / (t - t_old);
					sw_r_cm += r;
					sw_v_cm += v;
					sw_spin += spin;
					sw_l += r.x * v.y - r.y * v.x;
					sw_xx += r.x*r.x;
					sw_yy += r.y*r.y;
					sw_xy += r.x*r.y;
				}
		}
	double buffer[n_data];
//	double buffer_sum[n_data] = {0};
	double buffer_sum[n_data];

	for (int i = 0; i < n_data; i++)
		buffer_sum[i] =  0;

	buffer[0] = mb_r_cm.x;
	buffer[1] = mb_r_cm.y;
	buffer[2] = mb_v_cm.x;
	buffer[3] = mb_v_cm.y;
	buffer[4] = mb_xx;
	buffer[5] = mb_yy;
	buffer[6] = mb_xy;
	buffer[7] = mb_l;

	buffer[8] = sw_r_cm.x;
	buffer[9] = sw_r_cm.y;
	buffer[10] = sw_v_cm.x;
	buffer[11] = sw_v_cm.y;
	buffer[12] = sw_xx;
	buffer[13] = sw_yy;
	buffer[14] = sw_xy;
	buffer[15] = sw_l;
	buffer[16] = sw_spin;

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Reduce(&buffer, &buffer_sum, n_data, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);

	mb_r_cm.x = buffer_sum[0] / Nm;
	mb_r_cm.y = buffer_sum[1] / Nm;
	mb_v_cm.x = buffer_sum[2] / Nm;
	mb_v_cm.y = buffer_sum[3] / Nm;
	mb_xx = (buffer_sum[4] / Nm) - mb_r_cm.x*mb_r_cm.x;
	mb_yy = (buffer_sum[5] / Nm) - mb_r_cm.y*mb_r_cm.y;
	mb_xy = (buffer_sum[6] / Nm) - mb_r_cm.x*mb_r_cm.y;
	mb_l = (buffer_sum[7] / Nm) - mb_r_cm.x*mb_v_cm.y + mb_r_cm.y*mb_v_cm.x;

	sw_r_cm.x = buffer_sum[8] / Ns;
	sw_r_cm.y = buffer_sum[9] / Ns;
	sw_v_cm.x = buffer_sum[10] / Ns;
	sw_v_cm.y = buffer_sum[11] / Ns;
	sw_xx = (buffer_sum[12] / Ns) - sw_r_cm.x*sw_r_cm.x;
	sw_yy = (buffer_sum[13] / Ns) - sw_r_cm.y*sw_r_cm.y;
	sw_xy = (buffer_sum[14] / Ns) - sw_r_cm.x*sw_r_cm.y;

	sw_l = (buffer_sum[15] / Ns) - sw_r_cm.x*mb_v_cm.y + sw_r_cm.y*mb_v_cm.x - mb_r_cm.x*sw_v_cm.y + mb_r_cm.y*sw_v_cm.x + mb_r_cm.x*mb_v_cm.y - mb_r_cm.y*mb_v_cm.x;
	sw_spin = (buffer_sum[16] / Ns);

	mb_lambda1 = 0.5*(mb_xx + mb_yy) + 0.5*sqrt((mb_xx - mb_yy)*(mb_xx - mb_yy) + 4*mb_xy*mb_xy);
	mb_lambda2 = 0.5*(mb_xx + mb_yy) - 0.5*sqrt((mb_xx - mb_yy)*(mb_xx - mb_yy) + 4*mb_xy*mb_xy);

	mb_Rg = sqrt(mb_xx + mb_yy);
	mb_Delta = (mb_lambda1 - mb_lambda2) / (mb_lambda1 + mb_lambda2);
	mb_Delta = mb_Delta*mb_Delta;
	mb_omega = mb_l / (mb_xx + mb_yy);

	sw_lambda1 = 0.5*(sw_xx + sw_yy) + 0.5*sqrt((sw_xx - sw_yy)*(sw_xx - sw_yy) + 4*sw_xy*sw_xy);
	sw_lambda2 = 0.5*(sw_xx + sw_yy) - 0.5*sqrt((sw_xx - sw_yy)*(sw_xx - sw_yy) + 4*sw_xy*sw_xy);

	sw_Rg = sqrt(sw_xx + sw_yy);
	sw_Delta = (sw_lambda1 - sw_lambda2) / (sw_lambda1 + sw_lambda2);
	sw_Delta = sw_Delta * sw_Delta;
	sw_omega = sw_l / (sw_xx + sw_yy + (mb_r_cm - sw_r_cm).Square());

	Save_Particles_Positions();
}

void Box::Save_All_Variables(std::ostream& os) // Save center of mass position and speed, I, ...
{
	Compute_All_Variables();
	static bool first_time = true;
	if (thisnode->node_id == 0)
	{
		if (first_time)
			os << "#\ttime\tx_cm\ty_cm\tvx_cm\tvy_cm\tangular momentum\tangular freq.\tGyration raduis\tAsphericity\tQxx\tQxy\tQyy\tsw. ang. freq." << endl;
		os << t << "\t" << mb_r_cm << "\t" << mb_v_cm << "\t" << mb_l << "\t" << mb_omega << "\t" << mb_Rg << "\t" << mb_Delta << "\t" << mb_xx << "\t" << mb_xy << "\t" << mb_yy << "\t" << sw_omega << endl;
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
