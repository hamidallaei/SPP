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
	int N, wall_num; // N is the number of particles and wallnum is the number of walls in the system.
	Particle particle[max_N]; // Array of particles that we are going to simulate.
	Wall wall[8]; // Array of walls in our system.

	Real density;
	stringstream info; // information stream that contains the simulation information, like noise, density and etc. this will be used for the saving name of the system.

	Node* thisnode; // Node is a class that has information about the node_id and its boundaries, neighbores and etc.

	Box();
	void Init_Topology(); // Initialize the wall positions and numbers.
	void Init(Node* input_node, Real input_density); // Intialize the box, positioning particles, giving them velocities, updating cells and sending information to all nodes.
	bool Init(Node* input_node, const string name); // Intialize the box from a file, this includes reading particles information, updating cells and sending information to all nodes.

	void Load(const State_Hyper_Vector&); // Load new position and angles of particles and a gsl random generator from a state hyper vector
	void Save(State_Hyper_Vector&) const; // Save current position and angles of particles and a gsl random generator to a state hyper vector
	void Add_Deviation(const State_Hyper_Vector&); // Add argument state vector as a deviation to the state of the box

	void Interact(); // Here the intractio of particles are computed that is the applied tourque to each particle.
	void Move(); // Move all particles of this node.
	void One_Step(); // One full step, composed of interaction computation and move.
	void Multi_Step(int steps); // Several steps befor a cell upgrade.
	void Translate(C2DVector d); // Translate position of all particles with vector d

	Real Rlative_Change(int tau);
	Real Mean_Rlative_Change(int tau, int number_of_tries, Real& error);
	void Lyapunov_Exponent(); // Finding the largest lyapunov exponent

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
	N = (int) round(Lx2*Ly2*input_density);

	Init_Topology(); // Adding walls

	if (thisnode->node_id == 0)
	{
		cout << "number_of_particles = " << N << endl; // Printing number of particles.
// Positioning the particles
//		Polar_Formation(particle,N);
//		Triangle_Lattice_Formation(particle, N, 1);
		Random_Formation(particle, N, 0); // Positioning partilces Randomly, but distant from walls (the last argument is the distance from walls)
//		Random_Formation_Circle(particle, N, Lx-1); // Positioning partilces Randomly, but distant from walls
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


// Intialize the box from a file, this includes reading particles information, updating cells and sending information to all nodes. 
bool Box::Init(Node* input_node, const string input_name)
{
	#ifdef TRACK_PARTICLE
	track_p = &particle[track];
	#endif

	Real input_mu_plus;
	Real input_mu_minus;
	Real input_Dphi;
	Real input_L;

	string name = input_name;
	stringstream address(name);
	ifstream is;
	is.open(address.str().c_str());
	if (!is.is_open())
		return false;

	boost::replace_all(name, "-r-v.bin", "");
	boost::replace_all(name, "rho=", "");
	boost::replace_all(name, "-mu+=", "\t");
	boost::replace_all(name, "-mu-=", "\t");
	boost::replace_all(name, "-Dphi=", "\t");
	boost::replace_all(name, "-L=", "\t");

	stringstream ss_name(name);
	ss_name >> density;
	ss_name >> input_mu_plus;
	ss_name >> input_mu_minus;
	ss_name >> input_Dphi;
	ss_name >> input_L;
	if (input_L != Lx_int)
	{
		cout << "The specified box size " << input_L << " is not the same as the size in binary file which is " << Lx_int << " please recompile the code with the right Lx_int in parameters.h file." << endl;
		return false;
	}

	Particle::mu_plus = input_mu_plus;
	Particle::mu_minus = input_mu_minus;
	Particle::D_phi = input_Dphi;

	thisnode = input_node;

	is.read((char*) &N, sizeof(int) / sizeof(char));
	if (N < 0 || N > 1000000)
		return (false);

	for (int i = 0; i < N; i++)
	{
		is >> particle[i].r;
		is >> particle[i].v;
	}

	is.close();

	MPI_Barrier(MPI_COMM_WORLD);

	Init_Topology(); // Adding walls

	if (thisnode->node_id == 0)
	{
		cout << "number_of_particles = " << N << endl; // Printing number of particles.
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
	return (true);
}


// Loading a state to the box.
void Box::Load(const State_Hyper_Vector& sv)
{
	if (N != sv.N)
	{
		cout << "Error: Number of particles in state vectors differ from box" << endl;
		exit(0);
	}
	for (int i = 0; i < N; i++)
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
	if (N != sv.N)
	{
		cout << "Error: Number of particles in state vectors differ from box" << endl;
		exit(0);
	}
	thisnode->Root_Gather();
	thisnode->Root_Bcast();
	for (int i = 0; i < N; i++)
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

// Add argument state vector as a deviation to the state of the box
void Box::Add_Deviation(const State_Hyper_Vector& dsv)
{
	thisnode->Root_Gather();
	if (N != dsv.N)
	{
		cout << "Error: Number of particles in state vectors differ from box" << endl;
		exit(0);
	}
	for (int i = 0; i < N; i++)
	{
		particle[i].r += dsv.particle[i].r;
		particle[i].r.Periodic_Transform();
		particle[i].theta += dsv.particle[i].theta;
		particle[i].v.x = cos(particle[i].theta);
		particle[i].v.y = sin(particle[i].theta);
	}
	thisnode->Root_Bcast();
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

// Finding the relative change of perturbation
Real Box::Rlative_Change(int tau)
{
	static State_Hyper_Vector gamma(N);
	static State_Hyper_Vector gamma_prime(N);
	static State_Hyper_Vector gamma_0(N);
	static State_Hyper_Vector dgamma_0(N);
	static State_Hyper_Vector dgamma(N);
	
	Save(gamma_0);
	MPI_Barrier(MPI_COMM_WORLD);


	for (int i = 0; i < tau/20; i++)
		Multi_Step(20);
	Multi_Step(tau % 20);
	
	Save(gamma);
	Load(gamma_0);
	dgamma_0.Rand(0.0,1e-3);
	MPI_Barrier(MPI_COMM_WORLD);
	Add_Deviation(dgamma_0);

	for (int i = 0; i < tau/20; i++)
		Multi_Step(20);
	Multi_Step(tau % 20);

	Save(gamma_prime);
	Load(gamma);

	dgamma = gamma_prime - gamma;

	Real d1 = dgamma.Square();
	Real d0 = dgamma_0.Square();

	if (thisnode->node_id == 0)
	{
		int index = dgamma.Max_Index();
//		cout << tau << "\t" << d1 << endl;
//		cout << tau << "\t" << d0 << endl;
		cout << tau << "\t" << sqrt(d1/d0) << endl;
		cout << tau << "\t" << log(sqrt(d1/d0)) / (tau*dt) << endl;
		Real dtheta = dgamma.particle[index].theta;
		dtheta -= 2*M_PI*ceil((dtheta - M_PI) / (2*M_PI));
		cout << index << "\t" << dtheta << endl;
	}
	
	return (sqrt(d1/d0));
}

// Finding the mean relative change of perturbation
Real Box::Mean_Rlative_Change(int tau, int number_of_tries, Real& error)
{
	Real result = 0;
	error = 0;
	for (int i = 0; i < number_of_tries; i++)
	{
		Real rc = Rlative_Change(tau);
		result += rc;
		error += rc*rc;
	}
	result /= number_of_tries;
	error /= number_of_tries;
	error -= result*result;
	error /= sqrt(number_of_tries);
	return (result);
}

// Finding the largest lyapunov exponent
void Box::Lyapunov_Exponent()
{
	stringstream adress;
	adress.str("");
	adress << info.str() << "-deviation.dat";
	ofstream outfile(adress.str().c_str());
	for (Real t = 0.1; t < 4; t+=0.4)
	{
		if (thisnode->node_id == 0)
			cout << "tau is: " << t << endl;
		int i = (int) round(t/dt);
		Real error;
		Real mrc = Mean_Rlative_Change(i, 100, error);
		if (thisnode->node_id == 0)
			outfile << i*dt << "\t" << mrc << "\t" << error << endl;
	}
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
