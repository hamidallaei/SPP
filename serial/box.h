#ifndef _BOX_
#define _BOX_

#include "../shared/parameters.h"
#include "../shared/c2dvector.h"
#include "../shared/particle.h"
#include "../shared/cell.h"
#include "../shared/wall.h"
#include "../shared/set-up.h"
#include "../shared/state-hyper-vector.h"
#include "../shared/geometry.h"

#include <boost/algorithm/string.hpp>

class Box{
public:
	int N, Ns, Nm, wall_num; // N is the number of particles and wallnum is the number of walls in the system.
	Particle* particle; // Array of particles that we are going to simulate.
	Wall wall[8]; // Array of walls in our system.

	Real Lbx, Lby;
	Real t;
	Real density, packing_fraction;
	Real polarization;
	stringstream info; // information stream that contains the simulation information, like noise, density and etc. this will be used for the saving name of the system.

	Cell cell[divisor_x][divisor_y];


	Box();
	Box(const Real input_Lx, const Real input_Ly, const Real input_phi);
	// I believe that it is better to move these init functions to main files
	void Init_Topology(); // Initialize the wall positions and numbers.
	void Init(Node* input_node, Real input_density); // Intialize the box, positioning particles, giving them velocities, updating cells and sending information to all nodes.
	bool Positioning_Particles(const string name); // Intialize the box from a file, this includes reading particles information, updating cells and sending information to all nodes.

	void Load(const State_Hyper_Vector&); // Load new position and angles of particles and a gsl random generator from a state hyper vector
	void Save(State_Hyper_Vector&) const; // Save current position and angles of particles and a gsl random generator to a state hyper vector

	void Update_Neighbor_List(); // This will update verlet neighore list of each particle
	void Update_Cells();
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
	void Compute_Polarization();
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
	Cell::particle = particle;
}

Box::Box(const Real input_Lx, const Real input_Ly, const Real input_phi)
{
	Ns = 0;
	Nm = 0;
	density = 0;
	wall_num = 0;

	Lbx = input_Lx;
	Lby = input_Ly;

	Lx = input_Lx;
	Ly = input_Ly;
	Lx2 = 2*input_Lx;
	Ly2 = 2*input_Ly;

	Real input_rho = 4*input_phi / (M_PI*Particle::sigma_p*Particle::sigma_p);

	max_N = (int) floor(4*1.1*input_Lx*input_Ly*input_rho);
	particle = new Particle[max_N];
	Cell::particle = particle;
}

// Initialize the wall positions and numbers.
void Box::Init_Topology()
{
	#ifndef PERIODIC_BOUNDARY_CONDITION
		wall_num = 4;
		wall[0].Init(-Lx,-Ly,-Lx, Ly);
		wall[1].Init(-Lx, Ly, Lx, Ly);
		wall[2].Init( Lx, Ly, Lx,-Ly);
		wall[3].Init( Lx,-Ly,-Lx,-Ly);
	#endif
}

void Box::Update_Cells()
{
	for (int x = 0; x < divisor_x; x++)
		for (int y = 0; y < divisor_y; y++)
			cell[x][y].Delete();

	for (int i = 0; i < N; i++)
	{
		int x,y;
		x = (int) (particle[i].r.x + Lx)*divisor_x / Lx2;
		y = (int) (particle[i].r.y + Ly)*divisor_y / Ly2;

		#ifdef DEBUG
		if ((x >= divisor_x) || (x < 0) || (y >= divisor_y) || (y < 0))
		{
			cout << "\n Particle number " << i << " is Out of the box" << endl << flush;
			cout << "Particle Position is " << particle[i].r << endl;
			exit(0);
		}
		#endif

		cell[x][y].Add(i);
	}

	#ifdef verlet_list
		Update_Neighbor_List();
	#endif
}

// This function will update verlet neighore list of particles
void Box::Update_Neighbor_List()
{
// Each cell must interact with itself and 4 of its 8 neihbors that are right cell, up cell, righ up and right down. Because each intertion compute the torque to both particles we need to use 4 of the 8 directions.
	for (int x = 0; x < divisor_x; x++)
		for (int y = 0; y < divisor_y; y++)
		{
			cell[x][y].Clear_Neighbor_List();
			// Self interaction
			cell[x][y].Neighbor_List();
		}

#ifdef PERIODIC_BOUNDARY_CONDITION
// right, up and up right cells:
// The righmost cells and top cells must be excluded to avoid nieghbor node interactions.
	for (int x = 0; x < divisor_x; x++)
		for (int y = 0; y < divisor_y; y++)
		{
			cell[x][y].Neighbor_List(&cell[(x+1)%divisor_x][y]);
			cell[x][y].Neighbor_List(&cell[x][(y+1)%divisor_y]);
			cell[x][y].Neighbor_List(&cell[(x+1)%divisor_x][(y+1)%divisor_y]);
			cell[x][y].Neighbor_List(&cell[(x+1)%divisor_x][(y-1+divisor_y)%divisor_y]);
		}
#else
// right, up and up right cells:
// The righmost cells and top cells must be excluded to avoid nieghbor node interactions.
	for (int x = 0; x < divisor_x-1; x++)
		for (int y = 0; y < divisor_y-1; y++)
		{
			cell[x][y].Neighbor_List(&cell[x+1][y]);
			cell[x][y].Neighbor_List(&cell[x][y+1]);
			cell[x][y].Neighbor_List(&cell[x+1][y+1]);
		}
	for (int x = 0; x < (divisor_x-1); x++)
		cell[x][divisor_y-1].Neighbor_List(&cell[x+1][divisor_y-1]);
	for (int y = 0; y < (divisor_y-1); y++)
		cell[divisor_x-1][y].Neighbor_List(&cell[divisor_x-1][y+1]);
	// right down cell:
// The righmost cells and buttom cells must be excluded to avoid nieghbor node interactions.
	for (int x = 0; x < (divisor_x-1); x++)
		for (int y = 1; y < divisor_y; y++)
			cell[x][y].Neighbor_List(&cell[x+1][y-1]);
#endif
}

// Intialize the box, positioning particles, giving them velocities, updating cells and sending information to all nodes.
void Box::Init(Real input_density)
{
	#ifdef TRACK_PARTICLE
	track_p = &particle[track];
	#endif

	Lbx = Lx;
	Lby = Ly;

	density = input_density;
	Ns = (int) round(Lx2*Ly2*input_density);
	Nm = 0;
	N = Ns + Nm;

	Init_Topology(); // Adding walls

// Positioning the particles at first time. Note that, the positions can be tunned in the main file as well
	if (thisnode->node_id == 0)
	{
		cout << "number_of_particles = " << Ns << endl; // Printing number of particles.
//		Triangle_Lattice_Formation(particle, Ns, 1);
	}

// Buliding up info stream. In next versions we will take this part out of box, making our libraries more abstract for any simulation of SPP.
	info.str("");
}


// Intialize the box from a file, this includes reading particles information, updating cells and sending information to all nodes. 
bool Box::Positioning_Particles(const string input_name)
{
	#ifdef TRACK_PARTICLE
		track_p = &particle[track];
	#endif

	string name = input_name;
	stringstream address(name);

	stringstream command("");
	command << "./fix-file.out ";
	command << name;
	system(command.str().c_str());

	ifstream is;
	is.open(address.str().c_str(), fstream::in);
	if (!is.is_open())
		return false;

	is >> this;

	if (Ns < 0 || Ns > 1000000)
		return (false);

	is.seekg(0,ios_base::end);
	int end_of_is = is.tellg();
	is.seekg(0,ios_base::beg);
	while (is.tellg() < end_of_is && is.tellg() >= 0)
	{
		static int counter = 0;
		is >> this;
		counter++;
	}

	for (int i = 0; i < Nm; i++)
	{
		particle[i].r = particle[i].r_original;
		particle[i].r.Periodic_Transform();
	}

	for (int i = Nm; i < Nm+Ns; i++)
	{
		particle[i].r = particle[i].r_original;
		particle[i].r.Periodic_Transform();
		particle[i].v.x = cos(particle[i].theta);
		particle[i].v.y = sin(particle[i].theta);
		particle[i].Reset();
	}

//	cout << "number_of_particles = " << Ns << endl; // Printing number of particles.
	is.close();

//	density = Ns / (Lx2*Ly2);

	Update_Cells();

	#ifdef verlet_list
		Update_Neighbor_List();
	#endif

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
	for (int i = 0; i < Ns; i++)
	{
		particle[i].r = sv.particle[i].r;
		particle[i].theta = sv.particle[i].theta;
		particle[i].v.x = cos(particle[i].theta);
		particle[i].v.y = sin(particle[i].theta);
	}
	sv.Set_C2DVector_Rand_Generator();

	Update_Cells();

	#ifdef verlet_list
		Update_Neighbor_List();
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

	for (int i = 0; i < N; i++)
	{
		sv.particle[i].r = particle[i].r;
		sv.particle[i].theta = particle[i].theta;
	}
	sv.Get_C2DVector_Rand_Generator();
// We need to make sure that indexing of particles are the same to exactly recompute the same values. Therefor at a saving we update cells and neighore list therefore if we load the same sv and update cells and neighore list we will come to the same indexing
	Update_Cells();

	#ifdef verlet_list
		Update_Neighbor_List();
	#endif
}

// Here the intractio of particles are computed that is the applied tourque to each particle.
void Box::Interact()
{
	#ifdef verlet_list
		for (int x = 0; x < divisor_x; x++)
			for (int y = 0; y < divisor_y; y++)
				cell[x][y].Interact();
	#else
		#ifdef PERIODIC_BOUNDARY_CONDITION
			for (int x = 0; x < divisor_x; x++)
				for (int y = 0; y < divisor_y; y++)
				{
					cell[x][y].Self_Interact();
					cell[x][y].Interact(&cell[(x+1)%divisor_x][y]);
					cell[x][y].Interact(&cell[x][(y+1)%divisor_x]);
					cell[x][y].Interact(&cell[(x+1)%divisor_x][(y+1)%divisor_y]);
					cell[x][y].Interact(&cell[(x+1)%divisor_x][(y-1+divisor_y)%divisor_y]);
				}
		#else
			for (int x = 0; x < divisor_x-1; x++)
				for (int y = 1; y < divisor_y-1; y++)
				{
					cell[x][y].Self_Interact();
					cell[x][y].Interact(&cell[(x+1)%divisor_y][y]);
					cell[x][y].Interact(&cell[x][(y+1)%divisor_y]);
					cell[x][y].Interact(&cell[(x+1)%divisor_y][(y+1)%divisor_y]);
					cell[x][y].Interact(&cell[(x+1)%divisor_y][(y-1+divisor_y)%divisor_y]);
				}
			for (int y = 1; y < divisor_y-1; y++)
			{
				cell[divisor_x-1][y].Self_Interact();
				cell[divisor_x-1][y].Interact(&cell[divisor_x-1][y+1]);
			}
			for (int x = 0; x < divisor_x-1; x++)
			{
				cell[x][0].Self_Interact();
				cell[x][0].Interact(&cell[x+1][0]);
				cell[x][0].Interact(&cell[x][1]);
				cell[x][0].Interact(&cell[x+1][1]);

				cell[x][divisor_x-1].Self_Interact();
				cell[x][divisor_x-1].Interact(&cell[x+1][divisor_y-1]);
				cell[x][divisor_x-1].Interact(&cell[x+1][divisor_y-2]);
			}

			cell[divisor_x-1][divisor_y-1].Self_Interact();
			cell[divisor_x-1][0].Self_Interact();
			cell[divisor_x-1][0].Interact(&cell[divisor_y-1][1]);
		#endif
	#endif
}

// Move all particles of this node.
void Box::Move()
{
	for (int i = 0; i < N; i++)
		particle[i].Move();
}

#ifdef RUNGE_KUTTA2
// Do the first step of second order Runge Kutta
void Box::Move_Runge_Kutta2_1()
{
	for (int i = 0; i < N; i++)
		particle[i].Move_Runge_Kutta2_1();
}

// Do the second step of second order Runge Kutta
void Box::Move_Runge_Kutta2_2()
{
	for (int i = 0; i < N; i++)
		particle[i].Move_Runge_Kutta2_2();
}
#endif

#ifdef RUNGE_KUTTA4
// Do the first step of forth order Runge Kutta
void Box::Move_Runge_Kutta4_1()
{
	for (int i = 0; i < N; i++)
		particle[i].Move_Runge_Kutta4_1();
}

// Do the second step of forth Runge Kutta
void Box::Move_Runge_Kutta4_2()
{
	for (int i = 0; i < N; i++)
		particle[i].Move_Runge_Kutta4_2();
}

// Do the third step of forth Runge Kutta
void Box::Move_Runge_Kutta4_3()
{
	for (int i = 0; i < N; i++)
		particle[i].Move_Runge_Kutta4_3();
}

// Do the forth step of forth Runge Kutta
void Box::Move_Runge_Kutta4_4()
{
	for (int i = 0; i < N; i++)
		particle[i].Move_Runge_Kutta4_4();
}
#endif

// One full step, composed of interaction computation and move.
void Box::One_Step()
{
	#ifdef RUNGE_KUTTA2
		Interact();
		Move_Runge_Kutta2_1();

		Interact();
		Move_Runge_Kutta2_2();
	#else
		#ifdef RUNGE_KUTTA4
			Interact();
			Move_Runge_Kutta4_1();

			Interact();
			Move_Runge_Kutta4_2();

			Interact();
			Move_Runge_Kutta4_3();

			Interact();
			Move_Runge_Kutta4_4();
		#else
			Interact();
			Move();
		#endif
	#endif

	t += dt;d
}

// Several steps befor a cell upgrade.
void Box::Multi_Step(int steps)
{
	for (int i = 0; i < steps; i++)
	{
		Interact();
		Move();
	}
	Update_Cells();
}

// Several steps with a cell upgrade call after each interval.
void Box::Multi_Step(int steps)
{
	for (int i = 0; i < steps; i++)
	{
		#ifdef RUNGE_KUTTA2
			Interact();
			Move_Runge_Kutta2_1();

			Interact();
			Move_Runge_Kutta2_2();
		#else
			#ifdef RUNGE_KUTTA4
				Interact();
				Move_Runge_Kutta4_1();

				Interact();
				Move_Runge_Kutta4_2();

				Interact();
				Move_Runge_Kutta4_3();

				Interact();
				Move_Runge_Kutta4_4();
			#else
				Interact();
				Move();
			#endif
		#endif
	}

	t += dt*steps;
//	cout << "Updating cells" << endl;
	Update_Cells();
//	cout << "Updating finished" << endl;
	#ifdef verlet_list
		Update_Neighbor_List();
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
	for (int i = 0; i < N; i++)
	{
		particle[i].r += d;
		particle[i].r.Periodic_Transform();
	}
	Update_Cells();
	#ifdef verlet_list
		Update_Neighbor_List();
	#endif
}

void Box::Compute_Polarization()
{
	C2DVector vp_sum;
	for (int i = 0; i < Ns; i++)
		vp_sum += particle[i].v;
	vp_sum /= N;
	polarization = sqrt(vp_sum.Square());
}

// Save polarization of the particles inside the box
void Box::Save_Polarization(std::ostream& os)
{
	Compute_Polarization();
	if (thisnode->node_id == 0)
		os << t << "\t" << polarization << endl;
}

// Saving the particle information (position and velocities) to a standard output stream (probably a file). This must be called by only the root.
std::ostream& operator<<(std::ostream& os, Box* box)
{
	// binary output
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

	for (int i = box->Nm; i < (box->Nm+box->Ns); i++)
	{
		box->particle[i].Write(os);
	}
}

// Reading the particle information (position and velocities) from a standard input stream (probably a file).
std::istream& operator>>(std::istream& is, Box* box)
{
	// binary input
	is.read((char*) &(box->t), sizeof(Real) / sizeof(char));
	is.read((char*) &(box->Lbx), sizeof(Real) / sizeof(char));
	is.read((char*) &(box->Lby), sizeof(Real) / sizeof(char));
	is.read((char*) &(Particle::nb), sizeof(int) / sizeof(char));
	is.read((char*) &(box->Ns), sizeof(box->Ns) / sizeof(char));
	is.read((char*) &(box->Nm), sizeof(box->Nm) / sizeof(char));

	for (int i = 0; i < box->Nm; i++)
	{
		is >> box->particle[i].r_original;
	}

	for (int i = box->Nm; i < (box->Nm + box->Ns); i++)
	{
		box->particle[i].Read(is);
	}
}

#endif
