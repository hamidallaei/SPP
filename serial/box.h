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
	int N; // N is the number of particles and wallnum is total number of walls in the system.
	Particle particle[max_N]; // Array of particles that we are going to simulate.
	Geometry geometry; // the entire geometry that the particle interact with. 
	Cell cell[divisor_x][divisor_y];

	Real density;
	stringstream info; // information stream that contains the simulation information, like noise, density and etc. this will be used for the saving name of the system.

	Box();

	void Load(const State_Hyper_Vector&); // Load new position and angles of particles and a gsl random generator from a state hyper vector
	void Save(State_Hyper_Vector&) const; // Save current position and angles of particles and a gsl random generator to a state hyper vector

	void Update_Neighbor_List(); // This will update verlet neighore list of each particle
	void Update_Cells();
	void Interact(); // Here the intractio of particles are computed that is the applied tourque to each particle.
	void Move(); // Move all particles of this node.
	void One_Step(); // One full step, composed of interaction computation and move.
	void Multi_Step(int steps); // Several steps befor a cell upgrade.
	void Multi_Step(int steps, int interval); // Several steps with a cell upgrade call after each interval.
	void Translate(C2DVector d); // Translate position of all particles with vector d
	void Center();

	void Make_Traj(Real scale, std::ofstream& data_file);

	friend std::ostream& operator<<(std::ostream& os, Box* box); // Save
	friend std::istream& operator>>(std::istream& is, Box* box); // Input
};

Box::Box()
{
	N = 0;
	density = 0;
	#ifdef TRACK_PARTICLE
		track_p = &particle[track];
	#endif

	Cell::particle = particle;

	for (int i = 0; i < divisor_x; i++)
		for (int j = 0; j < divisor_y; j++)
			cell[i][j].Init((Real) Lx*(2*i-divisor_x + 0.5)/divisor_x, (Real) Ly*(2*j-divisor_y + 0.5)/divisor_y);
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

	Update_Cells();
}

// Saving state of the box.
void Box::Save(State_Hyper_Vector& sv) const
{
	if (N != sv.N)
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

	for(int i = 0 ; i < N; i++)
		geometry.Interact(&particle[i]);
}

// Move all particles of this node.
void Box::Move()
{
	for (int i = 0; i < N; i++)
		particle[i].Move();
}

// One full step, composed of interaction computation and move.
void Box::One_Step()
{
	Interact();
	Move();
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
void Box::Multi_Step(int steps, int interval)
{
	for (int i = 0; i < steps/interval; i++)
		Multi_Step(interval);
	Multi_Step(steps % interval);
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
}

void Box::Center()
{
	C2DVector cm;
	for (int i = 0; i < N; i++)
		cm -= particle[i].r;
	cm = cm / N;
	Translate(cm);
}



void Box::Make_Traj(Real scale, ofstream& data_file)
{
	data_file << N << endl;
	data_file << "something" << endl;
	for (int i = 0; i < N; i++)
		data_file << "H	" << particle[i].r * scale << "\t" << 0.0 << endl;
}

// Saving the particle information (position and velocities) to a standard output stream (probably a file). This must be called by only the root.
std::ostream& operator<<(std::ostream& os, Box* box)
{
	// binary output
	os.write((char*) &box->N, sizeof(box->N) / sizeof(char));
	for (int i = 0; i < box->N; i++)
	{
		box->particle[i].r.write(os);
		C2DVector v;
		v.x = cos(box->particle[i].theta);
		v.y = sin(box->particle[i].theta);
		v.write(os);
	}

//	// txt output
//	os.write((char*) &box->N, sizeof(box->N) / sizeof(char));
//	for (int i = 0; i < box->N; i++)
//	{
//		box->particle[i].r.write(os);
//		C2DVector v;
//		v.x = cos(box->particle[i].theta);
//		v.y = sin(box->particle[i].theta);
//		v.write(os);
//	}
}

// Reading the particle information (position and velocities) from a standard input stream (probably a file).
std::istream& operator>>(std::istream& is, Box* box)
{
	// binary input
	is.read((char*) &box->N, sizeof(int) / sizeof(char));
	for (int i = 0; i < box->N; i++)
	{
		is >> box->particle[i].r;
		is >> box->particle[i].v;
	}

//	// txt input
//	is.read((char*) &box->N, sizeof(int) / sizeof(char));
//	for (int i = 0; i < box->N; i++)
//	{
//		is >> box->particle[i].r;
//		is >> box->particle[i].v;
//	}
}

#endif
