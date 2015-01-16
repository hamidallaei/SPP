#ifndef _BOX_
#define _BOX_

#include "../shared/parameters.h"
#include "../shared/c2dvector.h"
#include "../shared/particle.h"
#include "../shared/cell.h"
#include "../shared/wall.h"
#include "../shared/set-up.h"



class Box{
public:
	int N, wall_num; // N is the number of particles and wallnum is the number of walls in the system.
	Particle particle[max_N]; // Array of particles that we are going to simulate.
	Wall wall[8]; // Array of walls in our system.
	Cell cell[divisor_x][divisor_y];

	Real density, volume_fraction;
	stringstream info;

	Box();
	void Update_Cells();
	void Interact();
	void Move();
	void One_Step();
	void Multi_Step(int steps);
	void Translate(C2DVector d);
	void Center();

	void Make_Traj(Real scale, std::ofstream& data_file);

	friend std::ostream& operator<<(std::ostream& os, Box* box);
	friend std::istream& operator>>(std::istream& is, Box* box);
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

	#pragma omp parallel for
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
}


void Box::Interact()
{
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

	for (int i = 0; i < wall_num; i++)
		for (int j = 0; j < N; j++)
			wall[i].Interact(&particle[j]);
	#endif

//	for (int i = 0; i < N; i++)
//		for (int j = 0; j < 4; j++)
//			particle[i].Interact(&wall[j]);
}


void Box::Move()
{
	for (int i = 0; i < N; i++)
		particle[i].Move();
}

void Box::One_Step()
{
	Interact();
	Move();
}

void Box::Multi_Step(int steps)
{
	for (int i = 0; i < steps; i++)
	{
		Interact();
		Move();
	}
	Update_Cells();
}

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

std::ostream& operator<<(std::ostream& os, Box* box)
{
	os.write((char*) &box->N, sizeof(box->N) / sizeof(char));
	for (int i = 0; i < box->N; i++)
	{
		box->particle[i].r.write(os);
		C2DVector v;
		v.x = cos(box->particle[i].theta);
		v.y = sin(box->particle[i].theta);
		v.write(os);
	}
}

std::istream& operator>>(std::istream& is, Box* box)
{
	is.read((char*) &box->N, sizeof(int) / sizeof(char));
	for (int i = 0; i < box->N; i++)
	{
		is >> box->particle[i].r;
		is >> box->particle[i].v;
	}
}

#endif
