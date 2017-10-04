#include <boost/algorithm/string.hpp>

#include<iostream>
#include<cstdlib>
#include<vector>

#include"analyze.h"

using namespace std;

const double threshold_r2 = 0.1;
SceneSet* sceneset;

VisualParticle* particle;

bool Cell::Check(Cell* c)
{
	for (int i = 0; i < pid.size(); i++)
		for (int j = 0; j < c->pid.size(); j++)
		{
			double d2 = (particle[pid[i]].r - particle[c->pid[j]].r).Square();
			if (d2 < threshold_r2)
				return true;
		}
	return false;
}

bool Cell::Self_Check()
{
	for (int i = 0; i < pid.size(); i++)
		for (int j = i+1; j < pid.size(); j++)
		{
			double d2 = (particle[pid[i]].r - particle[pid[j]].r).Square();
			if (d2 < threshold_r2)
				return true;
		}
	return false;
}

Cell** cell;

void Make_Cell_Array()
{
	cell = new Cell*[divisor_x];
	for (int x = 0; x < divisor_x; x++)
		cell[x] = new Cell[divisor_y];
}

void Clear_Cell_Array()
{
	for (int x = 0; x < divisor_x; x++)
		for (int y = 0; y < divisor_y; y++)
			cell[x][y].Delete();
}

void Delete_Cell_Array()
{
	Clear_Cell_Array();
	for (int x = 0; x < divisor_x; x++)
		delete[] cell[x];
	delete[] cell;
}

bool Check_Overlap_Oldstyle(Scene& sc)
{
	for (int j = 0; j < sc.Ns; j++)
		sc.sparticle[j].r.Periodic_Transform();
	for (int j = 0; j < sc.Ns; j++)
		for (int k = j+1; k < sc.Ns; k++)
		{
			double d2 = (sc.sparticle[j].r - sc.sparticle[k].r).Square();
			if (d2 < threshold_r2)
				return true;
		}
	return false;
}


bool Check_Overlap_Newstyle(Scene& sc)
{
	particle = sc.sparticle;
	for (int k = 0; k < sc.Ns; k++)
	{
		int x,y;
		particle[k].r.Periodic_Transform();
		x = (int) (particle[k].r.x + sc.L.x)*divisor_x / (sc.L.x*2);
		y = (int) (particle[k].r.y + sc.L.y)*divisor_y / (sc.L.y*2);
		cell[x][y].Add(k);
	}
	bool b = false;
	for (int x = 0; x < divisor_x; x++)
		for (int y = 0; y < divisor_y; y++)
		{
			b = b || cell[x][y].Self_Check();
			b = b || cell[x][y].Check(&cell[(x+1)%divisor_x][y]);
			b = b || cell[x][y].Check(&cell[(x-1+divisor_x)%divisor_x][y]);

			b = b || cell[x][y].Check(&cell[x][(y+1)%divisor_y]);
			b = b || cell[x][y].Check(&cell[x][(y-1+divisor_y)%divisor_y]);

			b = b || cell[x][y].Check(&cell[(x+1)%divisor_x][(y+1)%divisor_y]);
			b = b || cell[x][y].Check(&cell[(x+1)%divisor_x][(y-1+divisor_y)%divisor_y]);
			b = b || cell[x][y].Check(&cell[(x-1+divisor_x)%divisor_x][(y+1)%divisor_y]);
			b = b || cell[x][y].Check(&cell[(x-1+divisor_x)%divisor_x][(y-1+divisor_y)%divisor_y]);

			if (b)
				return b;
		}
	return b;
}

int main(int argc, char** argv)
{
	SavingVector::Init_Rand(321);

	for (int i = 1; i < argc; i++)
	{
		string name = argv[i];
		sceneset = new SceneSet(name);
		cout << "Reading the file" << endl;
		bool read_state = sceneset->Read();
		cout << "Finished" << endl;
		Lx = sceneset->scene[2].L.x;
		Ly = sceneset->scene[2].L.y;
		Lx2 = 2*Lx;
		Ly2 = 2*Ly;
		divisor_x = (int) Lx2;
		divisor_y = (int) Ly2;
		Make_Cell_Array();
		cout << "Check the file" << endl;
		bool overlap = false;
		if (read_state)
		{
			for (int frame = 0; frame < sceneset->Nf; frame++)
			{
				overlap = Check_Overlap_Newstyle(sceneset->scene[frame]);
				if (overlap)
				{
					cout << "The file is curropted at frame " << frame << " of total " << sceneset->Nf << " frames" << endl;
					cout << "The file is curropted at time " << sceneset->scene[frame].t << " of total " << sceneset->scene[sceneset->Nf-1].t << " units of time" << endl;
					exit(0);
				}
			}
		}
		else
			cout << "Can not read the file" << endl;
		cout << "Finished" << endl;

		if (overlap)
			cout << "The file is curropted." << endl;
		else
			cout << "The file is healthy." << endl;

		delete sceneset;
		Delete_Cell_Array();
	}

	return 0;
}

