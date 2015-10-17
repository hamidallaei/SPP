#ifndef _FIELD_
#define _FIELD_

#include "../shared/c2dvector.h"
#include "../shared/particle.h"
#include <boost/tuple/tuple.hpp>
#include "gnuplot-iostream.h"
#include <vector>

class Field_Cell{
public:
	vector<BasicParticle*> particle;
	vector<Real> dtheta;
	C2DVector r,v,W;
	Real density;
	Real cohesion; // cohesion shows the amount of velocity cohesion
	Real omega;
	Real curl;
	Real theta_ave;
	static Real dim_x;

	Field_Cell();
	~Field_Cell();
	void Init(Real,Real);
	void Reset();
	void Add(BasicParticle* p);
	void Compute_Fields(Real L);
	void Delta_Theta_Stat();
};

Real Field_Cell::dim_x = 0;

Field_Cell::Field_Cell()
{
}

Field_Cell::~Field_Cell()
{
	Reset();
}

void Field_Cell::Init(Real x, Real y)
{
	r.x = x;
	r.y = y;
	Reset();
}

void Field_Cell::Reset()
{
	particle.clear();
	v.Null();
	W.Null();
	density = 0;
	cohesion = 0;
	curl = 0;
	omega = 0;
}

void Field_Cell::Add(BasicParticle* p)
{
	particle.push_back((BasicParticle*) p);
}

void Field_Cell::Compute_Fields(Real L)
{
	v.Null();
	dtheta.clear();
	theta_ave = cohesion = omega = 0;
	density = particle.size()/(dim_x*dim_x);
	for (int i = 0; i < particle.size(); i++)
	{
		v += particle[i]->v;
		Real theta_i = atan2(particle[i]->v.y,particle[i]->v.x);
		dtheta.push_back(theta_i);
		theta_ave += theta_i;
		omega += (particle[i]->r - r).x*particle[i]->v.y - (particle[i]->r - r).y*particle[i]->v.x;
		for (int j = i+1; j < particle.size(); j++)
			cohesion += (particle[i]->v*particle[j]->v)/sqrt(particle[i]->v.Square()*particle[j]->v.Square());
	}
	v /= particle.size();
	if (particle.size() > 0)
		theta_ave /= particle.size();
	omega /= particle.size();
	cohesion /= particle.size()*(particle.size() - 1);
	cohesion *= 2;
	W = v*(density);
	if (particle.size() < 2)
	{
		cohesion = 0;
		omega = 0;
		v.Null();
		W.Null();
	}
	for (int i = 0; i < particle.size(); i++)
	{
		dtheta[i] -= theta_ave;
	}
}

class Field{
public:
	Field_Cell** cell;
	int grid_dim_x, sample; // sample is number of sampled field that we average. we need it to compute the average by dividing the summation of fields
	Real L;
	Field(int input_grid_dim_x, Real);
	~Field();
	void Init();
	void Compute(BasicParticle* particle, int N);
	void Save(ofstream& data_file);
	void Save_Theta_Deviation(ofstream& data_file);
	void Angle_Deviation(ofstream& data_file);
	void Draw(string);
	void Draw_Section(string, int);
	void Draw_Density_Contour(string, double);
	void Add(Field* f);
	void Average();
	void Reset();
};

Field::Field(int input_grid_dim_x, Real input_L): L(input_L), grid_dim_x(input_grid_dim_x)
{
	Field_Cell::dim_x = 2*L / grid_dim_x;
	cell = new Field_Cell*[grid_dim_x];
	for (int i = 0; i < grid_dim_x; i++)
	{
		cell[i] = new Field_Cell[grid_dim_x];
	}
	Init();
}

Field::~Field()
{
	for (int i = 0; i < grid_dim_x; i++)
		delete [] cell[i];
	delete [] cell;
}

void Field::Init()
{
	for (int x = 0; x < grid_dim_x; x++)
		for (int y = 0; y < grid_dim_x; y++)
		{
			cell[x][y].Reset();
			cell[x][y].Init(2*(x+0.5)*(L / grid_dim_x) - L, 2*(y+0.5)*(L / grid_dim_x) - L);
		}
	Reset();
}

void Field::Compute(BasicParticle* particle, int N)
{
	for (int i = 0; i < N; i++)
	{
		int x = (int) floor((particle[i].r.x + L)*grid_dim_x / (2*L));
		int y = (int) floor((particle[i].r.y + L)*grid_dim_x / (2*L));
		cell[x][y].Add(&particle[i]);
	}

	for (int i = 0; i < grid_dim_x; i++)
		for (int j = 0; j < grid_dim_x; j++)
			cell[i][j].Compute_Fields(L);
	for (int i = 0; i < grid_dim_x; i++)
		for (int j = 0; j < grid_dim_x; j++)
			cell[i][j].curl = (grid_dim_x/L)*(cell[(i+1)%grid_dim_x][j].v.y - cell[i][j].v.y) - (grid_dim_x/L)*(cell[i][(j+1)%grid_dim_x].v.x - cell[i][j].v.x);
}

void Field::Save(ofstream& data_file)
{
	for (int i = 0; i < grid_dim_x; i++)
	{
		for (int j = 0; j < grid_dim_x; j++)
			data_file << cell[i][j].r << "\t" << cell[i][j].v << "\t" << cell[i][j].density << "\t"<< cell[i][j].cohesion << "\t" << cell[i][j].curl << "\t" << cell[i][j].omega << "\t" << cell[i][j].W << endl;
		data_file << endl;
	}
}

void Field::Save_Theta_Deviation(ofstream& data_file)
{
	for (int i = 0; i < grid_dim_x; i++)
	{
		for (int j = 0; j < grid_dim_x; j++)
			for (int n = 0; n < cell[i][j].dtheta.size(); n++)
				data_file << std::setprecision(20) << cell[i][j].dtheta[n] << endl;
	}
}

void Field::Draw(string info)
{
	Gnuplot gp;

	ofstream temp_file("data.dat");

////	info="phi=".phi."-g=".g."-noise=".noise
////	filename="data-".info.".dat"

//	vector<boost::tuple<double, double, double, double, double, double, double, double> > temp;
//	vector<boost::tuple<double, double, double, double, double, double, double, double> > temp;
//	vector<boost::tuple<double, double, double, double> > pts;

	for (int x = 0; x < grid_dim_x; x++)
	{
		for (int y = 0; y < grid_dim_x; y++)
		{
//			pts.push_back(boost::make_tuple(cell[x][y].r.x,cell[x][y].r.y,cell[x][y].v.x,cell[x][y].v.y,cell[x][y].density,cell[x][y].cohesion,cell[x][y].curl,cell[x][y].omega));
//			pts.push_back(boost::make_tuple(cell[x][y].r.x,cell[x][y].r.y,cell[x][y].v.x,cell[x][y].v.y));
			temp_file << std::fixed << std::setprecision(8) << cell[x][y].r / (2*L) << "\t" << cell[x][y].v << "\t" << cell[x][y].density << "\t" << cell[x][y].cohesion << "\t" << cell[x][y].curl << "\t" << cell[x][y].omega << "\t" << cell[x][y].W << endl;
		}
			temp_file << endl;
	}

	gp << "L=" << 0.5 << "\n";

	gp << "set term postscript eps enhanced color\n";
	gp << "reset\n";
	gp << "l=2\n";


	gp << "set size square\n";

	gp << "set term postscript eps enhanced color font \"Times-Roman,30\"\n";

//	gp << "set contour base\n";
//	gp << "set cntrparam level incremental 0.4, 0.5\n";

	gp << "set output \"figures/" << info << "-density.eps\"\n";
	gp << "set pm3d map\n";
	gp << "set palette rgb 21,22,23\n";
	gp << "set cbrange [0:3]\n";
	gp << "set colorbox vertical size 0.1,0.4\n";
	gp << "set colorbox vertical user origin 0.72,0.28 size 0.035,0.5\n";
	gp << "set cbtics offset -1,0\n";
	gp << "set cblabel \"{/Symbol r}\" offset -6.5,3.5\n";
//	gp << "unset colorbox\n";
	gp << "set xrange [-L:L]\n";
	gp << "set yrange [-L:L]\n";
	gp << "set xlabel \"x/L\" offset 0,1.2\n";
//	gp << "set ylabel \"y/L\" offset 3,0\n";
	gp << "set xtics offset 0,0.7\n";
	gp << "set ytics format \" \" \n";
	gp << "set ytics offset 0.7,0\n";
	gp << "splot \"data.dat\" using 1:2:5 notitle\n";
//	gp.send1d(pts);




	gp << "set output \"figures/" << info << "-cohesion.eps\"\n";
	gp << "set pm3d map\n";
	gp << "set cbrange [0:1]\n";
	gp << "set cblabel \"{/Symbol F}\" offset -6.2,3.5\n";
	gp << "splot \"data.dat\" using 1:2:6 notitle\n";
//	gp.send1d(pts);

	gp << "set term postscript eps enhanced color font \"Times-Roman,35\"\n";
	gp << "set palette rgb 0,3,21\n";
	gp << "l=0.035\n";
	gp << "set output \"figures/"<< info << "-velocity-vector.eps\"\n";
	gp << "set xtics offset 0,0.4\n";
	gp << "set xlabel \"x/L\" offset 0,0.8\n";
	gp << "set ylabel \"y/L\" offset 4.5,0\n";
	gp << "unset ylabel\n";
	gp << "set ytics format \" \" \n";
	gp << "set colorbox vertical user origin 0.8,0.24 size 0.04,0.63\n";
	gp << "set cblabel \"v\" offset -6.4,4.2\n";
//	gp << "unset colorbox\n";
	gp << "set style arrow 1 head filled size screen 0.008,25,30 lw 2 lc palette\n";
	gp << "plot \"data.dat\" using 1:2:(l*($3)):(l*($4)):(($3**2 + $4**2)**0.5) with vectors arrowstyle 1 notitle\n";
//	gp.send1d(pts);

	gp << "l=1\n";
	gp << "set output \"figures/"<< info << "-W.eps\"\n";
	gp << "set xrange [-L:L]\n";
	gp << "set yrange [-L:L]\n";
	gp << "plot \"data.dat\" using 1:2:(l*($9)):(l*($10)):(($9**2 + $10**2)**0.5) with vectors arrowstyle 1 notitle\n";
//	gp.send1d(pts);
}

void Field::Draw_Section(string info, int y)
{
	Gnuplot gp;

	ofstream temp_file("data.dat");

////	info="phi=".phi."-g=".g."-noise=".noise
////	filename="data-".info.".dat"

//	vector<boost::tuple<double, double, double, double, double, double, double, double> > temp;
//	vector<boost::tuple<double, double, double, double, double, double, double, double> > temp;
//	vector<boost::tuple<double, double, double, double> > pts;

	for (int x = 0; x < grid_dim_x; x++)
	{

		temp_file << std::fixed << std::setprecision(8) << cell[x][y].r.x / (2*L) << "\t" << cell[x][y].v << "\t" << cell[x][y].density << "\t" << cell[x][y].cohesion << "\t" << cell[x][y].curl << "\t" << cell[x][y].omega << "\t" << cell[x][y].W << endl;
			temp_file << endl;
	}

	gp << "reset\n";

	gp << "L=" << 0.5 << "\n";

	gp << "set term postscript eps enhanced color \n";
//	gp << "l=2\n";

	gp << "set size square\n";
//	gp << "set size 0.85,1.05\n";

//	gp << "set lmargin at screen 0.05\n";
//	gp << "set rmargin at screen 0.75\n";
//	gp << "set bmargin at screen 0.05\n";
//	gp << "set tmargin at screen 1.1\n";

	gp << "set term postscript eps enhanced color font \"Times-Roman,35\"\n";

	gp << "set style line 1 lc rgb '#FF4500' lt 1 lw 2 pt 7 pi 0 ps 1.5\n";
	gp << "set style line 2 lc rgb '#008000' lt 1 lw 2 pt 4 pi 0 ps 1.5\n";
	gp << "set pointintervalbox 3\n";

	if (fabs(y) > 20)
		gp << "set key right top\n";
	else
		gp << "set key center top\n";

	gp << "set log y\n";
	gp << "set xlabel \"x/L\" offset 0,1 \n";
	gp << "set xtics offset 0,0.4 \n";
	gp << "set ytics offset 0.5,0 \n";
	gp << "set ylabel \"{/Symbol r}, {/Symbol f}\" offset 3,0 \n";
	gp << "set output \"figures/" << info << "-density-section-y=" << std::fixed << std::setprecision(0) << round((y+0.5)*Field_Cell::dim_x - L)  << ".eps\"\n";
	gp << "set xrange [-L:L]\n";
	gp << "plot \"data.dat\" using 1:4 w lp ls 1 ti \"{/Symbol r}\", \"data.dat\" using 1:5 w lp ls 2 ti \"{/Symbol f}\"\n";

	gp << "unset log y\n";
//	gp.send1d(pts);

//	gp << "set ylabel \"{/Symbol f}\"\n";
//	gp << "set output \"figures/" << info << "-cohesion-section.eps\"\n";
//	gp << "set xrange [-L:L]\n";
//	gp << "plot \"data.dat\" using 1:5 w lp ls 1\n";
//	gp.send1d(pts);

//	gp << "set output \"" << info << "-curl-section.eps\"\n";
//	gp << "plot \"data.dat\" using 1:6 w lp ls 1\n";
//	gp.send1d(pts);

	gp << "unset log \n";
	gp << "set ylabel \"v\" offset 3,0 \n";
	gp << "set output \"figures/"<< info << "-velocity-section-y=" << std::fixed << round((y+0.5)*Field_Cell::dim_x - L)  << ".eps\"\n";
	gp << "set xrange [-L:L]\n";
	gp << "plot \"data.dat\" using 1:3 w lp ls 1 ti \"v_y\", \"data.dat\" using 1:2 w lp ls 2 ti \"v_x\" \n";
//	gp.send1d(pts);

	gp << "set ylabel \"W\" offset 3,0 \n";
	gp << "set output \"figures/"<< info << "-W-section-y=" << std::fixed << std::setprecision(0) << round((y+0.5)*Field_Cell::dim_x - L)  << ".eps\"\n";
	gp << "set xrange [-L:L]\n";
	gp << "plot \"data.dat\" using 1:9 w lp ls 1 ti \"W_y\", \"data.dat\" using 1:8 w lp ls 2 ti \"W_x\" \n";
//	gp.send1d(pts);

	gp << "set nokey\n";
	gp << "set ylabel \"{/Symbol w}\" offset 3,0 \n";
	gp << "set output \"figures/"<< info << "-omega-section-y=" << std::fixed << std::setprecision(0) << round((y+0.5)*Field_Cell::dim_x - L)  << ".eps\"\n";
	gp << "set xrange [-L:L]\n";
	gp << "plot \"data.dat\" using 1:(-$3/$1) w lp ls 1\n";
//	gp.send1d(pts);


}


void Field::Draw_Density_Contour(string info, Real rho)
{

	C2DVector r[grid_dim_x];
	for (int x = 0; x < grid_dim_x; x++)
		for (int y = (grid_dim_x-3); y > 1; y--)
		{
			if (fabs(log(cell[x][y-1].density / cell[x][y].density)/(L/grid_dim_x)) > rho)
			{
				r[x] = cell[x][y].r / (2*L);
				r[x].y = (0.5 - r[x].y);
				r[x].x += 0.5;
				y = 0;
			}
		}

	for (int x = 0; x < grid_dim_x; x++)
		if (fabs(r[x].y) > 0.001)
			cout << std::fixed << std::setprecision(8) << r[x] << endl;
}

void Field::Reset()
{
	for (int x = 0; x < grid_dim_x; x++)
		for (int y = 0; y < grid_dim_x; y++)
			cell[x][y].Reset();
	sample = 0;
}

void Field::Add(Field* f)
{
	if (grid_dim_x == f->grid_dim_x)
	{
		for (int x = 0; x < grid_dim_x; x++)
		{
			for (int y = 0; y < grid_dim_x; y++)
			{
				cell[x][y].dtheta.insert(cell[x][y].dtheta.end(), f->cell[x][y].dtheta.begin(), f->cell[x][y].dtheta.end());
				cell[x][y].v += f->cell[x][y].v;
				cell[x][y].W += f->cell[x][y].W;
				cell[x][y].theta_ave += f->cell[x][y].theta_ave;
				cell[x][y].density += f->cell[x][y].density;
				cell[x][y].cohesion += f->cell[x][y].cohesion;
				cell[x][y].omega += f->cell[x][y].omega;
				cell[x][y].curl += f->cell[x][y].curl;
			}
		}
		sample++;
	}
	else
		cout << "Error! Size of two fields are different. The program can't add two fields with different sizes";
}

void Field::Average()
{
	for (int x = 0; x < grid_dim_x; x++)
	{
		for (int y = 0; y < grid_dim_x; y++)
		{
			cell[x][y].v /= sample;
			cell[x][y].W /= sample;
			cell[x][y].theta_ave /= sample;
			cell[x][y].density /= sample;
			cell[x][y].cohesion /= sample;
			cell[x][y].omega /= sample;
			cell[x][y].curl /= sample;
		}
	}
}

#endif
