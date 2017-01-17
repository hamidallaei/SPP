#include <boost/algorithm/string.hpp>

#include<iostream>
#include<cstdlib>
#include<vector>

#include"read.h"
#include"statistics.h"
#include"field.h"
#include"pair-set.h"

using namespace std;

double Local_Cohesion(SceneSet* s, double rc)
{
	long double phi = 0;
	long int counter = 0;
	for (int i = 0; i < s->scene.size(); i++)
	{
		for (int j = 0; j < Scene::Ns; j++)
			for (int k = j+1; k < Scene::Ns; k++)
				if ((s->scene[i].sparticle[j].r - s->scene[i].sparticle[k].r).Square() < (rc*rc))
				{
					phi += cos(s->scene[i].sparticle[j].theta - s->scene[i].sparticle[k].theta);
					counter++;
				}
	}
	phi /= counter;
	return(phi);
}

void Polarization_AutoCorr(SceneSet* s)
{
	Stat<double> p;
	for (int i = 0; i < s->scene.size(); i++)
	{
		C2DVector vp;
		vp.Null();
		for (int j = 0; j < Scene::Ns; j++)
		{
			C2DVector temp_vec;
			temp_vec.x = cos(s->scene[i].sparticle[j].theta);
			temp_vec.y = sin(s->scene[i].sparticle[j].theta);
			vp += temp_vec;
		}
		vp = vp / Scene::Ns;
		p.Add_Data(sqrt(vp.Square()));
	}
	p.Compute();
	p.Correlation();
}

void Polarization_Time(SceneSet* s)
{
	for (int i = 0; i < s->scene.size(); i++)
	{
		C2DVector p;
		p.Null();
		for (int j = 0; j < Scene::Ns; j++)
		{
			C2DVector temp_vec;
			temp_vec.x = cos(s->scene[i].sparticle[j].theta);
			temp_vec.y = sin(s->scene[i].sparticle[j].theta);
			p += temp_vec;
		}
		p = p / Scene::Ns;
		cout << i << "\t" << sqrt(p.Square()) << endl;
	}
}

void Quantities_Time(SceneSet* s, const int start_number, ostream& os)
{
	cout << "time\tp\tS\tdr2" << endl;
	for (int i = start_number; i < s->scene.size(); i++)
	{
		C2DVector p;
		C2DVector dr;
		p.Null();
		double c2 = 0;
		double s2 = 0;
		double dr2 = 0;
		for (int j = 0; j < Scene::Ns; j++)
		{
			C2DVector temp_vec;
			temp_vec.x = cos(s->scene[i].sparticle[j].theta);
			temp_vec.y = sin(s->scene[i].sparticle[j].theta);
			dr = s->scene[i].sparticle[j].r - s->scene[start_number].sparticle[j].r;
			p += temp_vec;
			c2 += 2*temp_vec.x*temp_vec.x - 1;
			s2 += 2*temp_vec.x*temp_vec.y;
			dr2 += dr.Square();
		}
		p = p / Scene::Ns;
		c2 /= Scene::Ns;
		s2 /= Scene::Ns;
		dr2 /= Scene::Ns;
		double S = sqrt(c2*c2 + s2*s2);
		os << s->scene[i].t << "\t" << sqrt(p.Square()) << "\t" << S << "\t" << dr2 << endl;
	}
}

void Compute_Polarization(SceneSet* s, Stat<double>* polarization)
{
	for (int i = 0; i < s->scene.size(); i++)
	{
		C2DVector p;
		p.Null();
		for (int j = 0; j < Scene::Ns; j++)
		{
			C2DVector temp_vec;
			temp_vec.x = cos(s->scene[i].sparticle[j].theta);
			temp_vec.y = sin(s->scene[i].sparticle[j].theta);
			p += temp_vec;
		}
		p = p / Scene::Ns;
		polarization->Add_Data(sqrt(p.Square()));
	}
	polarization->Compute();
}

void Compute_Order_Parameters(SceneSet* s, double& polarization, double& error_polarization, double& sigma2, double& G)
{
	Stat<double> p,p4;
	for (int i = 0; i < s->scene.size(); i++)
	{
		C2DVector vp;
		vp.Null();
		for (int j = 0; j < Scene::Ns; j++)
		{
			C2DVector temp_vec;
			temp_vec.x = cos(s->scene[i].sparticle[j].theta);
			temp_vec.y = sin(s->scene[i].sparticle[j].theta);
			vp += temp_vec;
		}
		vp = vp / Scene::Ns;
		p.Add_Data(sqrt(vp.Square()));
		p4.Add_Data(vp.Square()*vp.Square());
	}
	p.Compute();
	p4.Compute();
	polarization = p.mean;
	sigma2 = p.variance;
	error_polarization = p.error;
	sigma2 *= (4*s->L.x*s->L.y);
	G = 1 - (p4.mean / (3*p.mean_square));
}

void Compute_Angular_Momentum(SceneSet* s, Stat<double>* angular_momentum)
{
	for (int i = 0; i < s->scene.size(); i++)
	{
		double M = 0;
		for (int j = 0; j < Scene::Ns; j++)
		{
			C2DVector temp_vec;
			temp_vec.x = cos(s->scene[i].sparticle[j].theta);
			temp_vec.y = sin(s->scene[i].sparticle[j].theta);
			M += (s->scene[i].sparticle[j].r.x * temp_vec.y - s->scene[i].sparticle[j].r.y * temp_vec.x);
		}
		M /= Scene::Ns;
		angular_momentum->Add_Data(M);
	}
	angular_momentum->Compute();
}

double Time_AutoCorrelationPoint(SceneSet* s, int tau)
{
	double result = 0;
		for (int i = 0; i < s->scene.size() - tau; i++)
			for (int j = 0; j < Scene::Ns; j++)
				result += cos(s->scene[i].sparticle[j].theta - s->scene[i+tau].sparticle[j].theta);
	result /= (Scene::Ns * (s->scene.size() - tau));
	return result;
}

void Time_AutoCorrelation(SceneSet* s, int step)
{
	for (int i = 0; i < (s->scene.size() - 5); i+=step)
		cout << i << "\t" << Time_AutoCorrelationPoint(s, i) << endl;
}

void Spatial_AutoCorrelation(SceneSet* s, int size, double rc)
{
	double bin[size];
	int num[size];

	for (int x = 0; x < size; x++)
	{
		bin[x] = 0;
		num[x] = 0;
	}

	for (int i = s->scene.size()/2; i < s->scene.size(); i+=100)
	{
		for (int j = 0; j < Scene::Ns; j++)
			for (int k = j+1; k < Scene::Ns; k++)
			{
				C2DVector dr = s->scene[i].sparticle[j].r - s->scene[i].sparticle[k].r;
				double r = sqrt(dr.Square());
				if (r < rc)
				{
					int x = (int) round(size*(r / rc));
					num[x]++;
					bin[x] += cos(s->scene[i].sparticle[j].theta - s->scene[i].sparticle[k].theta);
				}
			}
	}

	for (int x = 1; x < size; x++)
	{
		double r = (x*rc)/size;
//		bin[x] /= Scene::Ns*((s->scene.size() - s->scene.size()/2)/100);
//		bin[x] /= 3.1415*((r+rc/ size)*(r+rc/ size) - r*r)/2;
//		bin[x] -= Scene::Ns / (4*s->L.x*s->L.y);
		if (num[x] != 0)
		  bin[x] /= num[x];
		else
		  bin[x] = 0;
		cout << r << "\t" << bin[x] << endl;
	}
}


void Trajectory(SceneSet* s, int index)
{
	for (int i = 0; i < s->scene.size(); i++)
		cout << s->scene[i].sparticle[index].r << endl;
}

double Find_Angle(C2DVector r)
{
	double angle = atan2(r.x , r.y);
//	if (r.x < 0)
	//	angle += 3.1415;
	return (angle);
}

void Angle_Time(SceneSet* s, int index)
{
	for (int i = 0; i < (s->scene.size() - 1); i++)
		cout << 4*i << "\t" << Find_Angle(s->scene[i].sparticle[index].r) << endl;
}

double Find_Angular_Velocity(SceneSet* s, int index, int t)
{
	float dt = 1.0;
	if ((t < s->scene.size()) && (t > 0))
	{
		double dtheta = Find_Angle(s->scene[t].sparticle[index].r);
		dtheta = dtheta - Find_Angle(s->scene[t-1].sparticle[index].r);
		if (abs(dtheta) > 3)
			dtheta = Find_Angle(s->scene[t+1].sparticle[index].r) - Find_Angle(s->scene[t].sparticle[index].r);
		dtheta /= dt;
		return dtheta;
	}
	else
		cout << "error\t time is out of range" << endl;
}

void Angular_Velocity_Time(SceneSet* s, int index)
{
	static bool b = true;
	if (b)
	{
		cout << "rho	noise	x	y	r	omega" << endl;
		b = false;
	}
	for (int i = 1; i < (s->scene.size() - 1); i+=100)
		cout << Scene::density << "\t" << Scene::noise << "\t" << s->scene[i].sparticle[index].r << "\t" << sqrt(s->scene[i].sparticle[index].r.Square()) << "\t" << Find_Angular_Velocity(s, index, i) << endl;
}

void Window_Fluctuation(SceneSet* s, int smaller_number_of_windows, double& mean, double& variance)
{
	int number_of_windows_x, number_of_windows_y;
	if (s->L.x > s->L.y)
	{
		number_of_windows_y = smaller_number_of_windows;
		number_of_windows_x = (int) round(s->L.x*smaller_number_of_windows / s->L.y);
	}
	else
	{
		number_of_windows_x = smaller_number_of_windows;
		number_of_windows_y = (int) round(s->L.y*smaller_number_of_windows / s->L.x);
	}
	Stat<int> window[number_of_windows_x][number_of_windows_y];
	for (int i = s->scene.size()/2; i < s->scene.size(); i++)
	{
		int Np[number_of_windows_x][number_of_windows_y];
		for (int x = 0; x < number_of_windows_x; x++)
			for (int y = 0; y < number_of_windows_y; y++)
				Np[x][y] = 0;
		for (int j = 0; j < Scene::Ns; j++)
		{
			int x = (int) floor(number_of_windows_x*(s->scene[i].sparticle[j].r.x / s->L.x + 1)/2);
			int y = (int) floor(number_of_windows_y*(s->scene[i].sparticle[j].r.y / s->L.y + 1)/2);
			Np[x][y]++;
		}
		for (int x = 0; x < number_of_windows_x; x++)
			for (int y = 0; y < number_of_windows_y; y++)
				window[x][y].Add_Data(Np[x][y]);
	}
	for (int x = 0; x < number_of_windows_x; x++)
		for (int y = 0; y < number_of_windows_y; y++)
			window[x][y].Compute();

	mean = variance = 0;
	for (int x = 0; x < number_of_windows_x; x++)
		for (int y = 0; y < number_of_windows_y; y++)
		{
			mean += window[x][y].mean;
			variance += window[x][y].variance;
		}
	mean /= (number_of_windows_x*number_of_windows_y);
	variance /= (number_of_windows_x*number_of_windows_y);
}

double Compute_Fluctuation(SceneSet* s)
{
	double mean, variance;
	for (int i = 5; i < s->L.y; i = (int) (i*1.3))
	{
		Window_Fluctuation(s, i, mean, variance);
		cout << mean << "\t" << variance << endl;
	}
}

void Angle_Distribution(SceneSet* s)
{
}

// Find radial density.
void Radial_Density(SceneSet* s, int number_of_points)
{
	double* rho = new double[number_of_points];
	for (int i = 0; i < number_of_points; i++)
		rho[i] = 0;
	int counter = 0;
	double radius[number_of_points];
	radius[0] = 1;
	double factor = pow(((s->L.x+1)/radius[0]-0),1.0/number_of_points);
	for (int i = 1; i < number_of_points; i++)
		radius[i] = factor*radius[i-1];
	for (int i = 0; i < s->scene.size(); i++)
	{
		counter++;
		for (int j = 0; j < Scene::Ns; j++)
		{
			Real r = sqrt(s->scene[i].sparticle[j].r.Square());
			int index = (int) (log(r/radius[0]) / log(factor));
			if (index >= 0)
				rho[index]++;
		}
	}
	for (int i = 1; i < number_of_points; i++)
	{
		rho[i] /= M_PI*(radius[i]*radius[i] - radius[i-1]*radius[i-1]);
		rho[i] /= s->scene.size();
		cout << sqrt(radius[i-1]*radius[i]) << "\t" << rho[i] << endl;
	}

	delete [] rho;
}

// Find radial density.
void Pair_Distribution(SceneSet* s, Real lx, Real ly,int smaller_grid_size)
{
	int grid_size_x, grid_size_y;
	if (s->L.x > s->L.y)
	{
		grid_size_y = smaller_grid_size;
		grid_size_x = (int) round(s->L.x*smaller_grid_size / s->L.y);
	}
	else
	{
		grid_size_x = smaller_grid_size;
		grid_size_y = (int) round(s->L.y*smaller_grid_size / s->L.x);
	}

	double bin[grid_size_x][grid_size_y];
	int counter = 0;

	for (int x = 0; x < grid_size_x; x++)
		for (int y = 0; y < grid_size_y; y++)
		{
			bin[x][y] = 0;
		}

	for (int i = s->scene.size()/2; i < s->scene.size(); i+=100)
	{
		for (int j = 0; j < Scene::Ns; j++)
			for (int k = 0; k < Scene::Ns; k++)
			{
				if (j != k)
				{
					C2DVector temp_vec;
					temp_vec.x = cos(s->scene[i].sparticle[j].theta);
					temp_vec.y = sin(s->scene[i].sparticle[j].theta);
					C2DVector dr = s->scene[i].sparticle[k].r - s->scene[i].sparticle[j].r;
					C2DVector tdr = dr;
					tdr.y = temp_vec * dr;
					tdr.x = (temp_vec.y * dr.x) - (temp_vec.x * dr.y);
					if (fabs(tdr.x) < lx && fabs(tdr.y) < ly)
					{
						int x = (int) (grid_size_x*((tdr.x / lx) + 1)/2);
						int y = (int) (grid_size_y*((tdr.y / ly) + 1)/2);
						bin[x][y]++;
					}
				}
			}
		counter++;
	}

	for (int x = 0; x < grid_size_x; x++)
	{
		for (int y = 0; y < grid_size_y; y++)
		{
			bin[x][y] /= (Scene::Ns*(4*lx*ly/grid_size_x/grid_size_y));
			bin[x][y] /= counter;
			cout << ((2.0*x)/grid_size_x-1)*lx << "\t" << ((2.0*y)/grid_size_y-1)*ly << "\t" << bin[x][y] << endl;
		}
		cout << endl;
	}
}

// Find distance growth in time (Diffusion)
void Mean_Squared_Distance_Growth(SceneSet* s, int frames, int number_of_points, int number_of_pair_sets, Real r_cut)
{
	int interval = (s->scene.size() - frames) / number_of_pair_sets;
	Pair_Set ps[number_of_pair_sets];
	Pair_Set::sceneset = s;
	for (int i = 0; i < number_of_pair_sets; i++)
		ps[i].Find_Close_Particle(r_cut, i*interval);

	Real* md2 = new Real[number_of_points];
	int tau[number_of_points];
	for (int i = 0; i < number_of_points; i++)
	{
		md2[i] = 0;
		tau[i] = i*frames / (number_of_points);
		for (int j = 0; j < number_of_pair_sets; j++)
			md2[i] += ps[j].Find_Mean_Square_Distance(tau[i]);
		md2[i] / number_of_pair_sets;
	}
	for (int i = 1; i < number_of_points; i++)
		cout << tau[i] << "\t" << md2[i] - md2[0] << endl;

	delete [] md2;
}

// Find diffusion 
void Mean_Squared_Displacement_Growth(SceneSet* s, int frames, int number_of_points)
{
	Real* md2 = new Real[number_of_points];
	int tau[number_of_points];
	int N = s->scene[0].Ns;
	for (int i = 0; i < number_of_points; i++)
	{
		md2[i] = 0;
		tau[i] = i*frames / (number_of_points);
		for (int j = 0; j < N; j++)
			md2[i] += (s->scene[tau[i]].sparticle[j].r - s->scene[0].sparticle[j].r).Square();
		md2[i] /= N;
	}
	for (int i = 1; i < number_of_points; i++)
		cout << tau[i] << "\t" << md2[i] << endl;

	delete [] md2;
}


// Find distance growth in time (Lyapanov)
bool Lyapunov_Exponent(SceneSet* s, int frames, int number_of_points, int number_of_pair_sets, Real r_min, Real r_max)
{
	int total_frame = s->scene.size();
	int interval = (total_frame - frames) / number_of_pair_sets;
	Pair_Set ps[number_of_pair_sets];
	Pair_Set::sceneset = s;
	if (interval < 0)
		return(false);
	for (int i = 0; i < number_of_pair_sets; i++)
		ps[i].Find_Particle(r_min, r_max, i*interval);

	Real* lambda = new Real[number_of_points];
	int tau[number_of_points];

	for (int i = 0; i < number_of_points; i++)
	{
		lambda[i] = 0;
		tau[i] = i*frames / (number_of_points);
		for (int j = 0; j < number_of_pair_sets; j++)
			lambda[i] += ps[j].Find_Short_Lyapunov_Exponent(tau[i]);
		lambda[i] / number_of_pair_sets;
	}
	for (int i = 1; i < number_of_points; i++)
		cout << tau[i] << "\t" << lambda[i] << endl;

	delete [] lambda;
	return(true);
}

