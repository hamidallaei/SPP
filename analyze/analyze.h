#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#include <boost/algorithm/string.hpp>

#include<iostream>
#include<cstdlib>
#include<vector>

#include"read.h"
#include"statistics.h"
#include"field.h"

using namespace std;

void Compute_Polarization(SceneSet* s, Stat<double>* polarization)
{
	double L = 0;
	for (int i = 0; i < s->scene.size(); i++)
	{
		C2DVector p;
		p.Null();
		for (int j = 0; j < Scene::number_of_particles; j++)
			p += s->scene[i].particle[j].v;
		p = p / Scene::number_of_particles;
		polarization->Add_Data(sqrt(p.Square()));
	}
	polarization->Compute();
}

void Compute_Angular_Momentum(SceneSet* s, Stat<double>* angular_momentum)
{
	for (int i = s->scene.size()/2; i < s->scene.size(); i++)
	{
		double L = 0;
		for (int j = 0; j < Scene::number_of_particles; j++)
			L += (s->scene[i].particle[j].r.x * s->scene[i].particle[j].v.y - s->scene[i].particle[j].r.y * s->scene[i].particle[j].v.x);
		L /= Scene::number_of_particles;
		angular_momentum->Add_Data(L);
	}
	angular_momentum->Compute();
}

double Time_AutoCorrelationPoint(SceneSet* s, int tau)
{
	double result = 0;
		for (int i = 0; i < s->scene.size() - tau; i++)
			for (int j = 0; j < Scene::number_of_particles; j++)
				result += (s->scene[i].particle[j].v * s->scene[i+tau].particle[j].v);
	result /= (Scene::number_of_particles * (s->scene.size() - tau));
	return result;
}

void Time_AutoCorrelation(SceneSet* s)
{
	for (int i = 0; i < (s->scene.size() - 5); i+=20)
		cout << i << "\t" << Time_AutoCorrelationPoint(s, i) << endl;
}

void Spatial_AutoCorrelation(SceneSet* s, int size, double rc)
{
	double bin[size];

	for (int x = 0; x < size; x++)
	{
		bin[x] = 0;
	}

	for (int i = s->scene.size()/2; i < s->scene.size(); i+=100)
	{
		for (int j = 0; j < Scene::number_of_particles; j++)
			for (int k = j+1; k < Scene::number_of_particles; k++)
			{
				C2DVector dr = s->scene[i].particle[j].r - s->scene[i].particle[k].r;
				double r = sqrt(dr.Square());
				if (r < rc)
				{
					int x = (int) round(size*(r / rc));
						bin[x]++;
				}
			}
	}

	for (int x = 1; x < size; x++)
	{
		double r = (x*rc)/size;
		bin[x] /= Scene::number_of_particles*((s->scene.size() - s->scene.size()/2)/100);
		bin[x] /= 3.1415*((r+rc/ size)*(r+rc/ size) - r*r)/2;
		bin[x] -= Scene::number_of_particles / (4*s->L*s->L);
		cout << r << "\t" << bin[x] << endl;
	}
}


void Trajectory(SceneSet* s, int index)
{
	for (int i = 0; i < s->scene.size(); i++)
		cout << s->scene[i].particle[index].r << endl;
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
		cout << 4*i << "\t" << Find_Angle(s->scene[i].particle[index].r) << endl;
}

double Find_Angular_Velocity(SceneSet* s, int index, int t)
{
	float dt = 1.0;
	if ((t < s->scene.size()) && (t > 0))
	{
		double dtheta = Find_Angle(s->scene[t].particle[index].r);
		dtheta = dtheta - Find_Angle(s->scene[t-1].particle[index].r);
		if (abs(dtheta) > 3)
			dtheta = Find_Angle(s->scene[t+1].particle[index].r) - Find_Angle(s->scene[t].particle[index].r);
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
		cout << Scene::density << "\t" << Scene::noise << "\t" << s->scene[i].particle[index].r << "\t" << sqrt(s->scene[i].particle[index].r.Square()) << "\t" << Find_Angular_Velocity(s, index, i) << endl;
}

void Window_Fluctuation(SceneSet* s, int number_of_windows, double& mean, double& variance)
{
	Stat<int> window[number_of_windows][number_of_windows];
	for (int i = s->scene.size()/2; i < s->scene.size(); i++)
	{
		int Np[number_of_windows][number_of_windows];
		for (int x = 0; x < number_of_windows; x++)
			for (int y = 0; y < number_of_windows; y++)
				Np[x][y] = 0;
		for (int j = 0; j < Scene::number_of_particles; j++)
		{
			int x = (int) floor(number_of_windows*(s->scene[i].particle[j].r.x / s->scene[i].L + 1)/2);
			int y = (int) floor(number_of_windows*(s->scene[i].particle[j].r.y / s->scene[i].L + 1)/2);
			Np[x][y]++;
		}
		for (int x = 0; x < number_of_windows; x++)
			for (int y = 0; y < number_of_windows; y++)
				window[x][y].Add_Data(Np[x][y]);
	}
	for (int x = 0; x < number_of_windows; x++)
		for (int y = 0; y < number_of_windows; y++)
			window[x][y].Compute();

	mean = variance = 0;
	for (int x = 0; x < number_of_windows; x++)
		for (int y = 0; y < number_of_windows; y++)
		{
			mean += window[x][y].mean;
			variance += window[x][y].variance;
		}
	mean /= (number_of_windows*number_of_windows);
	variance /= (number_of_windows*number_of_windows);
}

double Compute_Fluctuation(SceneSet* s)
{
	double mean, variance;
	for (int i = 5; i < s->L; i = (int) (i*1.3))
	{
		Window_Fluctuation(s, i, mean, variance);
		cout << mean << "\t" << variance << endl;
	}
}


