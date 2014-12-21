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
	for (int i = 0; i < s->scene.size(); i++)
	{
		double L = 0;
		for (int j = 0; j < Scene::number_of_particles; j++)
			L += (s->scene[i].particle[j].r.x * s->scene[i].particle[j].v.y - s->scene[i].particle[j].r.y * s->scene[i].particle[j].v.x);
		angular_momentum->Add_Data(L);
	}
	angular_momentum->Compute();
}

double AutoCorrelationPoint(SceneSet* s, int tau)
{
	double result = 0;
		for (int i = 0; i < s->scene.size() - tau; i++)
			for (int j = 0; j < Scene::number_of_particles; j++)
				result += (s->scene[i].particle[j].v * s->scene[i+tau].particle[j].v);
	result /= (Scene::number_of_particles * (s->scene.size() - tau));
	return result;
}

void AutoCorrelation(SceneSet* s)
{
	for (int i = 0; i < (s->scene.size() - 5); i++)
		cout << i*4 << "\t" << AutoCorrelationPoint(s, i) << endl;
}
