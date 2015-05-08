#ifndef _STATE_HYPER_VECTOR_
#define _STATE_HYPER_VECTOR_

#include "c2dvector.h"

class State_Hyper_Vector{
public:
	int N;
	BasicParticle0* particle;
	Real Square() const;
	Real Magnitude() const;
	State_Hyper_Vector operator+ (const State_Hyper_Vector& s1) const;
	State_Hyper_Vector operator- (const State_Hyper_Vector& s1) const;
	Real operator* (const State_Hyper_Vector& s1) const;
};


State_Hyper_Vector State_Hyper_Vector::operator+ (const State_Hyper_Vector& s1) const
{
	State_Hyper_Vector result;

	for (int i = 0; i < N; i++)
	{
		result.particle[i].r = (particle[i].r + s1.particle[i].r);
		result.particle[i].r.Periodic_Transform();
		result.particle[i].theta = particle[i].theta + s1.particle[i].theta;
		result.particle[i].theta -= (int) (result.particle[i].theta / (2*M_PI));
	}
	return result;
}


State_Hyper_Vector State_Hyper_Vector::operator- (const State_Hyper_Vector& s1) const
{
	State_Hyper_Vector result;

	for (int i = 0; i < N; i++)
	{
		result.particle[i].r = (particle[i].r - s1.particle[i].r);
		result.particle[i].r.Periodic_Transform();
		result.particle[i].theta = particle[i].theta - s1.particle[i].theta;
		result.particle[i].theta -= (int) (result.particle[i].theta / (2*M_PI));
	}
	return result;
}

Real State_Hyper_Vector::operator* (const State_Hyper_Vector& s1) const
{
	Real result = 0;

	for (int i = 0; i < N; i++)
	{
		result += particle[i].r * s1.particle[i].r;
		result += cos(particle[i].theta - s1.particle[i].theta);
	}
	return result;
}

Real State_Hyper_Vector::Square() const
{
	Real result = 0;

	for (int i = 0; i < N; i++)
	{
		result += particle[i].r * particle[i].r;
		result += 1;
	}
	return result;
}

Real State_Hyper_Vector::Magnitude() const
{
	Real result = 0;

	for (int i = 0; i < N; i++)
	{
		result += particle[i].r * particle[i].r;
		result += 1;
	}
	return sqrt(result);
}

#endif
