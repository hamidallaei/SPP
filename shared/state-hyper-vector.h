#ifndef _STATE_HYPER_VECTOR_
#define _STATE_HYPER_VECTOR_

#include "c2dvector.h"

class State_Hyper_Vector{
	void Init_Random_Generator();
	~State_Hyper_Vector();
public:
	gsl_rng* gsl_r;
	int N;
	BasicParticle0* particle;
	
	State_Hyper_Vector(int);
	
	State_Hyper_Vector& operator= ( const State_Hyper_Vector& sv);
	State_Hyper_Vector& operator+ (const State_Hyper_Vector& s1) const;
	State_Hyper_Vector& operator- (const State_Hyper_Vector& s1) const;
	Real operator* (const State_Hyper_Vector& s1) const;

	void Rand(const Real position_amplitude, const Real angle_amplitude);
	Real Square() const;
	Real Magnitude() const;
};

void State_Hyper_Vector::Init_Random_Generator()
{
	const gsl_rng_type * T;
	gsl_rng_env_setup();
	gsl_rng_default_seed = 0;
	T = gsl_rng_default;
	gsl_r = gsl_rng_alloc (T);
	gsl_rng_memcpy (gsl_r, C2DVector::gsl_r);
}

State_Hyper_Vector::State_Hyper_Vector(int particle_number) : N(particle_number)
{
	Init_Random_Generator();
	particle = new BasicParticle0[N];
}

State_Hyper_Vector::State_Hyper_Vector(const State_Hyper_Vector& sv) : N(sv.N)
{
	Init_Random_Generator();
	gsl_rng_memcpy (gsl_r, sv.gsl_r);
	particle = new BasicParticle0[N];
	for (int i = 0; i < N; i++)
		particle[i] = sv.particle[i];
}

State_Hyper_Vector::~State_Hyper_Vector()
{
	gsl_rng_free(gsl_r);
	delete [] particle;
}

State_Hyper_Vector& State_Hyper_Vector::operator= ( const State_Hyper_Vector& sv)
{
	gsl_rng_memcpy (gsl_r, sv.gsl_r);
	for (int i = 0; i < N; i++)
		particle[i] = sv.particle[i];
	return *this;
}

State_Hyper_Vector& State_Hyper_Vector::operator+ (const State_Hyper_Vector& s1) const
{
	for (int i = 0; i < N; i++)
	{
		particle[i].r = (particle[i].r + s1.particle[i].r);
		particle[i].r.Periodic_Transform();
		particle[i].theta = particle[i].theta + s1.particle[i].theta;
		particle[i].theta -= (int) (result.particle[i].theta / (2*M_PI));
	}
	return *this;
}


State_Hyper_Vector State_Hyper_Vector::operator- (const State_Hyper_Vector& s1) const
{
	for (int i = 0; i < N; i++)
	{
		particle[i].r = (particle[i].r - s1.particle[i].r);
		particle[i].r.Periodic_Transform();
		particle[i].theta = particle[i].theta - s1.particle[i].theta;
		particle[i].theta -= (int) (result.particle[i].theta / (2*M_PI));
	}
	return *this;
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

void State_Hyper_Vector::Rand(const Real position_amplitude, const Real angle_amplitude)
{
	for (int i = 0; i < N; i++)
	{
		particle[i].r.x = gsl_ran_flat(gsl_r, -position_amplitude, position_amplitude);
		particle[i].r.y = gsl_ran_flat(gsl_r, -position_amplitude, position_amplitude);
		particle[i].theta = gsl_ran_flat(gsl_r, -angle_amplitude, angle_amplitude);;
	}
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
