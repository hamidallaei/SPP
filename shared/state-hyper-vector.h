#ifndef _STATE_HYPER_VECTOR_
#define _STATE_HYPER_VECTOR_

#include "c2dvector.h"

class State_Hyper_Vector{
	gsl_rng* gsl_r;
	void Init_Random_Generator(int);
public:
	int N;
	BasicParticle0* particle;
	
	State_Hyper_Vector(int, int);
	State_Hyper_Vector(const State_Hyper_Vector&);
	~State_Hyper_Vector();

	State_Hyper_Vector& operator= ( const State_Hyper_Vector& sv);
	const State_Hyper_Vector operator+ (const State_Hyper_Vector& s1);
	const State_Hyper_Vector operator- (const State_Hyper_Vector& s1);
	State_Hyper_Vector& operator+= (const State_Hyper_Vector& s1);
	State_Hyper_Vector& operator-= (const State_Hyper_Vector& s1);
	const Real operator* (const State_Hyper_Vector& s1) const;

	void Set_C2DVector_Rand_Generator() const;
	void Get_C2DVector_Rand_Generator();

	void Rand(const Real position_amplitude, const Real angle_amplitude);
	Real Square() const;
	Real Magnitude() const;
};

void State_Hyper_Vector::Init_Random_Generator(int seed)
{
	const gsl_rng_type * T;
	gsl_rng_env_setup();
	gsl_rng_default_seed = seed;
	T = gsl_rng_default;
	gsl_r = gsl_rng_alloc (T);
	gsl_rng_memcpy (gsl_r, C2DVector::gsl_r);
}

State_Hyper_Vector::State_Hyper_Vector(int particle_number, int seed = 0) : N(particle_number)
{
	Init_Random_Generator(seed);
	particle = new BasicParticle0[N];
}

State_Hyper_Vector::State_Hyper_Vector(const State_Hyper_Vector& sv) : N(sv.N)
{
	Init_Random_Generator(0);
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

// Alwayse retun random generator of the lhs
const State_Hyper_Vector State_Hyper_Vector::operator+ (const State_Hyper_Vector& s1)
{
	State_Hyper_Vector result(*this);
	for (int i = 0; i < N; i++)
	{
		result.particle[i].r = (particle[i].r + s1.particle[i].r);
		result.particle[i].r.Periodic_Transform();
		result.particle[i].theta = particle[i].theta + s1.particle[i].theta;
		result.particle[i].theta -= 2*M_PI*floor(result.particle[i].theta / (2*M_PI));
	}
	return result;
}

// Alwayse retun random generator of the lhs
const State_Hyper_Vector State_Hyper_Vector::operator- (const State_Hyper_Vector& s1)
{
	State_Hyper_Vector result(*this);
	for (int i = 0; i < N; i++)
	{
		result.particle[i].r = (particle[i].r - s1.particle[i].r);
		result.particle[i].r.Periodic_Transform();
		result.particle[i].theta = particle[i].theta - s1.particle[i].theta;
		result.particle[i].theta -= 2*M_PI*floor(particle[i].theta / (2*M_PI));
	}
	return result;
}

// Alwayse retun random generator of the lhs
State_Hyper_Vector& State_Hyper_Vector::operator+= (const State_Hyper_Vector& s1)
{
	for (int i = 0; i < N; i++)
	{
		particle[i].r = (particle[i].r + s1.particle[i].r);
		particle[i].r.Periodic_Transform();
		particle[i].theta = particle[i].theta + s1.particle[i].theta;
		particle[i].theta -= 2*M_PI*floor(particle[i].theta / (2*M_PI));
	}
	return *this;
}

// Alwayse retun random generator of the lhs
State_Hyper_Vector& State_Hyper_Vector::operator-= (const State_Hyper_Vector& s1)
{
	for (int i = 0; i < N; i++)
	{
		particle[i].r = (particle[i].r - s1.particle[i].r);
		particle[i].r.Periodic_Transform();
		particle[i].theta = particle[i].theta - s1.particle[i].theta;
		particle[i].theta -= 2*M_PI*floor(particle[i].theta / (2*M_PI));
	}
	return *this;
}

const Real State_Hyper_Vector::operator* (const State_Hyper_Vector& s1) const
{
	Real result = 0;

	for (int i = 0; i < N; i++)
	{
		result += particle[i].r * s1.particle[i].r;
		result += cos(particle[i].theta - s1.particle[i].theta);
	}
	return result;
}

void State_Hyper_Vector::Set_C2DVector_Rand_Generator() const
{
	gsl_rng_memcpy (C2DVector::gsl_r, gsl_r);
}

void State_Hyper_Vector::Get_C2DVector_Rand_Generator()
{
	gsl_rng_memcpy (gsl_r, C2DVector::gsl_r);
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
//		result += particle[i].r * particle[i].r;
		Real dtheta = particle[i].theta - 2*M_PI*ceil((particle[i].theta - M_PI) / (2*M_PI));
		result += dtheta*dtheta;
	}
	return result;
}

Real State_Hyper_Vector::Magnitude() const
{
	return sqrt(Square());
}

#endif
