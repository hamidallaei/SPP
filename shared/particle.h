#ifndef _PARTICLE_
#define _PARTICLE_

#include "c2dvector.h"
#include "parameters.h"

class BasicParticle{
public:
	C2DVector r,v;
};

class BasicDynamicParticle: public BasicParticle {  // This is an abstract class for all of dynamic particles
public:
	C2DVector f;
	int neighbor_size;
	Real theta;
	static Real noise_amplitude;

	void Init();
	void Init(C2DVector);
	void Init(C2DVector, C2DVector);
	virtual void Reset();
	void Move();
	void Interact();
};

void BasicDynamicParticle::Init()
{
	r.Rand();
	v.Rand(1.0);
	theta = atan(v.y/v.x);
	if (v.x < 0)
		theta += PI;
	theta -= 2*PI * (int (theta / (2*PI)));
	v.x = cos(theta);
	v.y = sin(theta);
	Reset();
}

void BasicDynamicParticle::Init(C2DVector position)
{
	r = position;
	v.Rand(1.0);
	theta = atan(v.y/v.x);
	if (v.x < 0)
		theta += PI;
	theta -= 2*PI * (int (theta / (2*PI)));
	v.x = cos(theta);
	v.y = sin(theta);
	Reset();
}

void BasicDynamicParticle::Init(C2DVector position, C2DVector velocity)
{
	r = position;
	v = velocity;
	theta = atan(v.y/v.x);
	if (v.x < 0)
		theta += PI;
	theta -= 2*PI * (int (theta / (2*PI)));
	v.x = cos(theta);
	v.y = sin(theta);
	Reset();
}

void BasicDynamicParticle::Reset() {}

class VicsekParticle: public BasicDynamicParticle {
public:
	Real average_theta;
	void Move()
	{
		average_theta /= neighbor_size;
		theta = theta + average_theta + gsl_ran_gaussian(C2DVector::gsl_r,noise_amplitude);
		C2DVector old_v = v;
		v.x = cos(theta);
		v.y = sin(theta);
		r += v*dt;
		#ifdef PERIODIC_BOUNDARY_CONDITION
			r.Periodic_Transform();
		#endif
		Reset();
	}

	void Reset()
	{
		neighbor_size = 1;
		average_theta = 0;
		f.Null();
	}

	void Interact(VicsekParticle& p)
	{
		C2DVector dr = r - p.r;
		#ifdef PERIODIC_BOUNDARY_CONDITION
			dr.Periodic_Transform();
		#endif
		Real d2 = dr.Square();

		if (d2 < 1)
		{
			neighbor_size++;
			p.neighbor_size++;
			Real dtheta = p.theta - theta;
			dtheta -= 2*PI * (int (dtheta / PI));
			average_theta += dtheta;
			p.average_theta -= dtheta;
		}
	}
};


class ContinuousParticle: public BasicDynamicParticle {
public:
	Real torque;
	static Real g;
	static Real alpha;

	ContinuousParticle();
	void Move()
	{
		#ifdef COMPARE
			torque = round(digits*torque)/digits;
		#endif
		torque = g*torque + gsl_ran_gaussian(C2DVector::gsl_r,noise_amplitude);
		theta += torque*dt;
//		theta -= 2*PI * ((int) (theta / (PI)));
		C2DVector old_v = v;
		v.x = cos(theta);
		v.y = sin(theta);

		r += v*dt;
		#ifdef PERIODIC_BOUNDARY_CONDITION
			r.Periodic_Transform();
		#endif
		#ifdef TRACK_PARTICLE
			if (this == track_p && flag)
			{
				cout << "Particle:         " << setprecision(50)  << r << "\t" << theta << "\t" << torque << endl << flush;
			}
		#endif
		Reset();
	}

	virtual void Reset();
	void Interact(ContinuousParticle& p)
	{
		C2DVector dr = r - p.r;
		#ifdef PERIODIC_BOUNDARY_CONDITION
			dr.Periodic_Transform();
		#endif
		Real d2 = dr.Square();

		Real torque_interaction;
		if (d2 < 1)
		{
			neighbor_size++;
			p.neighbor_size++;

			Real d = sqrt(d2);

			torque_interaction = (1-alpha)*sin(p.theta - theta)/(PI);

			torque += torque_interaction;
			p.torque -= torque_interaction;

			torque -= alpha*(dr.x*v.y - dr.y*v.x) /(PI*d);
			p.torque += alpha*(dr.x*p.v.y - dr.y*p.v.x) /(PI*d);

			#ifdef TRACK_PARTICLE
				if (this == track_p && flag)
				{
//					if (abs(torque) > 0.1)
//						cout << "Intthis:     " << setprecision(100) << d2 << "\t" << torque_interaction << endl << flush;
//						cout << "Intthis:     " << setprecision(100) << r << "\t" << d2 << endl << flush;
				}
			#endif

			#ifdef TRACK_PARTICLE
				if (&p == track_p && flag)
				{
//						cout << "Intthat:     " << setprecision(100) <<  d2 << "\t" << torque_interaction << endl << flush;
//					if (abs(p.torque) > 0.1)
//						cout << "Intthat:     " << setprecision(100) << p.r << "\t" << d2 << "\t" << p.theta << "\t" << torque_interaction << endl << flush;
				}
			#endif
		}
	}
};

ContinuousParticle::ContinuousParticle()
{
	Init();
}

void ContinuousParticle::Reset()
{
	neighbor_size = 1;
	torque = 0;
	f.Null();
}

Real ContinuousParticle::g = 4;
Real ContinuousParticle::alpha = 1;


class RepulsiveParticle: public BasicDynamicParticle {
public:
	Real torque;
	static Real g;
	static Real kesi;

	RepulsiveParticle();
	void Move()
	{
		#ifdef COMPARE
			torque = round(digits*torque)/digits;
		#endif
		torque = torque + gsl_ran_gaussian(C2DVector::gsl_r,noise_amplitude);
		theta += torque*dt;
		C2DVector old_v = v;
		#ifdef COMPARE
			theta = round(digits*theta)/digits;
		#endif
		v.x = cos(theta);
		v.y = sin(theta);
		#ifdef COMPARE
			v.x = round(digits*v.x)/digits;
			v.y = round(digits*v.y)/digits;
			f.x = round(digits*f.x)/digits;
			f.y = round(digits*f.y)/digits;
		#endif
		v += f;
		r += v*dt;
		#ifdef PERIODIC_BOUNDARY_CONDITION
			r.Periodic_Transform();
		#endif
		#ifdef TRACK_PARTICLE
			if (this == track_p && flag)
			{
				cout << "Particle:         " << setprecision(50)  << v.y << endl << flush;
			}
		#endif
		Reset();
	}

	virtual void Reset();

	void Interact(RepulsiveParticle& p)
	{
		C2DVector dr = r - p.r;
		#ifdef PERIODIC_BOUNDARY_CONDITION
			dr.Periodic_Transform();
		#endif
		Real d2 = dr.Square();
		Real d = sqrt(d2);

		C2DVector interaction_force;
		if (d < r_c_p)
		{
			dr /= d; 
			Real r_c_p2 = r_c_p*r_c_p; 
			interaction_force = dr * A_p * ( exp(- d / sigma_p ) * ( 1. / d2 + 1. / (sigma_p * d)) - exp(- r_c_p / sigma_p ) * ( 1. / r_c_p2 + 1. / (sigma_p * r_c_p)) );

			f += interaction_force;
			p.f -= interaction_force;
		}

		Real torque_interaction;
		if (d < r_f_p)
		{
			neighbor_size++;
			p.neighbor_size++;

			torque_interaction = g*sin(p.theta - theta)/(PI);

			torque += torque_interaction;
			p.torque -= torque_interaction;

			#ifdef TRACK_PARTICLE
				if (this == track_p && flag)
				{
//					if (abs(torque) > 0.1)
//						cout << "Intthis:     " << setprecision(100) << d2 << "\t" << torque_interaction << endl << flush;
//						cout << "Intthis:     " << setprecision(100) << r << "\t" << d2 << endl << flush;
				}
			#endif

			#ifdef TRACK_PARTICLE
				if (&p == track_p && flag)
				{
//						cout << "Intthat:     " << setprecision(100) <<  d2 << "\t" << torque_interaction << endl << flush;
//					if (abs(p.torque) > 0.1)
//						cout << "Intthat:     " << setprecision(100) << p.r << "\t" << d2 << "\t" << p.theta << "\t" << torque_interaction << endl << flush;
				}
			#endif
		}
	}
};

RepulsiveParticle::RepulsiveParticle()
{
	Init();
}

void RepulsiveParticle::Reset()
{
	neighbor_size = 1;
	torque = 0;
	f.Null();
}

Real RepulsiveParticle::g = .5;
Real RepulsiveParticle::kesi = .5;

Real BasicDynamicParticle::noise_amplitude = .1;


#endif
