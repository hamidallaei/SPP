#ifndef _PARTICLE_
#define _PARTICLE_

#include "c2dvector.h"
#include "parameters.h"
#include "force-fields.h"
#include <vector>

//###################################################################
class BasicParticle00{
public:
	C2DVector r;
	virtual void Write(std::ostream&);
};

void BasicParticle00::Write(std::ostream& os)
{
	r.write(os);
}
//###################################################################

//###################################################################
class BasicParticle0{
public:
	C2DVector r;
	Real theta;
	virtual void Write(std::ostream&);
};

void BasicParticle0::Write(std::ostream& os)
{
	r.write(os);
	float temp_float = (float) theta;
	os.write((char*) &temp_float,sizeof(float) / sizeof(char));
}
//###################################################################

//###################################################################
class BasicParticle{
public:
	C2DVector r,v;
	#ifdef NonPeriodicCompute
		C2DVector r_original;
	#endif
	void Write(std::ostream&);
};
//###################################################################

//###################################################################
class BasicDynamicParticle: public BasicParticle {  // This is an abstract class for all of dynamic particles
public:
	int neighbor_size;
	Real theta;
	static Real Dr; // Dr is for continuous time particles. It is not usefull in discreate vicsek particles.
	static Real noise_amplitude;
	static Real speed;
	vector<int> neighbor_id; // id of neighboring particles

	void Init();
	void Init(C2DVector);
	void Init(C2DVector, C2DVector);
	virtual void Reset();

	static void Set_Dr(const Real input_Dr);
	void Set_speed(const Real input_speed);
	void Set_noise_amplitude(const Real input_speed);

	void Move();
	void Interact();
	void Write(std::ostream&);
};

void BasicDynamicParticle::Init()
{
	r.Rand();
	v.Rand(1.0);
	theta = atan(v.y/v.x);
	if (v.x < 0)
		theta += PI;
	theta -= 2*PI * (int (theta / (2*PI)));v.x = cos(theta);
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

void BasicDynamicParticle::Set_Dr(const Real input_Dr)
{
	Dr = input_Dr;
	noise_amplitude = sqrt(2*Dr*dt);
}

void BasicDynamicParticle::Set_noise_amplitude(const Real input_noise_amplitude)
{
	noise_amplitude = input_noise_amplitude;
	Dr = 0.5*noise_amplitude*noise_amplitude / dt;
}

void BasicDynamicParticle::Set_speed(const Real input_speed)
{
	speed = input_speed;
}

void BasicDynamicParticle::Write(std::ostream& os)
{
	r.write(os);
	C2DVector v_temp;
	v_temp.x = cos(theta);
	v_temp.y = sin(theta);
	v_temp.write(os);
}

Real BasicDynamicParticle::noise_amplitude = 0.0;
Real BasicDynamicParticle::Dr = 0.0;
Real BasicDynamicParticle::speed = 1;
//###################################################################


//###################################################################
class VicsekParticle: public BasicDynamicParticle {
public:
	Real average_theta;
	void Move()
	{
		if (neighbor_size > 0)
			average_theta /= neighbor_size;
		else
			average_theta = theta;
		theta = average_theta + noise_amplitude*gsl_ran_flat(C2DVector::gsl_r,-M_PI,M_PI);
		C2DVector old_v = v;
		v.x = cos(theta);
		v.y = sin(theta);
		r += v*(dt*speed);
		#ifdef PERIODIC_BOUNDARY_CONDITION
			r.Periodic_Transform();
		#endif
		Reset();
	}

	void Reset()
	{
		neighbor_size = 0;
		average_theta = 0;
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
			average_theta += p.theta;
			p.average_theta += theta;
		}
	}
};
//###################################################################

//###################################################################
class VicsekParticle2: public BasicDynamicParticle {
public:
	C2DVector average_v;
	static Real beta;
	static Real rc;
	void Move()
	{
		if (neighbor_size == 0)
			average_v = v;
		Real rand_angle = gsl_ran_flat(C2DVector::gsl_r,-M_PI,M_PI);
		average_v.x += neighbor_size*noise_amplitude*cos(rand_angle);
		average_v.y += neighbor_size*noise_amplitude*sin(rand_angle);
		theta = atan2(average_v.y,average_v.x);
		v.x = cos(theta);
		v.y = sin(theta);
		r += v*(dt*speed);
		#ifdef PERIODIC_BOUNDARY_CONDITION
			r.Periodic_Transform();
		#endif
		Reset();
	}

	void Reset()
	{
		neighbor_size = 0;
		average_v.Null();
	}

	void Interact(VicsekParticle2& p)
	{
		C2DVector dr = r - p.r;
		#ifdef PERIODIC_BOUNDARY_CONDITION
			dr.Periodic_Transform();
		#endif
		Real d2 = dr.Square();
		Real d = sqrt(d2);
		C2DVector ehat = dr / d;

		if (d2 < 1)
		{
			neighbor_size++;
			p.neighbor_size++;
			Real dtheta = p.theta - theta;
			dtheta -= 2*PI * (int (dtheta / PI));
			average_v += p.v;
			p.average_v += v;
			if (d < rc)
			{
				Real force_amplitude = beta/(1 + exp(d/rc-2));
				average_v += ehat*force_amplitude;
				p.average_v -= ehat*force_amplitude;
			}
		}
	}
};

Real VicsekParticle2::beta = 2.5;
Real VicsekParticle2::rc = 0.127;
//###################################################################

//###################################################################
class ContinuousParticle: public BasicDynamicParticle {
public:
	Real torque;
	static Real g, gw;
	static Real alpha;
	static Real Dr,K,Kamp;

	ContinuousParticle();

	static void Set_g(const Real);
	static void Set_gw(const Real);
	static void Set_alpha(const Real);
	static void Set_K(const Real);

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

		r += v*(speed*dt);
		r.x += gsl_ran_gaussian(C2DVector::gsl_r,Kamp);
		r.y += gsl_ran_gaussian(C2DVector::gsl_r,Kamp);
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
}

void ContinuousParticle::Set_g(const Real input_g)
{
	g = input_g;
}

void ContinuousParticle::Set_gw(const Real input_gw)
{
	gw = input_gw;
}

void ContinuousParticle::Set_alpha(const Real input_alpha)
{
	alpha = input_alpha;
}

void ContinuousParticle::Set_K(const Real input_K)
{
	K = input_K;
	Kamp = sqrt(2*K*dt);
}

Real ContinuousParticle::g = 4;
Real ContinuousParticle::alpha = 0.5;
Real ContinuousParticle::Dr = 0;
Real ContinuousParticle::K = 0;
Real ContinuousParticle::Kamp = 0;
Real ContinuousParticle::gw = 20;
//###################################################################



//###################################################################
class MarkusParticle: public BasicDynamicParticle {
public:
	Real torque;
	static Real mu_plus;
	static Real mu_minus;
	static Real kapa;
	static Real D_phi;
	static Real kisi_r, kisi_a, kisi;

	MarkusParticle();
	void Move()
	{
		#ifdef COMPARE
			torque = round(digits*torque)/digits;
		#endif
		torque = torque + gsl_ran_gaussian(C2DVector::gsl_r,noise_amplitude);
		theta += torque*dt;
//		theta -= 2*PI * ((int) (theta / (PI)));
		C2DVector old_v = v;
		v.x = cos(theta);
		v.y = sin(theta);

		r += v*(dt*speed);
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
	void Interact(MarkusParticle& p)
	{
		C2DVector dr = r - p.r;
		#ifdef PERIODIC_BOUNDARY_CONDITION
			dr.Periodic_Transform();
		#endif
		Real d2 = dr.Square();
		Real d = sqrt(d2);

		Real torque_interaction;
		if (d < kisi)
		{
			if (d < kisi_a)
			{
				neighbor_size++;
				p.neighbor_size++;
				torque_interaction = mu_plus*(1-(d2/(kisi_a*kisi_a)))*sin(p.theta - theta);
				torque += torque_interaction;
				p.torque -= torque_interaction;
			}
			else
			{
				neighbor_size++;
				p.neighbor_size++;
				torque_interaction = mu_minus*4*(d - kisi_a)*(1-d)*sin(theta - p.theta) / ((1-kisi_a)*(1-kisi_a));
				torque += torque_interaction;
				p.torque -= torque_interaction;
			}

			if (d < kisi_r)
			{
				Real alpha = atan2(-dr.y,-dr.x);
				Real factor = (1.0 - d / kisi_r);
				torque += factor*kapa*sin(theta - alpha);
				p.torque -= factor*kapa*sin(p.theta - alpha);
			}

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

MarkusParticle::MarkusParticle()
{
	Init();
}

void MarkusParticle::Reset()
{
	neighbor_size = 1;
	torque = 0;
}

Real MarkusParticle::mu_plus = 1;
Real MarkusParticle::mu_minus = 1;
Real MarkusParticle::kapa = 40;//40.0;
Real MarkusParticle::D_phi = 1;
Real MarkusParticle::kisi_r = 0.1;
Real MarkusParticle::kisi_a = 1;//0.2;
Real MarkusParticle::kisi = 1;
//###################################################################

//###################################################################
class RepulsiveParticle: public BasicDynamicParticle {
public:
	C2DVector f, r_old; // force, old position
	Real torque; // tourque acting on the chain
	Real theta_old; // self propullsion angle
	static Real F0; // self propullsion strength
	Real dtheta; // amount of noise added to theta
	static Real sigma_p;
	static Real repulsion_radius;
	static Real alignment_radius;		// flocking radius with particles
	static Real A_p;
	static Real g;
	static Real r_c_w;
	static Real r_f_w;
	static Real sigma_w;
	static Real A_w;
	static Real g_w;
	static int nb;
	#ifdef RUNGE_KUTTA4
		C2DVector k1_f,k2_f,k3_f,k4_f;
		Real k1_torque, k2_torque, k3_torque, k4_torque;
	#endif

	RepulsiveParticle();

	static void Set_F0(const Real);
	static void Set_sigma_p(const Real);
	static void Set_repulsion_radius(const Real);
	static void Set_alignment_radius(const Real);
	static void Set_A_p(const Real);
	static void Set_g(const Real);
	static void Set_nb(const int);

	void Reset();
	void Move();
	void Move_Runge_Kutta2_1();
	void Move_Runge_Kutta2_2();
	void Move_Runge_Kutta4_1();
	void Move_Runge_Kutta4_2();
	void Move_Runge_Kutta4_3();
	void Move_Runge_Kutta4_4();
	virtual void Noise_Gen();
	void Interact(RepulsiveParticle& ac);
	void Write(std::ostream& os);
};

RepulsiveParticle::RepulsiveParticle()
{
	Init();
}

void RepulsiveParticle::Set_F0(const Real input_F0) {F0 = input_F0;}
void RepulsiveParticle::Set_sigma_p(const Real input_sigma_p)  {sigma_p = input_sigma_p;}
void RepulsiveParticle::Set_repulsion_radius(const Real input_repulsion_radius) {repulsion_radius = input_repulsion_radius;}
void RepulsiveParticle::Set_alignment_radius(const Real input_alignment_radius) {alignment_radius = input_alignment_radius;}
void RepulsiveParticle::Set_A_p(const Real input_A_p) {A_p = input_A_p;}
void RepulsiveParticle::Set_g(const Real input_g) {g = input_g;}
void RepulsiveParticle::Set_nb(const int input_nb) {nb = input_nb;}

void RepulsiveParticle::Reset()
{
	neighbor_size = 1;
	torque = 0;
	f = v;
}

inline void RepulsiveParticle::Noise_Gen()
{
	dtheta =  gsl_ran_gaussian(C2DVector::gsl_r,noise_amplitude);
}

inline void RepulsiveParticle::Move()
{
	Noise_Gen();
	theta += dt*(torque);
	theta += dtheta;  // add noise
	v.x = cos(theta); // v is the direction of the particle
	v.y = sin(theta); // v is the direction of the particle

	r_original += f*dt;
	r = r_original;

	#ifdef PERIODIC_BOUNDARY_CONDITION
		r.Periodic_Transform();
	#endif

	Reset();
}

void RepulsiveParticle::Move_Runge_Kutta2_1() // half step forward
{
	Noise_Gen();

	theta_old = theta;
	r_old = r_original;

	theta += half_dt*torque;
	theta += dtheta/2;

	theta = theta - floor(theta/(2*M_PI) + 0.5)*2*M_PI; // -PI < theta < PI

	v.x = cos(theta); // v is the direction of the particle
	v.y = sin(theta); // v is the direction of the particle

	r += f*half_dt;

	#ifdef PERIODIC_BOUNDARY_CONDITION
		r.Periodic_Transform();
	#endif

	Reset();
}

void RepulsiveParticle::Move_Runge_Kutta2_2() // one step forward
{
	theta = theta_old + dt*torque;
	theta += dtheta;

	theta = theta - floor(theta/(2*M_PI) + 0.5)*2*M_PI;

	v.x = cos(theta); // v is the direction of the particle
	v.y = sin(theta); // v is the direction of the particle

	r_original = r_old + f*dt;
	r = r_original;

	#ifdef PERIODIC_BOUNDARY_CONDITION
		r.Periodic_Transform();
	#endif

	Reset();
}

#ifdef RUNGE_KUTTA4
void RepulsiveParticle::Move_Runge_Kutta4_1() // half step forward
{
	Noise_Gen();

	theta_old = theta;
	r_old = r_original;

	theta += half_dt*torque;
	theta += dtheta/2;

	theta = theta - floor(theta/(2*M_PI) + 0.5)*2*M_PI; // -PI < theta < PI

	v.x = cos(theta); // v is the direction of the particle
	v.y = sin(theta); // v is the direction of the particle

	r += f*half_dt;

	#ifdef PERIODIC_BOUNDARY_CONDITION
		r.Periodic_Transform();
	#endif

	k1_f = f;
	k1_torque = torque;

	Reset();
}

void RepulsiveParticle::Move_Runge_Kutta4_2() // half step forward correction
{
	theta = theta_old + half_dt*torque;
	theta += dtheta/2;

	theta = theta - floor(theta/(2*M_PI) + 0.5)*2*M_PI;

	v.x = cos(theta); // v is the direction of the particle
	v.y = sin(theta); // v is the direction of the particle

	r_original = r_old + f*half_dt;
	r = r_original;

	#ifdef PERIODIC_BOUNDARY_CONDITION
		r.Periodic_Transform();
	#endif

	k2_f = f;
	k2_torque = torque;

	Reset();
}

void RepulsiveParticle::Move_Runge_Kutta4_3() // full step forward
{
	theta = theta_old + dt*torque;
	theta += dtheta;

	theta = theta - floor(theta/(2*M_PI) + 0.5)*2*M_PI;

	v.x = cos(theta); // v is the direction of the particle
	v.y = sin(theta); // v is the direction of the particle

	r_original = r_old + f*dt;
	r = r_original;

	#ifdef PERIODIC_BOUNDARY_CONDITION
		r.Periodic_Transform();
	#endif

	k3_f = f;
	k3_torque = torque;

	Reset();
}

void RepulsiveParticle::Move_Runge_Kutta4_4() // full step forward corrected
{
	theta = theta_old + (k1_torque + k2_torque*2 + k3_torque*2 + torque)*dt_over_6;
	theta += dtheta;

	theta = theta - floor(theta/(2*M_PI) + 0.5)*2*M_PI;

	v.x = cos(theta); // v is the direction of the particle
	v.y = sin(theta); // v is the direction of the particle

	r_original = r_old + (k1_f + k2_f*2 + k3_f*2 + f)*dt_over_6;
	r = r_original;

	#ifdef PERIODIC_BOUNDARY_CONDITION
		r.Periodic_Transform();
	#endif

	Reset();
}
#endif

void RepulsiveParticle::Interact(RepulsiveParticle& p)
{
	C2DVector dr = r - p.r;
	#ifdef PERIODIC_BOUNDARY_CONDITION
		dr.Periodic_Transform();
	#endif
	Real d2 = dr.Square();
	Real d = sqrt(d2);

	C2DVector interaction_force;
	if (d < repulsion_radius)
	{
		interaction_force = R12_Repulsive_Truncated(dr,d,repulsion_radius,sigma_p,A_p);

		f += interaction_force;
		p.f -= interaction_force;
	}

	Real torque_interaction;
	if (d < alignment_radius)
	{
		neighbor_size++;
		p.neighbor_size++;

		torque_interaction = g*sin(p.theta - theta) / (PI*alignment_radius*alignment_radius);

		torque += torque_interaction;
		p.torque -= torque_interaction;
	}
}

void RepulsiveParticle::Write(std::ostream& os)
{
	r_original.write(os);
	float temp_float = (float) theta;
	os.write((char*) &temp_float,sizeof(float) / sizeof(char));
}

Real RepulsiveParticle::F0 = 1.0;
Real RepulsiveParticle::g = 1.0;
// Interactions for repulsive particles
Real RepulsiveParticle::A_p = 2.;		// interaction strength
Real RepulsiveParticle::sigma_p = 1.0;		// sigma in Potential
Real RepulsiveParticle::alignment_radius = 1.1;		// flocking radius with particles
Real RepulsiveParticle::repulsion_radius = 1.2;	// repulsive cutoff radius with particles
Real RepulsiveParticle::r_c_w;
Real RepulsiveParticle::r_f_w;
Real RepulsiveParticle::A_w;
Real RepulsiveParticle::sigma_w;
Real RepulsiveParticle::g_w;
int RepulsiveParticle::nb = 1;

//###################################################################


//###################################################################
class ActiveBrownianChain : public BasicDynamicParticle {
public:
//	const int dof = 3; // degree of freedom. It is required in the comminucation between nodes. The dof coordinates are sent and received.
	int nb; // number of beads
	Real m_parallel;  // Mobility parallel to the direction
	Real m_perpendicular; // Mobility perpendicular to the direction
	Real k_perpendicular; // Mobility perpendicular to the direction

	C2DVector f, r_old; // force, old position
	Real torque; // tourque acting on the chain
	Real theta; // self propullsion angle
	Real theta_old; // self propullsion angle
	Real F0; // self propullsion strength
	Real dtheta; // amount of noise added to theta
	static Real sigma_p;
	static Real repulsion_radius;
	static Real A_p;
	static Real torque0; // The intrinsinc torque in the particle. This make the motion chiral
	static Real R0; // The intrinsinc radius of motion. This make the motion chiral

	ActiveBrownianChain();

	void Set_F0(const Real);
	static void Set_sigma_p(const Real);
	static void Set_repulsion_radius(const Real);
	static void Set_A_p(const Real);
	void Set_nb(const int);
	void Set_R0(const Real); // Sould be called after Set_nb
	static void Set_torque0(const Real);

	void Reset();
	void Set_Parameters(int input_nb, Real input_F0);
	void Move();
	void Move_Runge_Kutta2_1();
	void Move_Runge_Kutta2_2();
	virtual void Noise_Gen();
	void Interact(ActiveBrownianChain& ac);
	void Write(std::ostream& os);
};

ActiveBrownianChain::ActiveBrownianChain()
{
	Init();
	F0 = 0;
}

void ActiveBrownianChain::Set_F0(const Real input_F0) {F0 = input_F0;}
void ActiveBrownianChain::Set_sigma_p(const Real input_sigma_p)  {sigma_p = input_sigma_p;}
void ActiveBrownianChain::Set_repulsion_radius(const Real input_repulsion_radius) {repulsion_radius = input_repulsion_radius;}
void ActiveBrownianChain::Set_A_p(const Real input_A_p) {A_p = input_A_p;}
void ActiveBrownianChain::Set_torque0(const Real input_torque0) {torque0 = input_torque0;}

void ActiveBrownianChain::Set_R0(const Real input_R0) // Should be called after set nb
{
	R0 = input_R0;
	torque0 = m_parallel*F0 / (k_perpendicular*R0);
}

void ActiveBrownianChain::Set_nb(const int input_nb)
{
	nb = input_nb;
	if (nb == 1)
	{	// case of a single spherial particle
		m_parallel = 1;  	 // Mobility parallel to the direction
		m_perpendicular = 1; // Mobility perpendicular to the direction
		k_perpendicular = 1; // Mobility perpendicular to the direction
	}
	else
	{	// these values are for two-sphere swimmers. TODO one has to import general relations.
		m_parallel = 1.0;
		m_perpendicular = 0.87;
		k_perpendicular = 4.8;
	}
}

void ActiveBrownianChain::Reset()
{
	neighbor_size = 0;
	torque = torque0;
	f = v*F0; // F0 is self-propullsion force, v is the direction of the particle
}

void ActiveBrownianChain::Set_Parameters(int input_nb, Real input_F0)
{
	Set_nb(input_nb);
	Set_F0(input_F0);
	Set_R0(R0);
	Set_Dr(Dr);
}

inline void ActiveBrownianChain::Noise_Gen()
{
	dtheta =  gsl_ran_gaussian(C2DVector::gsl_r,noise_amplitude);
}

inline void ActiveBrownianChain::Move()
{
	Noise_Gen();
	theta += k_perpendicular*dt*(torque);
	theta += dtheta;  // add noise
	v.x = cos(theta); // v is the direction of the particle
	v.y = sin(theta); // v is the direction of the particle

	r_original += f*(dt*m_perpendicular);
	r_original += v*((f*v)*(dt*(m_parallel - m_perpendicular)));
	r = r_original;

	#ifdef PERIODIC_BOUNDARY_CONDITION
		r.Periodic_Transform();
	#endif

	Reset();
}

void ActiveBrownianChain::Move_Runge_Kutta2_1() // half step forward
{
	Noise_Gen();

	theta_old = theta;
	r_old = r_original;

	theta += k_perpendicular*half_dt*(torque);
	theta += dtheta/2;

	theta = theta - floor(theta/(2*M_PI) + 0.5)*2*M_PI; // -PI < theta < PI

	r += f*(half_dt*m_perpendicular);
	r += v*((f*v)*(half_dt*(m_parallel - m_perpendicular)));

	v.x = cos(theta); // v is the direction of the particle
	v.y = sin(theta); // v is the direction of the particle

	#ifdef PERIODIC_BOUNDARY_CONDITION
		r.Periodic_Transform();
	#endif

	Reset();
}

void ActiveBrownianChain::Move_Runge_Kutta2_2() // one step forward
{
	theta = theta_old + k_perpendicular*dt*(torque);
	theta += dtheta;

	theta = theta - floor(theta/(2*M_PI) + 0.5)*2*M_PI;

	r_original = r_old + f*(dt*m_perpendicular);
	r_original += v*((f*v)*(dt*(m_parallel - m_perpendicular)));

	r = r_original;

	v.x = cos(theta); // v is the direction of the particle
	v.y = sin(theta); // v is the direction of the particle

	#ifdef PERIODIC_BOUNDARY_CONDITION
		r.Periodic_Transform();
	#endif

	Reset();
}

inline void ActiveBrownianChain::Interact(ActiveBrownianChain& ac)
{
	C2DVector dr0 = r - ac.r;
	#ifdef PERIODIC_BOUNDARY_CONDITION
		dr0.Periodic_Transform();
	#endif
	C2DVector dr;
	for (int i = 0; i < nb; i++)
		for (int j = 0; j < ac.nb; j++)
		{
			C2DVector s_this, s_that;	// s is position of the beads on each swimmer wrt swimmer middle point
			s_this = v*(sigma_p*((1-nb)/2.0+i)); 	   //  v is the direction of the particle
			s_that = ac.v*(sigma_p*((1-ac.nb)/2.0+j)); //  v is the direction of the particle

			dr = dr0 + s_this - s_that;

			Real d2 = dr.Square();
			Real d = sqrt(d2);

			C2DVector interaction_force;
			if (d < repulsion_radius)
			{
				interaction_force = R12_Repulsive_Truncated(dr,d,repulsion_radius,sigma_p,A_p);

				f += interaction_force;
				ac.f -= interaction_force;

				torque += s_this.x*interaction_force.y - s_this.y*interaction_force.x;
				ac.torque += -s_that.x*interaction_force.y + s_that.y*interaction_force.x;
			}
		}
}



void ActiveBrownianChain::Write(std::ostream& os)
{
/*	C2DVector v_temp, r_temp;*/
/*	v_temp.x = cos(theta);*/
/*	v_temp.y = sin(theta);*/

/*	for (int i = 0; i < nb; i++)*/
/*	{*/
/*		C2DVector s_i = v_temp*(sigma_p*((1-nb)/2.0 + i));*/

/*		r_temp = r + s_i;*/

/*		r_temp.write(os);*/
/*		v_temp.write(os);*/
/*	}*/

	if (F0 == 0)
		r_original.write(os);
	else
	{
		r_original.write(os);
		float temp_float = (float) theta;
		os.write((char*) &temp_float,sizeof(float) / sizeof(char));
	}
}

// Interactions parameters for active chain
Real ActiveBrownianChain::A_p = 1.;		// interaction strength
Real ActiveBrownianChain::sigma_p = 1.0;		// sigma in Yukawa Potential
Real ActiveBrownianChain::repulsion_radius = 1.1;		// repulsive cutoff radius
Real ActiveBrownianChain::torque0 = 0; // The intrinsinc torque in the particle. This make the motion chiral
Real ActiveBrownianChain::R0 = 0; // The intrinsinc radius of motion. This make the motion chiral
//###################################################################


//###################################################################
class RTPChain : public ActiveBrownianChain {
public:
	int tumble_flag; // in tumbling state = -+1, otherwise 0
	Real tumble_elapesed_time; // The time that the particle stays in tumbling state.
	static Real lambda; // tumbling rate
	static Real t_tumble; // tumbling duration
	static Real torque_tumble; // torque strength of a tumble

	RTPChain();

	void Set_lambda(const Real);
	void Set_t_tumble(const Real);
	void Set_torque_tumble(const Real);

	virtual void Noise_Gen();
};

RTPChain::RTPChain()
{
	Init();
	F0 = 0;
	tumble_flag = 0; // in tumbling state = true
	tumble_elapesed_time = 0; // The time that the particle stays in tumbling state.
}

void RTPChain::Set_lambda(const Real input_lambda) {lambda = input_lambda;}
void RTPChain::Set_t_tumble(const Real input_t_tumble) {t_tumble = input_t_tumble;}
void RTPChain::Set_torque_tumble(const Real input_torque_tumble) {torque_tumble = input_torque_tumble;}

inline void RTPChain::Noise_Gen()
{
	if (tumble_flag == 0)
	{
		Real random_number = gsl_rng_uniform(C2DVector::gsl_r);
		if (random_number < lambda*dt)
		{
			tumble_flag = 2*gsl_rng_uniform_int(C2DVector::gsl_r,2) - 1;
			tumble_elapesed_time += dt;
		}
	}
	else
	{
		if (tumble_elapesed_time < t_tumble)
			tumble_elapesed_time += dt;
		else
		{
			tumble_elapesed_time = 0;
			tumble_flag = 0;
		}
	}
	dtheta = k_perpendicular*dt*tumble_flag*torque_tumble;
}


Real RTPChain::lambda = 0.1; // tumbling rate
Real RTPChain::t_tumble = 10*RTPChain::lambda; // tumbling duration
Real RTPChain::torque_tumble = 1; // torque strength of a tumble

//###################################################################


//####################################################################
class EjtehadiParticle: public ContinuousParticle {
public:
	Real torque, torque_phi;
	Real phi;
	Real speed;
	static Real g, gw, g_phi;
	static Real alpha;
	static Real Dr,K,Kamp,D_phi;
	static Real vmin,vmax,vmid,vamp;
	static Real Rphi;
	static Real omega;
	static Real noise_amplitude_phi;

	EjtehadiParticle();
	void Move()
	{
		#ifdef COMPARE
			torque = round(digits*torque)/digits;
		#endif
		torque = g*torque + gsl_ran_gaussian(C2DVector::gsl_r,noise_amplitude);
		torque_phi = torque_phi + gsl_ran_gaussian(C2DVector::gsl_r,noise_amplitude_phi);
		theta += torque*dt;
		phi += torque_phi*dt;

		C2DVector old_v = v;
		v.x = cos(theta);
		v.y = sin(theta);

		speed = vmid + vamp*cos(phi);
		r += v*(speed*dt);
		r.x += gsl_ran_gaussian(C2DVector::gsl_r,Kamp);
		r.y += gsl_ran_gaussian(C2DVector::gsl_r,Kamp);
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
	void Interact(EjtehadiParticle& p)
	{
		C2DVector dr = r - p.r;
		#ifdef PERIODIC_BOUNDARY_CONDITION
			dr.Periodic_Transform();
		#endif
		Real d2 = dr.Square();
		Real d = sqrt(d2);

		Real torque_interaction;
		if (d2 < 1)
		{
			neighbor_size++;
			p.neighbor_size++;

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
		if (d < Rphi)
		{
			torque_interaction = g_phi*sin(p.phi - phi)/(PI);
			torque_phi += torque_interaction;
			p.torque_phi -= torque_interaction;
		}
	}
	static void Set_Variables();
};

EjtehadiParticle::EjtehadiParticle()
{
	Init();
}

void EjtehadiParticle::Reset()
{
	neighbor_size = 1;
	torque = 0;
	torque_phi = omega;
}

void EjtehadiParticle::Set_Variables()
{
	Kamp = sqrt(2*K*dt);
	noise_amplitude = sqrt(2*Dr/dt);
	noise_amplitude_phi = sqrt(2*D_phi/dt);
	vamp = (vmax - vmin) / 2.0;
	vmid = (vmax + vmin) / 2.0;
}

Real EjtehadiParticle::g = 1;
Real EjtehadiParticle::gw = 20;
Real EjtehadiParticle::g_phi = 1;
Real EjtehadiParticle::alpha = 0.5;
Real EjtehadiParticle::Rphi = 1;
Real EjtehadiParticle::D_phi = 0;
Real EjtehadiParticle::Dr = 0;
Real EjtehadiParticle::K = 0.2;
Real EjtehadiParticle::Kamp = 0.2;
Real EjtehadiParticle::noise_amplitude_phi = 0;
Real EjtehadiParticle::omega = 1;
Real EjtehadiParticle::vmin = 0.1;
Real EjtehadiParticle::vamp = 1.0;
Real EjtehadiParticle::vmax = 1.1;
Real EjtehadiParticle::vmid = 0.5;
//###################################################################

#endif
