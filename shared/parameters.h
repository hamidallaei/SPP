#ifndef _PARAMETERS_
#define _PARAMETERS_

#define RUNGE_KUTTA
#define PERIODIC_BOUNDARY_CONDITION
#define NonPeriodicCompute
//#define CIRCULAR_BOX
//#define verlet_list
// This is for checking particles outside of the box
#define DEBUG
// This tracks a specific particle. The id of the tracking particle is given below.
//#define TRACK_PARTICLE
// This will round torques to avoid any difference of this program and other versions caused by truncation of numbers (if we change order of a sum, the result will change because of the truncation error)
//#define COMPARE

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>
#include <ctime>
#include <algorithm>

using namespace std;

class BasicDynamicParticle;
class VicsekParticle;
class VicsekParticle2;
class ContinuousParticle;
class MarkusParticle;
class RepulsiveParticle;
class ActiveBrownianChain;
class RTPChain;
class EjtehadiParticle;

//#define ejtehadi

typedef double Real;

#ifdef ejtehadi
	typedef EjtehadiParticle Particle;
	const int dof = 4; // degree of freedom. It is required in the comminucation between nodes. The dof coordinates are sent and received.
#else
//	typedef ContinuousParticle Particle;
//	typedef VicsekParticle2 Particle;
//	typedef MarkusParticle Particle;
//	typedef RepulsiveParticle Particle;
	typedef ActiveBrownianChain Particle;
//	typedef RTPChain Particle;
	#ifdef NonPeriodicCompute
		const int dof = 5;
	#else
		const int dof = 3;
	#endif

#endif

long int seed = 10;

const Real PI = M_PI;

const int max_wall_num = 8;
const int max_N = 64000;

// Box
const int Lx_int = 64;
const int Ly_int = Lx_int;
const int L_int = Lx_int;
const Real Lx = Lx_int;
const Real Ly = Ly_int;
const Real Lx2 = 2*Lx;
const Real Ly2 = 2*Ly;

// Time
Real dt = 1.0/512/2;
Real half_dt = dt/2;
const int cell_update_period = 32;
const int saving_period = 16;
const long int equilibrium_step = 0;
const long int total_step = 33554432; // 2^23 = 8388608 ,  2^25 = 33554432

const Real speed = 1;
Real Dc = 0.5; // The noise above which the initial condition is disordered, and below it is polar ordered.
const Real K = 0;

const Real lx_min = (2 + 2*speed*cell_update_period*dt);
const int max_divisor_x = static_cast<int> (Lx_int / lx_min);// must be smaller than Lx2*(1 - 2*cell_update_period*dt);
const int max_divisor_y = static_cast<int> (Ly_int / lx_min);// must be smaller than Ly2*(1 - 2*cell_update_period*dt);
const int divisor_x = max_divisor_x;
const int divisor_y = max_divisor_y;

const Real rv = 1 + (2*speed*dt*(cell_update_period)); // Radius cut off for verlet list

// Parallel Use only
const int tag_max = 32767; // For parallel use only

#ifdef TRACK_PARTICLE
const int track = 59;
BasicDynamicParticle* track_p;
bool flag = false;
#endif

#ifdef COMPARE
//double digits = 100000000000000000;
double digits =   100000000000000000;
#endif

#ifdef DEBUG
int the_node_id = 0;
#endif

#endif
