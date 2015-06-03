#ifndef _PARAMETERS_
#define _PARAMETERS_

#define PERIODIC_BOUNDARY_CONDITION
//#define CIRCULAR_BOX
#define verlet_list
// This is for checking particles outside of the box
//#define DEBUG
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
class ContinuousParticle;
class MarkusParticle;
class RepulsiveParticle;

typedef double Real;
//typedef VicsekParticle Particle;
//typedef ContinuousParticle Particle;
typedef MarkusParticle Particle;
//typedef RepulsiveParticle Particle;

long int seed = 9348;

const Real PI = M_PI;

const int max_wall_num = 8;
const int max_N = 80000;

// Box
const int Lx_int = 5;
const int Ly_int = 5;
const int L_int = Lx_int;
const Real Lx = Lx_int;
const Real Ly = Ly_int;
const Real Lx2 = 2*Lx;
const Real Ly2 = 2*Ly;
//const Real L = L_int; No need to these variables
//const Real L2 = 2*L; No need to these variables

// Trap
const Real r_big = 15;
const Real r_small = 6;

// Cell division
const int max_divisor_x = 20*Lx_int/11;
const int max_divisor_y = 20*Ly_int/11;
const int divisor_x = max_divisor_x;
const int divisor_y = max_divisor_y;

// Parallel Use only
const int npx = 2; // For parallel use only. This number must be even to avoid dead locks
const int npy = 2; // For parallel use only. This number must be even to avoid dead locks
const int tag_max = 32767; // For parallel use only

// Time
const Real dt = 0.001;
const Real half_dt = dt/2;
const int cell_update_period = 20;
const int saving_period = 20;
const long int equilibrium_step = 10000;
const long int total_step = 1000;

// Interactions
const Real A_p = 1.;		// interaction strength
const Real A_w = 50.;
const Real sigma_p = .5;		// sigma in Yukawa Potential
const Real sigma_w = 1.;
const Real r_f_p = 1.;		// flocking radius with particles
const Real r_f_w = 1.;		// aligning radius with walls
const Real r_c_p = .5; 		// repulsive cutoff radius with particles
const Real r_c_w = 1.; 		// repulsive cutoff radius with walls

#ifdef TRACK_PARTICLE
const int track = 2;
BasicDynamicParticle* track_p;
bool flag = false;
#endif

#ifdef COMPARE
double digits = 10000;
#endif


#endif
