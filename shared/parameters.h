#ifndef _PARAMETERS_
#define _PARAMETERS_

#define PERIODIC_BOUNDARY_CONDITION
//#define NonPeriodicCompute
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
class VicsekParticle2;
class ContinuousParticle;
class MarkusParticle;
class RepulsiveParticle;

typedef double Real;
//typedef VicsekParticle2 Particle;
typedef ContinuousParticle Particle;
//typedef MarkusParticle Particle;
//typedef RepulsiveParticle Particle;

long int seed = 10;

const Real PI = M_PI;

const int max_wall_num = 8;
const int max_N = 20000;

// Box
const int Lx_int = 16;
const int Ly_int = Lx_int;
const int L_int = Lx_int;
const Real Lx = Lx_int;
const Real Ly = Ly_int;
const Real Lx2 = 2*Lx;
const Real Ly2 = 2*Ly;

// Time
Real dt = 1.0/64;
Real half_dt = dt/2;
const int cell_update_period = 4;
const int saving_period = 64;
const long int equilibrium_step = 0;//10000;
const long int total_step = 1024;//120000;

const Real speed = 2;
Real Dc = 0.5; // The noise above which the initial condition is disordered, and below it is polar ordered.

const Real lx_min = (1 + 2*speed*cell_update_period*dt);
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


#endif
