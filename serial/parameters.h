#ifndef _PARAMETERS_
#define _PARAMETERS_

//#define PERIODIC_BOUNDARY_CONDITION

// This is for checking particles outside of the box
//#define DEBUG
// This track an specific particle. At the buttom the id of the tracking particle is given
//#define TRACK_PARTICLE
// This will round torques to avoid any difference of this program and other versions caused by truncation of numbers (when the order of a sum is changed the result will change because of the truncation error)
//#define COMPARE

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>
#include <ctime>
#include <algorithm>


using namespace std;

class VicsekParticle;
class ContinuousParticle;

typedef double Real;
typedef ContinuousParticle Particle;
//typedef VicsekParticle Particle;

int seed = 1241;

const Real PI = 3.14159265359;

const int max_wall_num = 100;
const int max_N = 50000;
const int Lx_int = 60;
const int Ly_int = 60;
const int L_int = Lx_int;
const Real Lx = Lx_int;
const Real Ly = Ly_int;
const Real L = L_int;
const Real Lx2 = 2*Lx;
const Real Ly2 = 2*Ly;
const Real L2 = 2*L;
const int max_divisor_x = 20*Lx_int/11;
const int max_divisor_y = 20*Ly_int/11;
const int divisor_x = max_divisor_x;
const int divisor_y = max_divisor_y;




const Real dt = 0.005;
const Real half_dt = dt/2;
const int cell_update_period=20;
const long int equilibrium_step = 5000;
const long int total_step = 5000;
const int saving_period = 20;

const Real repulsion_strength = 5;
const Real wall_repulsion_strength = 10*repulsion_strength;
const Real sigma = 1;
const Real wall_sigma = wall_sigma;

#ifdef TRACK_PARTICLE
const int track = 2;
Particle* track_p;
bool flag = false;
#endif

#endif
