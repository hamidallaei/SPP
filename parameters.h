#ifndef _PARAMETERS_
#define _PARAMETERS_

//#define PERIODIC_BOUNDARY_CONDITION
#define DEBUG

using namespace std;

class VicsekParticle;
class ContinuousParticle;

typedef double Real;
typedef ContinuousParticle Particle;
//typedef VicsekParticle Particle;

int seed = 1241;

const Real PI = 3.14159265359;

const int max_wall_num = 100;
const int max_N = 20000;
const int L_int = 20;
const Real L = L_int;
const Real L2 = 2*L;
const int divisor = 20*L_int/11;
const int grid_size = 40;


const Real dt = 0.005;
const Real half_dt = dt/2;
const int cell_update_period=20;
const int equilibrium_step = 0;
const int total_step = 100000;
const int saving_period = 10;

const Real repulsion_strength = 5;
const Real wall_repulsion_strength = 10*repulsion_strength;
const Real sigma = 1;
const Real wall_sigma = wall_sigma;

const Real high_velocity = 5;
const Real low_velocity = 0.5;

#endif
