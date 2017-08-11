#ifndef _BOX_
#define _BOX_

#include "../shared/parameters.h"
#include "../shared/c2dvector.h"
#include "../shared/particle.h"
#include "../shared/cell.h"
#include "../shared/wall.h"
#include "../shared/set-up.h"
#include "../shared/state-hyper-vector.h"

#include <boost/algorithm/string.hpp>

class Box{
public:
	int N, Ns, Nm; // N = Ns+Nm is the number of particles, Nm is number of membrane particles, Ns is number of active particles
	Particle* particle; // Array of particles that we are going to simulate.

	Real Lbx, Lby;
	Real t;
	Real packing_fraction;
	Real cover_fraction;
	Real membrane_elasticity;
	Real membrane_radius;

	Real t_old; // old time, the previous time that coordinates are saved for further analysis
	Real* theta_old; // old angle of self-propelsion
	C2DVector* r_old; // old position of beads
	C2DVector mb_r_cm, mb_v_cm, sw_r_cm, sw_v_cm;
	Real mb_xx, mb_yy, mb_xy, sw_xx, sw_yy, sw_xy; // Qxx, Qyy, Qxy. Q = 1/N sum (\vec{r} - \vec{r}_cm) (\vec{r} - \vec{r}_cm)
	Real mb_l, sw_l; // Angular momentum per particle
	Real mb_omega, sw_omega, sw_spin; // Angular velocity
	Real mb_lambda1, mb_lambda2, sw_lambda1, sw_lambda2; // Eigenvalues of Q tensor
	Real mb_Rg, sw_Rg; // Gyration radius
	Real mb_Delta, sw_Delta; // Asphericity

	stringstream info; // information stream that contains the simulation information, like noise, density and etc. this will be used for the saving name of the system.

	Cell** cell;

	Box();
	Box(const Real);
	// I believe that it is better to move these init functions to main files
	void Init(int, int); // Intialize the box, positioning particles, giving them velocities, and updating cells.
	bool Positioning_Particles(const string name); // Intialize the box from a file, this includes reading particles information, and updating cells.
	void Update_Cells();
	void Update_Neighbor_List(); // This function will update verlet neighore list of particles
	void Sync();

	void Interact_Membrane_Beads(); // Here the wall particles interact via a spring
	void Interact(); // Here the intractio of particles are computed that is the applied tourque to each particle.
	void Move(); // Move all particles of the box.
	void Move_Runge_Kutta2_1(); // Do the first step of Runge Kutta
	void Move_Runge_Kutta2_2(); // Do the second step of Runge Kutta
	void One_Step(); // One full step, composed of interaction computation and move.
	void Multi_Step(int steps); // Several steps befor a cell upgrade.
	void Multi_Step(int steps, int interval); // Several steps with a cell upgrade call after each interval.
	void Translate(C2DVector d); // Translate position of all particles with vector d

	void RootGather(); // Gather all needed data for computing variables such as center of mass, angular momentum, angular velocity, and ...
	void Save_Particles_Positions(); // Save the position of membrane beads to compute their velocities in future.
	void Compute_All_Variables(); // Compute center of mass position and speed, I, ...
	void Save_All_Variables(std::ostream& os); // Save center of mass position and speed, I, ...

	friend std::ostream& operator<<(std::ostream& os, Box* box); // Save
	friend std::istream& operator>>(std::istream& is, Box* box); // Input
};

Box::Box()
{
	Ns = 0;
	Nm = 0;
	N = 0;
	particle = new Particle[max_N];
	r_old = new C2DVector[max_N];
	theta_old = new Real[max_N];

	Cell::particle = particle;
	cell = new Cell*[divisor_x];
	for (int x = 0; x < divisor_x; x++)
		cell[x] = new Cell[divisor_y];
}

Box::Box(const Real input_R)
{
	Ns = 0;
	Nm = 0;

	Lbx = input_R*2;
	Lby = Lbx;

	Lx = Lbx;
	Ly = Lbx;
	Lx2 = 2*Lx;
	Ly2 = 2*Ly;

	Real input_rho = 4 / (M_PI*Particle::sigma_p*Particle::sigma_p);
	max_N = (int) floor(M_PI*input_R*input_R);

	particle = new Particle[max_N];
	r_old = new C2DVector[max_N];
	theta_old = new Real[max_N];

	lx_min = (1 + 2*speed*cell_update_period*dt);
	rv = 2 + 2*speed*cell_update_period*dt;

	max_divisor_x = static_cast<int> (floor(Lx / lx_min));// must be smaller than Lx2*(1 - 2*cell_update_period*dt);
	max_divisor_y = static_cast<int> (floor(Ly / lx_min));// must be smaller than Ly2*(1 - 2*cell_update_period*dt);
	divisor_x = max_divisor_x;
	divisor_y = max_divisor_y;

	cell = new Cell*[divisor_x];
	for (int i = 0; i < divisor_x; i++)
		cell[i] = new Cell[divisor_y];

	for (int i = 0; i < divisor_x; i++)
		for (int j = 0; j < divisor_y; j++)
			cell[i][j].Init((Real) Lx*(2*i-divisor_x + 0.5)/divisor_x, (Real) Ly*(2*j-divisor_y + 0.5)/divisor_y); // setting the center position of each cell

	Cell::particle = particle;

	t = 0;
}

void Box::Update_Cells()
{
	for (int x = 0; x < divisor_x; x++)
		for (int y = 0; y < divisor_y; y++)
			cell[x][y].Delete();

	for (int i = 0; i < N; i++)
	{
		int x,y;
		x = (int) (particle[i].r.x + Lx)*divisor_x / Lx2;
		y = (int) (particle[i].r.y + Ly)*divisor_y / Ly2;

		#ifdef DEBUG
		if ((x >= divisor_x) || (x < 0) || (y >= divisor_y) || (y < 0))
		{
			cout << "\n Particle number " << i << " is Out of the box" << endl << flush;
			cout << "Particle Position is " << particle[i].r << endl;
			exit(0);
		}
		#endif

		cell[x][y].Add(i);
	}

	#ifdef verlet_list
		Update_Neighbor_List();
	#endif
}

// This function will update verlet neighore list of particles
void Box::Update_Neighbor_List()
{
// Each cell must interact with itself and 4 of its 8 neihbors that are right cell, up cell, righ up and right down. Because each intertion compute the torque to both particles we need to use 4 of the 8 directions.
	for (int x = 0; x < divisor_x; x++)
		for (int y = 0; y < divisor_y; y++)
		{
			cell[x][y].Clear_Neighbor_List();
			// Self interaction
			cell[x][y].Neighbor_List();
		}

#ifdef PERIODIC_BOUNDARY_CONDITION
// right, up and up right cells:
// The righmost cells and top cells must be excluded to avoid nieghbor cells interactions.
	for (int x = 0; x < divisor_x; x++)
		for (int y = 0; y < divisor_y; y++)
		{
			cell[x][y].Neighbor_List(&cell[(x+1)%divisor_x][y]);
			cell[x][y].Neighbor_List(&cell[x][(y+1)%divisor_y]);
			cell[x][y].Neighbor_List(&cell[(x+1)%divisor_x][(y+1)%divisor_y]);
			cell[x][y].Neighbor_List(&cell[(x+1)%divisor_x][(y-1+divisor_y)%divisor_y]);
		}
#else
// right, up and up right cells:
// The righmost cells and top cells must be excluded to avoid nieghbor cells interactions.
	for (int x = 0; x < divisor_x-1; x++)
		for (int y = 0; y < divisor_y-1; y++)
		{
			cell[x][y].Neighbor_List(&cell[x+1][y]);
			cell[x][y].Neighbor_List(&cell[x][y+1]);
			cell[x][y].Neighbor_List(&cell[x+1][y+1]);
		}
	for (int x = 0; x < (divisor_x-1); x++)
		cell[x][divisor_y-1].Neighbor_List(&cell[x+1][divisor_y-1]);
	for (int y = 0; y < (divisor_y-1); y++)
		cell[divisor_x-1][y].Neighbor_List(&cell[divisor_x-1][y+1]);
	// right down cell:
// The righmost cells and buttom cells must be excluded to avoid nieghbor cells interactions.
	for (int x = 0; x < (divisor_x-1); x++)
		for (int y = 1; y < divisor_y; y++)
			cell[x][y].Neighbor_List(&cell[x+1][y-1]);
#endif
}


void Box::Sync()
{
	Update_Cells();
}


// Intialize the box, positioning particles, giving them velocities, and updating cells.
void Box::Init(const int input_Ns, const int input_Nm)
{
	#ifdef TRACK_PARTICLE
		track_p = &particle[track];
	#endif

	Ns = input_Ns;
	Nm = input_Nm;
	N = Ns + Nm;

// Positioning the particles at first time. Note that, the positions can be tunned in the main file as well
	cout << "number_of_swimmer = " << Ns << endl; // Printing number of particles.
	cout << "number_of_membrane = " << Nm << endl; // Printing number of particles.
//		Triangle_Lattice_Formation(particle, N, 1);

// Buliding up info stream. In next versions we will take this part out of box, making our libraries more abstract for any simulation of SPP.
	info.str("");
}


// Intialize the box from a file, this includes reading particles information, and updating cells. 
bool Box::Positioning_Particles(const string input_name)
{
	#ifdef TRACK_PARTICLE
		track_p = &particle[track];
	#endif

	string name = input_name;
	stringstream address(name);

	stringstream command("");
	command << "./fix-file.out ";
	command << name;
	system(command.str().c_str());

	ifstream is;
	is.open(address.str().c_str(), fstream::in);
	if (!is.is_open())
		return false;

	is >> this;

	if (Ns < 0 || Ns > 1000000)
		return (false);

	is.seekg(0,ios_base::end);
	int end_of_is = is.tellg();
	is.seekg(0,ios_base::beg);
	while (is.tellg() < end_of_is && is.tellg() >= 0)
	{
		static int counter = 0;
		Save_Particles_Positions();
		is >> this;
		counter++;
	}

	for (int i = 0; i < Nm; i++)
	{
		particle[i].r = particle[i].r_original;
		particle[i].r.Periodic_Transform();
	}

	for (int i = Nm; i < Nm+Ns; i++)
	{
		particle[i].r = particle[i].r_original;
		particle[i].r.Periodic_Transform();
		particle[i].v.x = cos(particle[i].theta);
		particle[i].v.y = sin(particle[i].theta);
		particle[i].Reset();
	}

	is.close();

	Sync();

	return (true);
}


void Box::Interact_Membrane_Beads()
{
	for (int i = 0; i < Nm; i++)
	{
		C2DVector dr = particle[i].r - particle[(i+1)%Nm].r;
		dr.Periodic_Transform();
		Real d = sqrt(dr.Square());
		C2DVector f = Spring(dr, d, Particle::repulsion_radius_m, membrane_elasticity);
		particle[i].f += f;
		particle[(i+1)%Nm].f -= f;
	}
}

// Here the intractio of particles are computed that is the applied tourque to each particle.
void Box::Interact()
{
	#ifdef verlet_list
// with verlet list:
		for (int x = 0; x < divisor_x; x++)
			for (int y = 0; y < divisor_y; y++)
				cell[x][y].Interact();
	#else
// without verlet list:
		for (int x = 0; x < divisor_x; x++)
			for (int y = 0; y < divisor_y; y++)
			{
				cell[x][y].Self_Interact();
				cell[x][y].Interact(&cell[(x+1)%divisor_x][y]);
				cell[x][y].Interact(&cell[x][(y+1)%divisor_y]);
				cell[x][y].Interact(&cell[(x+1)%divisor_x][(y+1)%divisor_y]);
				cell[x][y].Interact(&cell[(x+1)%divisor_x][(y-1+divisor_y)%divisor_y]);
			}
	#endif
}

// Move all particles of the box.
void Box::Move()
{
	for (int i = 0; i < N; i++)
		particle[i].Move();
}

#ifdef RUNGE_KUTTA2
// Do the first step of Runge Kutta
void Box::Move_Runge_Kutta2_1()
{
	for (int i = 0; i < N; i++)
		particle[i].Move_Runge_Kutta2_1();
}

// Do the second step of Runge Kutta
void Box::Move_Runge_Kutta2_2()
{
	for (int i = 0; i < N; i++)
		particle[i].Move_Runge_Kutta2_2();
}
#endif

// One full step, composed of interaction computation and move.
void Box::One_Step()
{
	#ifndef RUNGE_KUTTA2
		Interact_Membrane_Beads();
		Interact();
		Move();
		MPI_Barrier(MPI_COMM_WORLD);
	#else
		Interact_Membrane_Beads();
		Interact();
		Move_Runge_Kutta2_1();

		Interact_Membrane_Beads();
		Interact();
		Move_Runge_Kutta2_2();
	#endif
	t += dt;
}

// Several steps befor a cell upgrade.
void Box::Multi_Step(int steps)
{
	for (int i = 0; i < steps; i++)
	{
		#ifndef RUNGE_KUTTA2
			Interact_Membrane_Beads();
			Interact();
			Move();
		#else
			Interact_Membrane_Beads();
			Interact();
			Move_Runge_Kutta2_1();

			Interact_Membrane_Beads();
			Interact();
			Move_Runge_Kutta2_2();
		#endif
	}
	t += dt*steps;
//	cout << "Updating cells" << endl;
	Update_Cells();
//	cout << "Updating finished" << endl;
}

// Several steps with a cell upgrade call after each interval. Warning, I see no cell update function call! I have to fix it!
void Box::Multi_Step(int steps, int interval)
{
	for (int i = 0; i < steps/interval; i++)
		Multi_Step(interval);
	Multi_Step(steps % interval);
	t += interval;
}

// Translate position of all particles with vector d
void Box::Translate(C2DVector d)
{
	for (int i = 0; i < N; i++)
	{
		particle[i].r += d;
		particle[i].r.Periodic_Transform();
	}
	Update_Cells();
}

void Box::Save_Particles_Positions()
{
	t_old = t;
	for (int i = 0; i < N; i++)
	{
		r_old[i] = particle[i].r_original;
		theta_old[i] = particle[i].theta;
	}
}

void Box::Compute_All_Variables() // Compute center of mass position and speed, I, ...
{
	mb_r_cm.Null();
	mb_v_cm.Null();
	sw_r_cm.Null();
	sw_v_cm.Null();
	mb_xx = mb_yy = mb_xy = mb_l = 0;
	sw_xx = sw_yy = sw_xy = sw_l = 0;

	for (int i = 0; i < Nm; i++)
	{
			C2DVector r, v;
			r = particle[i].r_original;
			v = (r - r_old[i]) / (t - t_old);
			mb_r_cm += r;
			mb_v_cm += v;
			mb_l += r.x * v.y - r.y * v.x;
			mb_xx += r.x*r.x;
			mb_yy += r.y*r.y;
			mb_xy += r.x*r.y;
	}
	for (int i = Nm; i < Ns+Nm; i++)
	{
		C2DVector r, v;
		Real spin, dtheta;
		r = particle[i].r_original;
		v = (r - r_old[i]) / (t - t_old);
		dtheta = particle[i].theta - theta_old[i];
		dtheta -= 2*M_PI*((int) floor(dtheta / (2*M_PI)));
		dtheta -= 2*M_PI*((int) floor(dtheta / (M_PI)));
		spin = (dtheta) / (t - t_old);
		sw_r_cm += r;
		sw_v_cm += v;
		sw_spin += spin;
		sw_l += r.x * v.y - r.y * v.x;
		sw_xx += r.x*r.x;
		sw_yy += r.y*r.y;
		sw_xy += r.x*r.y;
	}

	mb_r_cm.x /= Nm;
	mb_r_cm.y /= Nm;
	mb_v_cm.x /= Nm;
	mb_v_cm.y /= Nm;
	mb_xx = (mb_xx / Nm) - mb_r_cm.x*mb_r_cm.x;
	mb_yy = (mb_yy / Nm) - mb_r_cm.y*mb_r_cm.y;
	mb_xy = (mb_xy / Nm) - mb_r_cm.x*mb_r_cm.y;
	mb_l = (mb_l / Nm) - mb_r_cm.x*mb_v_cm.y + mb_r_cm.y*mb_v_cm.x;

	sw_r_cm.x /= Ns;
	sw_r_cm.y /= Ns;
	sw_v_cm.x /= Ns;
	sw_v_cm.y /= Ns;
	sw_xx = (sw_xx / Ns) - sw_r_cm.x*sw_r_cm.x;
	sw_yy = (sw_yy / Ns) - sw_r_cm.y*sw_r_cm.y;
	sw_xy = (sw_xy / Ns) - sw_r_cm.x*sw_r_cm.y;

	sw_l = (sw_l / Ns) - sw_r_cm.x*mb_v_cm.y + sw_r_cm.y*mb_v_cm.x - mb_r_cm.x*sw_v_cm.y + mb_r_cm.y*sw_v_cm.x + mb_r_cm.x*mb_v_cm.y - mb_r_cm.y*mb_v_cm.x;
	sw_spin /= Ns;

	mb_lambda1 = 0.5*(mb_xx + mb_yy) + 0.5*sqrt((mb_xx - mb_yy)*(mb_xx - mb_yy) + 4*mb_xy*mb_xy);
	mb_lambda2 = 0.5*(mb_xx + mb_yy) - 0.5*sqrt((mb_xx - mb_yy)*(mb_xx - mb_yy) + 4*mb_xy*mb_xy);

	mb_Rg = sqrt(mb_xx + mb_yy);
	mb_Delta = (mb_lambda1 - mb_lambda2) / (mb_lambda1 + mb_lambda2);
	mb_Delta = mb_Delta*mb_Delta;
	mb_omega = mb_l / (mb_xx + mb_yy);

	sw_lambda1 = 0.5*(sw_xx + sw_yy) + 0.5*sqrt((sw_xx - sw_yy)*(sw_xx - sw_yy) + 4*sw_xy*sw_xy);
	sw_lambda2 = 0.5*(sw_xx + sw_yy) - 0.5*sqrt((sw_xx - sw_yy)*(sw_xx - sw_yy) + 4*sw_xy*sw_xy);

	sw_Rg = sqrt(sw_xx + sw_yy);
	sw_Delta = (sw_lambda1 - sw_lambda2) / (sw_lambda1 + sw_lambda2);
	sw_Delta = sw_Delta * sw_Delta;
	sw_omega = sw_l / (sw_xx + sw_yy + (mb_r_cm - sw_r_cm).Square());

	Save_Particles_Positions();
}

void Box::Save_All_Variables(std::ostream& os) // Save center of mass position and speed, I, ...
{
	Compute_All_Variables();
	static bool first_time = true;
	if (first_time && os.tellp()==0)
		os << "#\ttime\tx_cm\ty_cm\tvx_cm\tvy_cm\tangular momentum\tangular freq.\tGyration raduis\tAsphericity\tQxx\tQxy\tQyy\tsw. ang. freq." << endl;
	os << t << "\t" << mb_r_cm << "\t" << mb_v_cm << "\t" << mb_l << "\t" << mb_omega << "\t" << mb_Rg << "\t" << mb_Delta << "\t" << mb_xx << "\t" << mb_xy << "\t" << mb_yy << "\t" << sw_omega << endl;
	first_time = false;
}

// Saving the particle information (position and velocities) to a standard output stream (probably a file). This must be called by only the root.
std::ostream& operator<<(std::ostream& os, Box* box)
{
	os.write((char*) &(box->t), sizeof(Real) / sizeof(char));
	os.write((char*) &(Lx), sizeof(Real) / sizeof(char));
	os.write((char*) &(Ly), sizeof(Real) / sizeof(char));
	os.write((char*) &(box->particle[box->Nm].nb), sizeof(int) / sizeof(char));
	os.write((char*) &(box->Ns), sizeof(box->Ns) / sizeof(char));
	os.write((char*) &(box->Nm), sizeof(box->Nm) / sizeof(char));

	for (int i = 0; i < box->Nm; i++)
	{
		box->particle[i].r_original.write(os);
	}

	for (int i = box->Nm; i < box->N; i++)
	{
		box->particle[i].Write(os);
//			cout << i << "\t" << setprecision(100) << box->particle[i].theta << endl;
	}
}

// Reading the particle information (position and velocities) from a standard input stream (probably a file).
std::istream& operator>>(std::istream& is, Box* box)
{
	int nb; // nb is number of beads for swimmers
	is.read((char*) &(box->t), sizeof(Real) / sizeof(char));
	is.read((char*) &(box->Lbx), sizeof(Real) / sizeof(char));
	is.read((char*) &(box->Lby), sizeof(Real) / sizeof(char));
	is.read((char*) &(nb), sizeof(int) / sizeof(char));
	is.read((char*) &(box->Ns), sizeof(box->Ns) / sizeof(char));
	is.read((char*) &(box->Nm), sizeof(box->Nm) / sizeof(char));

	for (int i = 0; i < box->Nm; i++)
	{
		is >> box->particle[i].r_original;
	}

	for (int i = box->Nm; i < (box->Nm + box->Ns); i++)
	{
		box->particle[i].Read(is);
	}
}

#endif
