#ifndef _C2DVECTOR_H_
#define _C2DVECTOR_H_

#include <cmath>
#include <iostream>
#include "parameters.h"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

class C2DVector
{
	public :
		static const gsl_rng_type * T;
		static gsl_rng * gsl_r;
		Real x,y;

		static void Init_Rand (long int seed)
		{
			gsl_rng_env_setup();
			T = gsl_rng_default;
			gsl_r = gsl_rng_alloc (T);
			gsl_rng_set(gsl_r, seed);
		}

		C2DVector () // This constructor will generate a null vector
		{
			x = y = 0;
		}
		C2DVector (const Real amplitude) // This constructor will initialize the vector as random vector with given magnitude (amplitude)
		{
			x = gsl_ran_flat (gsl_r, -amplitude, amplitude);
			y = gsl_ran_flat (gsl_r, -amplitude, amplitude);
		}

		void Null() // Make all components of the vector zero
		{
			x = y = 0;
		}

		void Rand() // Choose random numbers for components of the vector with x in interval [-Lx, Lx] and y in interval [-Ly, Ly]
		{
			x = gsl_ran_flat (gsl_r, -Lx, Lx);
			y = gsl_ran_flat (gsl_r, -Ly, Ly);
		}

		void Rand(const Real amplitude) // Choose random numbers for components of the vector in [-L, L]
		{
			x = gsl_ran_flat (gsl_r, -amplitude, amplitude);
			y = gsl_ran_flat (gsl_r, -amplitude, amplitude);
		}

		void Rand(const Real amplitude_x, const Real amplitude_y) // Choose random numbers for components of the vector in [-Lx, Lx] for x and [-Ly, Ly] for y
		{
			x = gsl_ran_flat (gsl_r, -amplitude_x, amplitude_x);
			y = gsl_ran_flat (gsl_r, -amplitude_y, amplitude_y);
		}

		void Rand_Lattice()
		{
			x = (Real) gsl_rng_uniform_int(gsl_r, (int) 2*Lx) - Lx;
			y = (Real) gsl_rng_uniform_int(gsl_r, (int) 2*Ly) - Ly;
		}


		void Periodic_Transform()
		{
//			x -= Lx2*((int) (x / Lx));
//			y -= Ly2*((int) (y / Ly));
			x -= Lx2*((int) floor(x / Lx2 + 0.5));
			y -= Ly2*((int) floor(y / Ly2 + 0.5));
		}

		Real Square() const// returns the magnitude of the vector
		{
			return (x*x + y*y);
		}

		void Unit() // return unit vector parallel to the vector
		{
			Real magnitude = sqrt(x*x + y*y);
			x /= magnitude;
			y /= magnitude;
		}

		C2DVector Rotate(Real phi)
		{
			C2DVector result;
			result.x = cos(phi)*x - sin(phi)*y;
			result.y = sin(phi)*x + cos(phi)*y;
			return (result);
		}
		
		C2DVector operator+ (const C2DVector p1) const
		{
			C2DVector result;

			result.x = x + p1.x;
			result.y = y + p1.y;

			return result;
		}
		
		C2DVector operator- (const C2DVector p1) const
		{
			C2DVector result;

			result.x = x - p1.x;
			result.y = y - p1.y;

			return result;
		}

		C2DVector operator/ (const Real lambda) const
		{
			C2DVector result;

			result.x = x / lambda;
			result.y = y / lambda;

			return result;
		}
		
		C2DVector operator* (const Real lambda) const
		{
			C2DVector result;

			result.x = x * lambda;
			result.y = y * lambda;

			return result;
		}

		C2DVector operator+= (const C2DVector p1)
		{
			x += p1.x;
			y += p1.y;

			return *this;
		}

		C2DVector operator-= (const C2DVector p1)
		{

			x -= p1.x;
			y -= p1.y;

			return *this;
		}

		C2DVector operator/= (const Real lambda)
		{

			x /= lambda;
			y /= lambda;

			return *this;
		}

		C2DVector operator*= (const Real lambda)
		{
			x *= lambda;
			y *= lambda;

			return *this;
		}

		Real operator* (const C2DVector p1) const
		{
			return (x*p1.x + y*p1.y);
		}

		// for binary output only (comment this function for txt output)
		void write(std::ostream& os)
		{
			float temp_float;

			temp_float = (float) x;
			os.write((char*) &temp_float,sizeof(float) / sizeof(char));

			temp_float = (float) y;
			os.write((char*) &temp_float,sizeof(float) / sizeof(char));
		}

		// for either binary or txt output
		friend std::ostream& operator<<(std::ostream& os, const C2DVector t)
		{
			os << t.x << "\t" << t.y;
			return (os);
		}

		friend std::istream& operator>>(std::istream& is, C2DVector& t)
		{
			// for binary input 
			float temp_float;

			is.read((char*) &temp_float,sizeof(float) / sizeof(char));
			t.x = temp_float;

			is.read((char*) &temp_float,sizeof(float) / sizeof(char));
			t.y = temp_float;
			return is;

//			// for txt input
//			is >> t.x;
//			is >> t.y;

//			return is;

		}

		~C2DVector()
		{
		}
};

const gsl_rng_type * C2DVector::T;
gsl_rng * C2DVector::gsl_r;

class Index{
public:
	int x,y;
	void Find(C2DVector r)
	{
		x = (int) (r.x + Lx)*divisor_x / (Lx2);
		y = (int) (r.y + Ly)*divisor_y / (Ly2);
	}
};

#endif
