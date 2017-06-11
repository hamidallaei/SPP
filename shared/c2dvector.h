#ifndef _VEC_TEMPLATE_H_
#define _VEC_TEMPLATE_H_

#include <cmath>
#include <iostream>
#include "parameters.h"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

template<typename Type = Real>
class Vec_template
{
	public :
		static const gsl_rng_type * T;
		static gsl_rng * gsl_r;
		Type x,y;

		static void Init_Rand (long int seed)
		{
			gsl_rng_env_setup();
			T = gsl_rng_default;
			gsl_r = gsl_rng_alloc (T);
			gsl_rng_set(gsl_r, seed);
		}

		Vec_template () // This constructor will generate a null vector
		{
			x = y = 0;
		}
		Vec_template (const Type amplitude) // This constructor will initialize the vector as random vector with given magnitude (amplitude)
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

		void Rand(const Type amplitude) // Choose random numbers for components of the vector in [-L, L]
		{
			x = gsl_ran_flat (gsl_r, -amplitude, amplitude);
			y = gsl_ran_flat (gsl_r, -amplitude, amplitude);
		}

		void Rand(const Type amplitude_x, const Type amplitude_y) // Choose random numbers for components of the vector in [-Lx, Lx] for x and [-Ly, Ly] for y
		{
			x = gsl_ran_flat (gsl_r, -amplitude_x, amplitude_x);
			y = gsl_ran_flat (gsl_r, -amplitude_y, amplitude_y);
		}

		void Rand_Lattice()
		{
			x = (Type) gsl_rng_uniform_int(gsl_r, (int) 2*Lx) - Lx;
			y = (Type) gsl_rng_uniform_int(gsl_r, (int) 2*Ly) - Ly;
		}


		void Periodic_Transform()
		{
//			x -= Lx2*((int) (x / Lx));
//			y -= Ly2*((int) (y / Ly));
			x -= Lx2*((int) floor(x / Lx2 + 0.5));
			y -= Ly2*((int) floor(y / Ly2 + 0.5));
		}

		Type Square() const// returns the magnitude of the vector
		{
			return (x*x + y*y);
		}

		void Unit() // return unit vector parallel to the vector
		{
			Type magnitude = sqrt(x*x + y*y);
			x /= magnitude;
			y /= magnitude;
		}

		Vec_template Rotate(Type phi)
		{
			Vec_template result;
			result.x = cos(phi)*x - sin(phi)*y;
			result.y = sin(phi)*x + cos(phi)*y;
			return (result);
		}
		
		Vec_template operator+ (const Vec_template p1) const
		{
			Vec_template result;

			result.x = x + p1.x;
			result.y = y + p1.y;

			return result;
		}
		
		Vec_template operator- (const Vec_template p1) const
		{
			Vec_template result;

			result.x = x - p1.x;
			result.y = y - p1.y;

			return result;
		}

		Vec_template operator/ (const Type lambda) const
		{
			Vec_template result;

			result.x = x / lambda;
			result.y = y / lambda;

			return result;
		}
		
		Vec_template operator* (const Type lambda) const
		{
			Vec_template result;

			result.x = x * lambda;
			result.y = y * lambda;

			return result;
		}

		Vec_template operator+= (const Vec_template p1)
		{
			x += p1.x;
			y += p1.y;

			return *this;
		}

		Vec_template operator-= (const Vec_template p1)
		{

			x -= p1.x;
			y -= p1.y;

			return *this;
		}

		Vec_template operator/= (const Type lambda)
		{

			x /= lambda;
			y /= lambda;

			return *this;
		}

		Vec_template operator*= (const Type lambda)
		{
			x *= lambda;
			y *= lambda;

			return *this;
		}

		Type operator* (const Vec_template p1) const
		{
			return (x*p1.x + y*p1.y);
		}

		// for binary output only (comment this function for txt output)
		void write(std::ostream& os)
		{
			Saving_Real temp_float;

			temp_float = (Saving_Real) x;
			os.write((char*) &temp_float,sizeof(Saving_Real) / sizeof(char));

			temp_float = (Saving_Real) y;
			os.write((char*) &temp_float,sizeof(Saving_Real) / sizeof(char));
		}

		// for either binary or txt output
		friend std::ostream& operator<<(std::ostream& os, const Vec_template t)
		{
			os << t.x << "\t" << t.y;
			return (os);
		}

		friend std::istream& operator>>(std::istream& is, Vec_template& t)
		{
			// for binary input 
			Saving_Real temp_float;

			is.read((char*) &temp_float,sizeof(Saving_Real) / sizeof(char));
			t.x = temp_float;

			is.read((char*) &temp_float,sizeof(Saving_Real) / sizeof(char));
			t.y = temp_float;
			return is;

//			// for txt input
//			is >> t.x;
//			is >> t.y;

//			return is;

		}

		~Vec_template()
		{
		}
};

template<typename Type> const gsl_rng_type * Vec_template<Type>::T;
template<typename Type> gsl_rng * Vec_template<Type>::gsl_r;

/*class Index{*/
/*public:*/
/*	int x,y;*/
/*	void Find(Vec_template<> r)*/
/*	{*/
/*		x = (int) (r.x + Lx)*divisor_x / (Lx2);*/
/*		y = (int) (r.y + Ly)*divisor_y / (Ly2);*/
/*	}*/
/*};*/


typedef Vec_template<Real> C2DVector;
typedef Vec_template<Saving_Real> SavingVector;

#endif
