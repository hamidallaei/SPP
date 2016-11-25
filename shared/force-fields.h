#ifndef _FORCEFILED_
#define _FORCEFILED_

#include "../shared/c2dvector.h"

/*
	The result returns the force acting on "this" bead from "other" beads.
	sigma is the characteristic length scale of the force field
	dr = r_this -r_other
	d = magnitude of dr
	cutoff is the cutoff radius
*/

C2DVector R12_Repulsive_Truncated(const C2DVector& dr, const Real& d, const Real& cutoff, const Real& sigma, const Real& amplitude)
{/*
	"Truncated and repulsive part of the Lennard-Jones force between overlapping beads"
	This force is a derivative of the Lennard-Jones potential which is equal to V(r) = A * sigma /12 * (sigma/r)^12
	The force is truncated at cutoff distance which means f = 0 for r >= cutoff
*/
	Real sigma_over_d = sigma / d;
	Real sigma_over_d2 = sigma_over_d*sigma_over_d;
	Real sigma_over_d4 = sigma_over_d2*sigma_over_d2;
	Real sigma_over_d8 = sigma_over_d4*sigma_over_d4;
	Real sigma_over_d13 = sigma_over_d8*sigma_over_d4*sigma_over_d;

	Real sigma_over_cutoff = sigma / cutoff;
	Real sigma_over_cutoff2 = sigma_over_cutoff*sigma_over_cutoff;
	Real sigma_over_cutoff4 = sigma_over_cutoff2*sigma_over_cutoff2;
	Real sigma_over_cutoff8 = sigma_over_cutoff4*sigma_over_cutoff4;
	Real sigma_over_cutoff13 = sigma_over_cutoff8*sigma_over_cutoff4*sigma_over_cutoff;

	Real strength = amplitude* (sigma_over_d13 - sigma_over_cutoff13);

	C2DVector result;
	result.x = strength*dr.x/d;
	result.y = strength*dr.y/d;

	return result;
}

C2DVector R12_Repulsive(const C2DVector& dr, const Real& d, const Real& sigma, const Real& amplitude)
{/*
	"Repulsive part of the Lennard-Jones potential between overlapping beads"
*/
	Real sigma_over_d = sigma / d;
	Real sigma_over_d2 = sigma_over_d*sigma_over_d;
	Real sigma_over_d4 = sigma_over_d2*sigma_over_d2;
	Real sigma_over_d8 = sigma_over_d4*sigma_over_d4;
	Real sigma_over_d13 = sigma_over_d8*sigma_over_d4*sigma_over_d;

	Real strength = amplitude*sigma_over_d13;

	C2DVector result;
	result.x = strength*dr.x/d;
	result.y = strength*dr.y/d;

	return result;
}

C2DVector Yukawa_Truncated(const C2DVector& dr, const Real& d, Real cutoff, const Real& sigma, const Real& amplitude)
{/*
	"Truncated Repulsive Yukawa force field between overlapping beads"
	This force is a derivative of the Yukawa potential which is equal to V(r) = A/r * exp(-r/sigma)
	The force is truncated at cutoff distance which means f = 0 for r >= cutoff
*/
	Real d2 = d*d;
	Real cutoff2 = cutoff*cutoff;

	C2DVector result;
	result.x = amplitude * ( exp(- d / sigma ) * ( 1. / d2 + 1. / (sigma * d)) - exp(- cutoff / sigma ) * ( 1. / cutoff2 + 1. / (sigma * cutoff)) ) * dr.x /d;
	result.y = amplitude * ( exp(- d / sigma ) * ( 1. / d2 + 1. / (sigma * d)) - exp(- cutoff / sigma ) * ( 1. / cutoff2 + 1. / (sigma * cutoff)) ) * dr.y /d;
	return result;
}

C2DVector Yukawa(const C2DVector& dr, const Real& d, const Real& sigma, const Real& amplitude)
{/*
	"Repulsive Yukawa force field between overlapping beads"
*/
	Real d2 = d*d;

	C2DVector result;
	result.x = amplitude * ( exp(- d / sigma ) * ( 1. / d2 + 1. / (sigma * d)) ) * dr.x /d;
	result.y = amplitude * ( exp(- d / sigma ) * ( 1. / d2 + 1. / (sigma * d)) ) * dr.y /d;
	return result;
}

C2DVector Spring(const C2DVector& dr, const Real& d, const Real& sigma, const Real& amplitude)
{/*
	"Spring force field between neighboring beads"
	Here sigma is the relaxed length of the spring.
*/
	C2DVector result;
	result.x = amplitude * (sigma - d) * dr.x/d;
	result.y = amplitude * (sigma - d) * dr.y/d;

	return result;
}



#endif
