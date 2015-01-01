#ifndef _WALL_
#define _WALL_

#include "c2dvector.h"

class Wall{
public:
	C2DVector point_1,point_2;
	C2DVector direction, normal;
	Real length, theta;

	Wall();
	void Init(C2DVector input_p_1, C2DVector intput_p_2);
	void Init(Real input_p_1_x, Real input_p_1_y, Real input_p_2_x, Real input_p_2_y);
	Real Distance(C2DVector p);
	Real Intercept_x(Real y_value);
	Real Intercept_y(Real x_value);
	void Translate(C2DVector delta);
	void Rotate(Real phi, C2DVector r);

	void Interact(VicsekParticle* p);
	void Interact(ContinuousParticle* p);
};


Wall::Wall()
{
	length = 0;
}

void Wall::Init(C2DVector input_p_1, C2DVector input_p_2)
{
	point_1 = input_p_1;
	point_2 = input_p_2;
	direction = point_2 - point_1;
	length = sqrt(direction.Square());
	direction = direction / length;
	normal.x = -direction.y;
	normal.y = direction.x;
	theta = atan(direction.y / direction.x);
	if (direction.x < 0)
		theta += PI;
	theta -= 2*PI * (int (theta / (PI)));
}


void Wall::Init(Real input_p_1_x, Real input_p_1_y, Real input_p_2_x, Real input_p_2_y)
{
	point_1.x = input_p_1_x;
	point_1.y = input_p_1_y;
	point_2.x = input_p_2_x;
	point_2.y = input_p_2_y;
	direction = point_2 - point_1;
	length = sqrt(direction.Square());
	direction = direction / length;
	normal.x = -direction.y;
	normal.y = direction.x;
	theta = atan(direction.y / direction.x);
	if (direction.x < 0)
		theta += PI;
	theta -= 2*PI * (int (theta / (PI)));
}

Real Wall::Intercept_x(Real y_value)
{
	Real dy = y_value - point_1.y;
	#ifdef PERIODIC_BOUNDARY_CONDITION
		dy -= L2*((int) (dy / L));
	#endif
	Real s = dy / direction.y;
	if ((s > 0) && (s < length))
	{
		Real result = point_1.x + s * direction.x;
		#ifdef PERIODIC_BOUNDARY_CONDITION
			result -= L2*((int) (result / L));
		#endif
		return (result);
	}
	else
		return (4*L);
}


Real Wall::Intercept_y(Real x_value)
{
	Real dx = x_value - point_1.x;
	#ifdef PERIODIC_BOUNDARY_CONDITION
		dx -= L2*((int) (dx / L));
	#endif
	Real s = dx / direction.x;
	if ((s > 0) && (s < length))
	{
		Real result = point_1.y + s * direction.y;
		#ifdef PERIODIC_BOUNDARY_CONDITION
			result -= L2*((int) (result / L));
		#endif
		return (result);
	}
	else
		return (4*L);
}

void Wall::Translate(C2DVector delta)
{
	point_1 += delta;
	point_2 += delta;
	#ifdef PERIODIC_BOUNDARY_CONDITION
		point_1.Periodic_Transform();
		point_2.Periodic_Transform();
	#endif
}

void Wall::Rotate(Real phi, C2DVector r)
{
	C2DVector dr;
	dr = point_1 - r;
	point_1 = dr.Rotate(phi) + r;
	dr = point_2 - r;
	point_2 = dr.Rotate(phi) + r;

	#ifdef PERIODIC_BOUNDARY_CONDITION
		point_1.Periodic_Transform();
		point_2.Periodic_Transform();
	#endif

	direction = direction.Rotate(phi);
	normal.x = -direction.y;
	normal.y = direction.x;
}


void Wall::Interact(VicsekParticle* p)
{
		C2DVector dr_1,dr_2,dr;
		Real d2_1, d2_2, d2;
		dr_1 = p->r - point_1;
		dr_2 = p->r - point_2;
		#ifdef PERIODIC_BOUNDARY_CONDITION
			dr_1.Periodic_Transform();
			dr_2.Periodic_Transform();
		#endif

		if ((dr_1*direction < length) && (dr_1*direction > 0))
		{
			dr = dr_1 - direction*(dr_1*direction);
			d2 = dr.Square();
			Real dtheta = p->theta - theta;
			dtheta -= 2*PI * (int (dtheta / PI));
			if ((dtheta > PI/2) || (dtheta < -PI/2))
			{
				dtheta += PI;
				dtheta -= 2*PI * (int (dtheta / PI));
			}

			if (d2 < 1)
			{
				p->neighbor_size++;
				p->average_theta += dtheta;
			}
		}
}

void Wall::Interact(ContinuousParticle* p)
{
		C2DVector dr_1,dr;
		Real d2;
		dr_1 = p->r - point_1;
		#ifdef PERIODIC_BOUNDARY_CONDITION
			dr_1.Periodic_Transform();
		#endif

			dr = dr_1 - direction*(dr_1*direction);
			d2 = dr.Square();
			if (d2 < 1)
			{
				p->torque += (40.0/Particle::g)*((p->v.x*dr.y - p->v.y*dr.x)/(2*PI*(d2)));
				#ifdef TRACK_PARTICLE
					if (p == track_p)
					{
		//				if (abs(p->torque) > 0.1)
	//						cout << "Intwall:     " << p->r << "\t" << d2 << "\t" << p->theta << "\t" << p->torque << endl << flush;
					}
				#endif
			}
}

#endif
