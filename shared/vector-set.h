#ifndef _VECTOR_SET_
#define _VECTOR_SET_

#include "c2dvector.h"
#include "state-hyper-vector.h"


class VectorSet{
public:
	int direction_num, particle_num;
	static Real amplitude;
	vector<State_Hyper_Vector> v;

	VectorSet(int vector_number, int particle_per_vector);
	VectorSet(VectorSet& vs);
	~VectorSet();

	void Rand();
	void Renormalize();
	void Renormalize(VectorSet& us);
	void Unit_Vector(VectorSet& us);
	void Init();
	void Scale(); // Scale all vectors by amplitude
	
	VectorSet& operator= ( const VectorSet& vs);
	VectorSet& operator*= ( const Real& factor);
	const VectorSet operator* ( const Real& factor);
};

VectorSet::VectorSet(int vector_number, int particle_per_vector) : direction_num(vector_number), particle_num(particle_per_vector)
{
	// Now the seed is a pice of trash! It will be removed in future
	int* seed;
	seed = new int[direction_num];
	bool b = false;
	while (!b)
	{
		b = true;
		for (int i = 0; i < direction_num; i++)
			seed[i] = i + time(NULL);
		for (int i = 0; i < direction_num; i++)
			for (int j = i+1; j < direction_num; j++)
				b = b && (seed[i] != seed[j]);
	}
	for (int i = 0; i < direction_num; i++)
	{
		State_Hyper_Vector temp(particle_num, seed[i]);
		v.push_back(temp);
	}
	delete [] seed;
}

VectorSet::VectorSet(VectorSet& vs)
{
	direction_num = vs.direction_num;
	particle_num = vs.particle_num;
	amplitude = vs.amplitude;
	for (int i = 0; i < direction_num; i++)
		v.push_back(vs.v[i]);
}

VectorSet::~VectorSet()
{
	v.clear();
}

VectorSet& VectorSet::operator= ( const VectorSet& vs)
{
	if (vs.direction_num != direction_num)
	{
		cout << "Error, number of vectors in set is not the same" << endl;
	}
	for (int i = 0; i < direction_num; i++)
		v[i] = vs.v[i];
	return *this;
}

VectorSet& VectorSet::operator*= ( const Real& factor)
{
	for (int i = 0; i < direction_num; i++)
		v[i] = v[i]*factor;
	return *this;
}

const VectorSet VectorSet::operator* ( const Real& factor)
{
	VectorSet result(*this);
	for (int i = 0; i < direction_num; i++)
		result.v[i] = v[i]*factor;
	return result;
}

void VectorSet::Scale()
{
	for (int i = 0; i < direction_num; i++)
		v[i] = v[i]*amplitude;
}

void VectorSet::Rand()
{
	State_Hyper_Vector temp(particle_num, time(NULL));
	for (int i = 0; i < v.size(); i++)
	{
		temp.Rand(0.001, 0.001);
		v[i] = temp;
	}
	Renormalize();
	Scale();
}

void VectorSet::Renormalize()
{
	State_Hyper_Vector temp(particle_num);
	for (int i = 0; i < v.size(); i++)
	{
		for (int j = 0; j < i; j++)
			v[i] -= v[j]*(v[i]*v[j]);
		v[i].Unit();
	}
}

void VectorSet::Renormalize(VectorSet& us)
{
	for (int i = 0; i < v.size(); i++)
	{
		us.v[i] = v[i];
		for (int j = 0; j < i; j++)
 			us.v[i] -= us.v[j]*(us.v[j]*us.v[i]);
		us.v[i] = us.v[i] / us.v[i].Magnitude();
	}
}

Real VectorSet::amplitude = 1e-7;

#endif

