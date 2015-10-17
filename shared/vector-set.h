#ifndef _VECTOR_SET_
#define _VECTOR_SET_

#include "c2dvector.h"
#include "state-hyper-vector.h"


class VectorSet{
public:
	static int direction_num, particle_num;
	static Real amplitude;
	vector<State_Hyper_Vector> v;

	VectorSet();
	VectorSet(const int vector_number, const int particle_per_vector);
	VectorSet(VectorSet& vs);
	~VectorSet();

	void Init();
	void Null();
	void Rand();
	void Renormalize();
	void Renormalize(VectorSet& us);
	void Unit_Vector(VectorSet& us);
	void Scale(); // Scale all vectors by amplitude
	
	VectorSet& operator= ( const VectorSet& vs);
	VectorSet& operator*= ( const Real& factor);
	const VectorSet operator* ( const Real& factor);

	friend std::ostream& operator<<(std::ostream& os, const VectorSet& vs); // Save
};

VectorSet::VectorSet()
{
	Init();
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

void VectorSet::Init()
{
	v.clear();
	for (int i = 0; i < direction_num; i++)
	{
		State_Hyper_Vector temp(particle_num);
		v.push_back(temp);
	}
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
	v[0].Null();
	for (int i = 1; i < v.size(); i++)
	{
		temp.Rand(0.001, 0.001);
		v[i] = temp;
	}
	Renormalize();
}

void VectorSet::Renormalize()
{
	State_Hyper_Vector temp(particle_num);
	for (int i = 1; i < v.size(); i++)
	{
		for (int j = 1; j < i; j++)
			v[i] -= v[j]*(v[i]*v[j]);
		v[i].Unit();
	}
}

void VectorSet::Renormalize(VectorSet& us)
{
	for (int i = 1; i < v.size(); i++)
	{
		us.v[i] = v[i];
		for (int j = 1; j < i; j++)
 			us.v[i] -= us.v[j]*(us.v[j]*us.v[i]);
		us.v[i] = us.v[i] / us.v[i].Magnitude();
	}
}

int VectorSet::direction_num = 0;
int VectorSet::particle_num = 0;
Real VectorSet::amplitude = 1e-7;

class GrowthRatio{
public:
	static int direction_num;
	int num;
	vector<Real> r;
	vector<Real> r2;

	GrowthRatio();
	GrowthRatio(const GrowthRatio& gr);

	~GrowthRatio();

	void Init();
	

	
};

GrowthRatio::GrowthRatio()
{
	num = 0;
	for (int i = 0; i < direction_num; i++)
	{
		r.push_back(0);
		r2.push_back(0);
	}
}

GrowthRatio::GrowthRatio(const GrowthRatio& gr)
{
	num = gr.num;
	for (int i = 0; i < direction_num; i++)
	{
		r.push_back(gr.r[i]);
		r2.push_back(gr.r2[i]);
	}
}

GrowthRatio::~GrowthRatio()
{
	r.clear();
	r2.clear();
}

void GrowthRatio::Init()
{
	r.clear();
	r2.clear();
	num = 0;
	for (int i = 0; i < direction_num; i++)
	{
		r.push_back(0);
		r2.push_back(0);
	}
}

std::ostream& operator<<(std::ostream& os, const VectorSet& vs) // Save
{
	os << "np";
	for (int i = 1; i < vs.direction_num; i++)
		os << "\t" << i;
	os << endl;
	long int digits = 10000000;
	for (int k = 0; k < vs.v[0].N; k++)
	{
		os << "p " << k;
		for (int i = 1; i < vs.direction_num; i++)
			os << "\t" << setprecision(10) << round(digits*vs.v[i].particle[k].r.x)/digits;
		os << endl;
		os << "p " << k;
		for (int i = 1; i < vs.direction_num; i++)
			os << "\t" << setprecision(10) << round(digits*vs.v[i].particle[k].r.y)/digits;
		os << endl;
		os << "p " << k;
		for (int i = 1; i < vs.direction_num; i++)
			os << "\t" << setprecision(10) << round(digits*vs.v[i].particle[k].theta)/digits;
		os << endl;
	}
	return (os);
}

int GrowthRatio::direction_num = 0;

#endif
