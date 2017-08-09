#include<iostream>
#include<cstdlib>
#include<vector>

#include"analyze.h"

using namespace std;

SavingVector mb_r_cm, mb_v_cm, sw_r_cm, sw_v_cm;
double mb_xx, mb_xy, mb_yy, sw_xx, sw_xy, sw_yy;
double mb_l, sw_l, sw_spin;

double mb_lambda1, mb_lambda2, mb_Rg, mb_Delta, mb_omega;
double sw_lambda1, sw_lambda2, sw_Rg, sw_Delta, sw_omega;

void Find_Quantities(Scene* s1, Scene* s0)
{
	mb_r_cm.Null();
	mb_v_cm.Null();
	sw_r_cm.Null();
	sw_v_cm.Null();
	mb_xx = mb_yy = mb_xy = mb_l = 0;
	sw_xx = sw_yy = sw_xy = sw_l = 0;

	for (int i = 0; i < s1->Nm; i++)
	{
			SavingVector r, v;
			r = s1->mparticle[i].r;
			v = (r - s0->mparticle[i].r) / (s1->t - s0->t);
			mb_r_cm += r;
			mb_v_cm += v;
			mb_l += r.x * v.y - r.y * v.x;
			mb_xx += r.x*r.x;
			mb_yy += r.y*r.y;
			mb_xy += r.x*r.y;
	}
	for (int i = 0; i < s1->Ns; i++)
	{
		SavingVector r, v;
		Real spin, dtheta;
		r = s1->sparticle[i].r;
		v = (r - s0->sparticle[i].r) / (s1->t - s0->t);
		dtheta = s1->sparticle[i].theta - s0->sparticle[i].theta;
		dtheta -= 2*M_PI*((int) floor(dtheta / (2*M_PI)));
		dtheta -= 2*M_PI*((int) floor(dtheta / (M_PI)));
		spin = (dtheta) / (s1->t - s0->t);
		sw_r_cm += r;
		sw_v_cm += v;
		sw_spin += spin;
		sw_l += r.x * v.y - r.y * v.x;
		sw_xx += r.x*r.x;
		sw_yy += r.y*r.y;
		sw_xy += r.x*r.y;
	}

	mb_r_cm.x /= s1->Nm;
	mb_r_cm.y /= s1->Nm;
	mb_v_cm.x /= s1->Nm;
	mb_v_cm.y /= s1->Nm;
	mb_xx = (mb_xx / s1->Nm) - mb_r_cm.x*mb_r_cm.x;
	mb_yy = (mb_yy / s1->Nm) - mb_r_cm.y*mb_r_cm.y;
	mb_xy = (mb_xy / s1->Nm) - mb_r_cm.x*mb_r_cm.y;
	mb_l = (mb_l / s1->Nm) - mb_r_cm.x*mb_v_cm.y + mb_r_cm.y*mb_v_cm.x;

	sw_r_cm.x /= s1->Ns;
	sw_r_cm.y /= s1->Ns;
	sw_v_cm.x /= s1->Ns;
	sw_v_cm.y /= s1->Ns;
	sw_xx = (sw_xx / s1->Ns) - sw_r_cm.x*sw_r_cm.x;
	sw_yy = (sw_yy / s1->Ns) - sw_r_cm.y*sw_r_cm.y;
	sw_xy = (sw_xy / s1->Ns) - sw_r_cm.x*sw_r_cm.y;

	sw_l = (sw_l / s1->Ns) - sw_r_cm.x*mb_v_cm.y + sw_r_cm.y*mb_v_cm.x - mb_r_cm.x*sw_v_cm.y + mb_r_cm.y*sw_v_cm.x + mb_r_cm.x*mb_v_cm.y - mb_r_cm.y*mb_v_cm.x;
	sw_spin /= s1->Ns;

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
}

void Write(SceneSet* sceneset, ofstream& os)
{
	os << "#\ttime\tx_cm\ty_cm\tvx_cm\tvy_cm\tangular momentum\tangular freq.\tGyration raduis\tAsphericity\tQxx\tQxy\tQyy\tsw. ang. freq." << endl;
		for (int i = 1; i < sceneset->Nf; i++)
	{
		Find_Quantities(&sceneset->scene[i], &sceneset->scene[i-1]);
		os << sceneset->scene[i].t << "\t" << mb_r_cm << "\t" << mb_v_cm << "\t" << mb_l << "\t" << mb_omega << "\t" << mb_Rg << "\t" << mb_Delta << "\t" << mb_xx << "\t" << mb_xy << "\t" << mb_yy << "\t" << sw_omega << endl;	
	}
}

void Read_Write(string name)
{
		SceneSet* sceneset = new SceneSet(name);
		bool read_state = sceneset->Read();
		if (read_state)
		{
			boost::replace_all(name, ".bin", ".dat");
			boost::replace_all(name, "r-v-", "quantities-");
			stringstream ss("");
			ss << name;
			ofstream out_file;
			out_file.open(ss.str().c_str());
			Write(sceneset, out_file);

			delete sceneset;
		}
		else
			cout << "Was not able to open file: " << name << endl;
}

int main(int argc, char** argv)
{
	SavingVector::Init_Rand(321);
	srand(time(NULL));
	for (int i = 1; i < argc; i++)
	{
		string name = argv[i];
		Read_Write(name);
	}

	return 0;
}

