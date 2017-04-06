#include<iostream>
#include<cstdlib>
#include<vector>

#include"analyze.h"

using namespace std;

/* 
This code reads position data of the membrane's particles from -r-v.bin files and finds the curvature of each particle/point using its neighbors. "Menger method"
*/


bool negative_sign(C2DVector r1, C2DVector r2, C2DVector r3){
	double sign_value = (r2.x - r1.x)*(r3.y - r1.y) - (r3.x - r1.x)*(r2.y - r1.y);

	bool sign;

	if (sign_value < 0){
		sign = true;
	}
	else
		sign = false;

	return sign;
}

int main(int argc, char** argv)
{
	C2DVector::Init_Rand(321);
	srand(time(NULL));

	stringstream ss("");
	ss << "curvatures.dat";
	ofstream out_file;
	out_file.open(ss.str().c_str());

	for (int i = 1; i < argc; i++)
	{
		string name = argv[i];
		SceneSet* sceneset = new SceneSet(name);
		bool read_state = sceneset->Read();
		if (read_state)
		{
			C2DVector box_dim(sceneset->L);

//			boost::replace_all(name, "-r-v.bin", "");
//			stringstream ss("");
//			ss << "curvatures-" << name << ".dat";
//			ofstream out_file;
//			out_file.open(ss.str().c_str());

			int every = 3;
			for (int j = 200; j < sceneset->scene.size(); j++)
			{
//				int j = 200;
//				cout << sceneset->scene[j].Nm << endl;
				for (int k = 0; k < sceneset->scene[j].Nm; k++)
				{
					if( k == 0){
						double length_a, length_b, length_c, Area, curvature = 0;
						length_a = sqrt( (sceneset->scene[j].mparticle[every].r - sceneset->scene[j].mparticle[k].r).Square() );

						length_b = sqrt( (sceneset->scene[j].mparticle[every*(sceneset->scene[j].Nm/every)-1].r - sceneset->scene[j].mparticle[every].r).Square() );

						length_c = sqrt( (sceneset->scene[j].mparticle[k].r - sceneset->scene[j].mparticle[every*(sceneset->scene[j].Nm/every)-1].r).Square() );

						Area = abs(0.5*( sceneset->scene[j].mparticle[every*(sceneset->scene[j].Nm/every)-1].r.x*sceneset->scene[j].mparticle[k].r.y + sceneset->scene[j].mparticle[k].r.x*sceneset->scene[j].mparticle[every].r.y + sceneset->scene[j].mparticle[every].r.x*sceneset->scene[j].mparticle[every*(sceneset->scene[j].Nm/every)-1].r.y - sceneset->scene[j].mparticle[every*(sceneset->scene[j].Nm/every)-1].r.x*sceneset->scene[j].mparticle[every].r.y - sceneset->scene[j].mparticle[k].r.x*sceneset->scene[j].mparticle[every*(sceneset->scene[j].Nm/every)-1].r.y - sceneset->scene[j].mparticle[every].r.x*sceneset->scene[j].mparticle[k].r.y ));

						curvature = 4*Area/(length_a*length_b*length_c);
						if ( negative_sign(sceneset->scene[j].mparticle[every*(sceneset->scene[j].Nm/every)-1].r, sceneset->scene[j].mparticle[0].r, sceneset->scene[j].mparticle[every].r) )
						{
							curvature *=-1;
						}
						out_file << curvature << endl;
//						out_file << k << "\t" << sceneset->scene[j].mparticle[k].r << "\t" << curvature << endl;
					}

					else if( k >= every && k< (sceneset->scene[j].Nm-every) && k%every == 0){
						double length_a, length_b, length_c, Area, curvature = 0;
						length_a = sqrt( (sceneset->scene[j].mparticle[k+every].r - sceneset->scene[j].mparticle[k].r).Square() );

						length_b = sqrt( (sceneset->scene[j].mparticle[k-every].r - sceneset->scene[j].mparticle[k+every].r).Square() );

						length_c = sqrt( (sceneset->scene[j].mparticle[k].r - sceneset->scene[j].mparticle[k-every].r).Square() );

						Area = abs(0.5*( sceneset->scene[j].mparticle[k-every].r.x*sceneset->scene[j].mparticle[k].r.y + sceneset->scene[j].mparticle[k].r.x*sceneset->scene[j].mparticle[k+every].r.y + sceneset->scene[j].mparticle[k+every].r.x*sceneset->scene[j].mparticle[k-every].r.y - sceneset->scene[j].mparticle[k-every].r.x*sceneset->scene[j].mparticle[k+every].r.y - sceneset->scene[j].mparticle[k].r.x*sceneset->scene[j].mparticle[k-every].r.y - sceneset->scene[j].mparticle[k+every].r.x*sceneset->scene[j].mparticle[k].r.y ));

						curvature = 4*Area/(length_a*length_b*length_c);
						if ( negative_sign(sceneset->scene[j].mparticle[k-every].r, sceneset->scene[j].mparticle[k].r, sceneset->scene[j].mparticle[k+every].r) )
						{
							curvature *=-1;
						}
						out_file << curvature << endl;
//						out_file << k << "\t" << sceneset->scene[j].mparticle[k].r << "\t" << curvature << endl;
					}

					else if( k == (every*(sceneset->scene[j].Nm/every)-1)){
						double length_a, length_b, length_c, Area, curvature = 0;
						length_a = sqrt( (sceneset->scene[j].mparticle[0].r - sceneset->scene[j].mparticle[k].r).Square() );

						length_b = sqrt( (sceneset->scene[j].mparticle[k-every].r - sceneset->scene[j].mparticle[0].r).Square() );

						length_c = sqrt( (sceneset->scene[j].mparticle[k].r - sceneset->scene[j].mparticle[k-every].r).Square() );

						Area = abs(0.5*( sceneset->scene[j].mparticle[k-every].r.x*sceneset->scene[j].mparticle[k].r.y + sceneset->scene[j].mparticle[k].r.x*sceneset->scene[j].mparticle[0].r.y + sceneset->scene[j].mparticle[0].r.x*sceneset->scene[j].mparticle[k-every].r.y - sceneset->scene[j].mparticle[k-every].r.x*sceneset->scene[j].mparticle[0].r.y - sceneset->scene[j].mparticle[k].r.x*sceneset->scene[j].mparticle[k-every].r.y - sceneset->scene[j].mparticle[0].r.x*sceneset->scene[j].mparticle[k].r.y ));

						curvature = 4*Area/(length_a*length_b*length_c);
						if ( negative_sign(sceneset->scene[j].mparticle[k-every].r, sceneset->scene[j].mparticle[k].r, sceneset->scene[j].mparticle[0].r) )
						{
							curvature *=-1;
						}
						out_file << curvature << endl;
//						out_file << k << "\t" << sceneset->scene[j].mparticle[k].r << "\t" << curvature << endl;
					}
				}
				out_file << "\n" << endl;
			}

//			cout << name << endl;

			delete sceneset;
		}
		else
			cout << "Was not able to open file: " << name << endl;
	}


	return 0;
}

