#include <boost/algorithm/string.hpp>

#include<iostream>
#include<cstdlib>
#include<vector>

#include"analyze.h"

using namespace std;

int main(int argc, char** argv)
{
	SavingVector::Init_Rand(321);

	for (int i = 1; i < argc; i++)
	{
		string name = argv[i];
		SceneSet* sceneset = new SceneSet(name);
		cout << "Reading the file" << endl;
		bool read_state = sceneset->Read();
		cout << "Finished" << endl;
		cout << "Check the file" << endl;
		bool overlap = false;
		if (read_state)
		{
			for (int frame = 0; frame < sceneset->Nf; frame++)
				for (int j = 0; j < sceneset->scene[frame].Ns; j++)
					for (int k = j+1; k < sceneset->scene[frame].Ns; k++)
					{
						double d = (sceneset->scene[frame].sparticle[j].r - sceneset->scene[frame].sparticle[k].r).Square();
						if (d < 0.1)
						{
							overlap = true;
							cout << "The file is curropted." << endl;
							exit(0);
						}
					}
		}
		else
			cout << "Can not read the file" << endl;
		cout << "Finished" << endl;

		if (overlap)
			cout << "The file is curropted." << endl;
		else
			cout << "The file is healthy." << endl;

		delete sceneset;
	}

	return 0;
}

