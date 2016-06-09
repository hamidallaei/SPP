#include <boost/algorithm/string.hpp>

#include<iostream>
#include<cstdlib>
#include<vector>

#include"analyze.h"

using namespace std;

int main(int argc, char** argv)
{
	C2DVector::Init_Rand(321);

	for (int i = 1; i < argc; i++)
	{
		string name = argv[i];
		SceneSet* sceneset = new SceneSet(name);
		cout << "Reading the file" << endl;
		bool read_state = sceneset->Read(1);
		cout << "Finished" << endl;
		cout << "Write the file" << endl;
		if (read_state)
			sceneset->Write(0, sceneset->scene.size()-1,0);
		else
			cout << "Can not read the file" << endl;
		cout << "Finished" << endl;

		delete sceneset;
	}

	return 0;
}

