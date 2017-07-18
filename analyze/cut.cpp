#include <boost/algorithm/string.hpp>

#include<iostream>
#include<cstdlib>
#include<vector>

#include"analyze.h"

using namespace std;

void Cut_Every(string& name, int jump)
{
		SceneSet* sceneset = new SceneSet(name);
		cout << "Reading the file" << endl;
		bool read_state = sceneset->Read();
		cout << "Finished" << endl;
		cout << "Write the file" << endl;
		if (read_state)
		{
			sceneset->Write_Every(jump);
		}
		else
			cout << "Can not read the file" << endl;
		cout << "Finished" << endl;

		delete sceneset;
}

void Cut_To(string& name, float final_time)
{
		SceneSet* sceneset = new SceneSet(name);
		cout << "Reading the file" << endl;
		bool read_state = sceneset->Read();
		cout << "Finished" << endl;
		cout << "Write the file" << endl;
		if (read_state)
		{
			sceneset->Write(0.0, final_time);
		}
		else
			cout << "Can not read the file" << endl;
		cout << "Finished" << endl;

		delete sceneset;
}

void Cut_From(string& name, float start_time)
{
		SceneSet* sceneset = new SceneSet(name);
		cout << "Reading the file" << endl;
		bool read_state = sceneset->Read();
		cout << "Finished" << endl;
		cout << "Write the file" << endl;
		if (read_state)
		{
			sceneset->Write(start_time, sceneset->scene[sceneset->Nf-1].t);
		}
		else
			cout << "Can not read the file" << endl;
		cout << "Finished" << endl;

		delete sceneset;
}

int main(int argc, char** argv)
{
	SavingVector::Init_Rand(321);

	string name = argv[1];
	float input_number = atof(argv[2]);
	Cut_To(name, input_number);

	return 0;
}

