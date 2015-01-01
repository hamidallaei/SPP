#ifndef _READ_
#define _READ_

#include "../shared/c2dvector.h"
#include "visualparticle.h"
#include "field.h"
#include <boost/algorithm/string.hpp>

class Scene{
public:
	VisualParticle* particle;
	static float L;
	static float density;
	static float noise;
	static int number_of_particles;
	Scene();
	Scene(const Scene&);
	~Scene();
	void Init(int);
	void Reset();
	void Draw();
	void Auto_Correlation();
	void Skip_File(std::istream& is, int n);
	friend std::istream& operator>>(std::istream& is, Scene& scene);
	friend std::ostream& operator<<(std::ostream& os, Scene& scene);
};

float Scene::L;
float Scene::density;
float Scene::noise;
int Scene::number_of_particles;

Scene::Scene()
{
	particle = NULL;
}

Scene::Scene(const Scene& s)
{
	particle = s.particle;
}

Scene::~Scene()
{
}

void Scene::Init(int num)
{
	number_of_particles = num;
	particle = new VisualParticle[number_of_particles];
}

void Scene::Reset()
{
	delete [] particle;
	particle = NULL;
}

void Scene::Draw()
{
	glColor3f(VisualParticle::color.red, VisualParticle::color.green, VisualParticle::color.blue);
	for (int i = 0; i < number_of_particles; i++)
		particle[i].Draw();
}

void Scene::Skip_File(std::istream& in, int n)
{
	in.read((char*) &number_of_particles, sizeof(int) / sizeof(char));
	for (int i = 0; i < number_of_particles; i++)
	{
		in >> particle[i].r;
		in >> particle[i].v;
	}
	in.ignore(numeric_limits<streamsize>::max(), '\n');
	for (int i = 0; i < (n*(number_of_particles+1)); i++)
		in.ignore(numeric_limits<streamsize>::max(), '\n');
}

std::istream& operator>>(std::istream& is, Scene& scene)
{
	is.read((char*) &(scene.number_of_particles), sizeof(int) / sizeof(char));
	scene.particle = new VisualParticle[scene.number_of_particles];
	for (int i = 0; i < scene.number_of_particles; i++)
	{
		is >> scene.particle[i].r;
		is >> scene.particle[i].v;
	}
}

std::ostream& operator<<(std::ostream& os, Scene& scene)
{
	os.write((char*) &(scene.number_of_particles), sizeof(int) / sizeof(char));
	for (int i = 0; i < scene.number_of_particles; i++)
	{
		scene.particle[i].r.write(os);
		scene.particle[i].v.write(os);
	}
}

class SceneSet{
public:
	vector<Scene> scene;
	Real L;
	stringstream address;
	string info;
	ifstream input_file;
	ofstream output_file;
	Field* field;

	SceneSet(string input_address);
	~SceneSet();
	void Read();
	void Write(int, int); // write from a time to the end
	void Plot_Fields(int, int, string);
	void Plot_Averaged_Fields(int grid_dim, string name);
	void Plot_Averaged_Fields_Section(int grid_dim, int y, string info);
	void Plot_Density_Contour(int grid_dim, double rho, string info);
};

SceneSet::SceneSet(string input_address) : L(0)
{
	address.str("");
	address << input_address;
	string name = input_address;
	boost::replace_all(name, "-r-v.bin", "");
	info = name;
	boost::replace_all(name, "rho=", "");
	boost::replace_all(name, "-noise=", "\t");

	stringstream ss;
	ss.str("");
	ss << name;
	ss >> Scene::density;
	ss >> Scene::noise;
}

SceneSet::~SceneSet()
{
	for (int i = 0; i < scene.size(); i++)
		scene[i].Reset();
	scene.clear();
}

void SceneSet::Read()
{
	int counter = 0;
	static Scene temp_scene;
	L = 0;

//	temp_scene.Skip_File(input_file, 3800);
	input_file.open(address.str().c_str());

	while (!input_file.eof())
	{
		counter++;
		input_file >> temp_scene;
		scene.push_back(temp_scene);
		for (int i = 0; i < Scene::number_of_particles; i++)
		{
			if (abs(temp_scene.particle[i].r.x) > L)
				L = abs(temp_scene.particle[i].r.x);
			if (abs(temp_scene.particle[i].r.y) > L)
				L = abs(temp_scene.particle[i].r.y);
		}
	}
	input_file.close();

	scene.pop_back();

//	for (int i = 0; i < (counter-1); i++)
//	{
//		for (int j = 0; j < scene[i].number_of_particles; j++)
//		{
//			scene[i].particle[j].v = (scene[i+1].particle[j].r - scene[i].particle[j].r);
//			scene[i].particle[j].v.Periodic_Transform();
//		}
//	}



	L = round(L+0.1);
	L += 0.5;
}

void SceneSet::Write(int start, int limit)
{
	if ((scene.size() - start) > limit)
	{
		output_file.open(address.str().c_str());
		for (int i = start; i < scene.size(); i++)
			output_file << scene[i];
		output_file.close();
	}
	else
		cout << "I will not cut the file because it is short enough!" << endl;
}

void SceneSet::Plot_Fields(int grid_dim, int t, string name)
{
	field = new Field(grid_dim,L);
	field->Compute((BasicParticle*) scene[t].particle, Scene::number_of_particles);
	field->Draw(name);
}

void SceneSet::Plot_Averaged_Fields(int grid_dim, string info)
{
	Field averaged_field(grid_dim, L);
	for (int i = 0; i < scene.size(); i++)
	{
		Field f(grid_dim, L);
		f.Compute(scene[i].particle, Scene::number_of_particles);
		averaged_field.Add(&f);
	}
	averaged_field.Average();
	averaged_field.Draw(info);
	averaged_field.Reset();
}

void SceneSet::Plot_Averaged_Fields_Section(int grid_dim, int y, string info)
{
	Field averaged_field(grid_dim, L);
	for (int i = 0; i < scene.size(); i++)
	{
		Field f(grid_dim, L);
		f.Compute(scene[i].particle, Scene::number_of_particles);
		averaged_field.Add(&f);
	}
	averaged_field.Average();
	averaged_field.Draw_Section(info, y);
	averaged_field.Reset();
}

void SceneSet::Plot_Density_Contour(int grid_dim, double rho, string info)
{
	Field averaged_field(grid_dim, L);
	Field f(grid_dim, L);
	for (int i = 0; i < scene.size(); i++)
	{
		f.Compute(scene[i].particle, Scene::number_of_particles);
		averaged_field.Add(&f);
		f.Reset();
	}
	averaged_field.Average();
	averaged_field.Draw_Density_Contour(info, rho);
}

#endif
