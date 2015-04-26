#ifndef _READ_
#define _READ_

#include "../shared/c2dvector.h"
#include "visualparticle.h"
#include "field.h"
#include <boost/algorithm/string.hpp>

class Scene{
public:
	VisualParticle* particle;
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
}

void Scene::Draw()
{
	glColor3f(VisualParticle::color.red, VisualParticle::color.green, VisualParticle::color.blue);
	for (int i = 0; i < number_of_particles; i++)
		particle[i].Draw();
}

void Scene::Skip_File(std::istream& in, int n)
{
	for (int i = 0; i < n; i++)
	{
		in.read((char*) &number_of_particles, sizeof(int) / sizeof(char));
		C2DVector temp_r, temp_v;
		for (int j = 0; j < number_of_particles; j++)
		{
			in >> temp_r;
			in >> temp_v;
		}
	}
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
	Real L_min;
	stringstream address;
	string info;
	ifstream input_file;
	ofstream output_file;
	Field* field;

	SceneSet(string input_address);
	~SceneSet();
	bool Read(int skip = 0);
	void Write(int, int); // write from a time to the end
	void Save_Theta_Deviation(int, int, int, string);
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

bool SceneSet::Read(int skip)
{
	int counter = 0;
	static Scene temp_scene;
	L = 0;

	input_file.open(address.str().c_str());
	if (input_file.is_open())
	{
		if (skip > 0)
			temp_scene.Skip_File(input_file, skip);

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
			if (counter*sizeof(temp_scene) > (3000000000))
				return(false);
		}
	}
	else
		return (false);
	input_file.close();

	scene.pop_back();

	L = round(L+0.1);
	L_min = L;
	L += 0.5;

	return (true);
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

void SceneSet::Save_Theta_Deviation(int grid_dim, int start_t, int end_t, string info)
{
	Field averaged_field(grid_dim, L);
	for (int i = start_t; i < scene.size() && i < end_t; i++)
	{
		Field f(grid_dim, L);
		f.Compute(scene[i].particle, Scene::number_of_particles);
		averaged_field.Add(&f);
	}
	averaged_field.Average();
	stringstream address("");
	address << info;
	ofstream outfile(address.str().c_str());
	averaged_field.Save_Theta_Deviation(outfile);
	outfile.close();
	averaged_field.Reset();
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
