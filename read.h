#ifndef _READ_
#define _READ_

#include "c2dvector.h"
#include "visualparticle.h"
#include "field.h"

class Scene{
public:
	VisualParticle* particle;
	static float L;
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
};

float Scene::L;
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

class SceneSet{
public:
	vector<Scene> scene;
	Real L;
	stringstream address;
	ifstream input_file;
	Field* field;

	SceneSet(string input_address);
	~SceneSet();
	void Read();
	void Plot_Fields(int);
	void Plot_Averaged_Fields(int grid_dim, string name);
};

SceneSet::SceneSet(string input_address) : L(0)
{
	address.str("");
	address << input_address;
	input_file.open(address.str().c_str());
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
//		if (counter % 100 == 0)
//			cout << counter << endl;
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
	L = round(1000*L)/1000;
	L += 0.5;
}

void SceneSet::Plot_Fields(int t)
{
	field = new Field(20,L);
	field->Compute((BasicParticle*) scene[t].particle, Scene::number_of_particles);
	field->Draw("Hello");
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

#endif
