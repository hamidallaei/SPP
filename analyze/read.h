#ifndef _READ_
#define _READ_

#include "../shared/c2dvector.h"
#include "visualparticle.h"
#include "field.h"
#include <boost/algorithm/string.hpp>

class Scene{
public:
	VisualChain* sparticle;
	VisualMembrane* mparticle;
	Real t;
	bool health_status;
	static Real density;
	static Real noise;
	int Ns;
	int Nm;
	static int chain_length;
	C2DVector L;
	Scene();
	Scene(const Scene&);
	~Scene();
	void Init(int, int);
	void Reset();
	void Draw();
	void Magnify(C2DVector r0, float d0, C2DVector r1, float d1);
	void Auto_Correlation();
	void Skip_File(std::istream& is, int n);
	friend std::istream& operator>>(std::istream& is, Scene& scene);
	friend std::ostream& operator<<(std::ostream& os, Scene& scene);
};

Real Scene::density;
Real Scene::noise;
int Scene::chain_length;

Scene::Scene()
{
	Ns = Nm = 0;
	sparticle = NULL;
	mparticle = NULL;
}

Scene::Scene(const Scene& s)
{
	health_status = s.health_status;
	t = s.t;
	chain_length = s.chain_length;
	Init(s.Ns, s.Nm);
	for (int i = 0; i < Ns; i++)
		sparticle[i] = s.sparticle[i];
	for (int i = 0; i < Nm; i++)
		mparticle[i] = s.mparticle[i];
	L = s.L;
}

Scene::~Scene()
{
	Reset();
}

void Scene::Init(int input_Ns, int input_Nm)
{
	Ns = input_Ns;
	Nm = input_Nm;
	if (Ns > 0)
		sparticle = new VisualChain[Ns];
	if (Nm > 0)
		mparticle = new VisualMembrane[Nm];
}

void Scene::Reset()
{
	if (Ns != 0)
		delete [] sparticle;
	if (Nm != 0)
		delete [] mparticle;
	Ns = Nm = 0;
}

#ifdef VISUAL
void Scene::Draw()
{
	for (int i = 0; i < Ns; i++)
		sparticle[i].Draw();

	for (int i = 0; i < Nm; i++)
		mparticle[i].Draw();

	if (Nm > 0)
		mparticle[0].Draw(0,0,1,true);

	cout << "Time is at:\t" << t << "\tR/v_0" << endl;
}

void Scene::Magnify(C2DVector r0, float d0, C2DVector r1, float d1)
{
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glLineWidth(2);
	glColor4f(1, 1, 1,1);
	glBegin(GL_POLYGON);
	glVertex2f(r1.x-d1, r1.y-d1);
	glVertex2f(r1.x+d1, r1.y-d1);
	glVertex2f(r1.x+d1, r1.y+d1);
	glVertex2f(r1.x-d1, r1.y+d1);
	glEnd();
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	glColor4f(0, 0, 0,1);
	glBegin(GL_POLYGON);
	glVertex2f(r1.x-d1, r1.y-d1);
	glVertex2f(r1.x+d1, r1.y-d1);
	glVertex2f(r1.x+d1, r1.y+d1);
	glVertex2f(r1.x-d1, r1.y+d1);
	glEnd();
	glBegin(GL_POLYGON);
	glVertex2f(r0.x-d0, r0.y-d0);
	glVertex2f(r0.x+d0, r0.y-d0);
	glVertex2f(r0.x+d0, r0.y+d0);
	glVertex2f(r0.x-d0, r0.y+d0);
	glEnd();
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	glLineWidth(1);
//	glBegin(GL_LINES);
//	glVertex2f(r1.x+d1, r1.y-d1);
//	glVertex2f(r0.x+d0, r0.y-d0);
//	glEnd();

//	glBegin(GL_LINES);
//	glVertex2f(r1.x-d1, r1.y+d1);
//	glVertex2f(r0.x-d0, r0.y+d0);
//	glEnd();

	for (int i = 0; i < Nm; i++)
		mparticle[i].Draw_Magnified(r0,d0,r1,d1);
	for (int i = 0; i < Ns; i++)
		sparticle[i].Draw_Magnified(r0,d0,r1,d1);
}
#endif

void Scene::Skip_File(std::istream& in, int n)
{
	for (int i = 0; i < n; i++)
	{
		in.read((char*) &t, sizeof(double) / sizeof(char));
		in.read((char*) &L.x, sizeof(double) / sizeof(char));
		in.read((char*) &L.y, sizeof(double) / sizeof(char));
		in.read((char*) &chain_length, sizeof(int) / sizeof(char));
		in.read((char*) &Nm, sizeof(int) / sizeof(char));
		in.read((char*) &Ns, sizeof(int) / sizeof(char));
		C2DVector temp_r;
		Real temp_theta;
		for (int j = 0; j < Nm; j++)
		{
			in >> temp_r;
		}
		for (int j = 0; j < Ns; j++)
		{
			in >> temp_r;
			Saving_Real temp_float;
			in.read((char*) &temp_float,sizeof(Saving_Real) / sizeof(char));
		}
	}
}

std::istream& operator>>(std::istream& is, Scene& scene)
{
	scene.health_status = true;
	is.read((char*) &(scene.t), sizeof(double) / sizeof(char));
	is.read((char*) &(scene.L.x), sizeof(double) / sizeof(char));
	is.read((char*) &(scene.L.y), sizeof(double) / sizeof(char));
	is.read((char*) &(scene.chain_length), sizeof(int) / sizeof(char));
	is.read((char*) &(scene.Ns), sizeof(int) / sizeof(char));
	is.read((char*) &(scene.Nm), sizeof(int) / sizeof(char));

	if (is.eof() || is.tellg() < 0)
	{
		scene.health_status = false;
		is.seekg(0,ios_base::end);
		return is;
	}

	scene.mparticle = new VisualMembrane[scene.Nm];
	scene.sparticle = new VisualChain[scene.Ns];

	Real Lx = scene.L.x;
	Real Ly = scene.L.y;
	for (int i = 0; i < scene.Nm; i++)
	{
		is >> scene.mparticle[i].r;
		#ifdef Periodic_Show
			scene.mparticle[i].r.x -= Lx*2*((int) floor(scene.mparticle[i].r.x / (Lx*2) + 0.5));
			scene.mparticle[i].r.y -= Ly*2*((int) floor(scene.mparticle[i].r.y / (Ly*2) + 0.5));
		#endif

		if (i == scene.Nm - 1 && (is.eof() || is.tellg() < 0))
		{
			scene.health_status = false;
			is.seekg(0,ios_base::end);
			return is;
		}
	}
	for (int i = 0; i < scene.Ns; i++)
	{
		is >> scene.sparticle[i].r;
		#ifdef Periodic_Show
			scene.sparticle[i].r.x -= Lx*2*((int) floor(scene.sparticle[i].r.x / (Lx*2) + 0.5));
			scene.sparticle[i].r.y -= Ly*2*((int) floor(scene.sparticle[i].r.y / (Ly*2) + 0.5));
		#endif
		Saving_Real temp_float;
		is.read((char*) &temp_float,sizeof(Saving_Real) / sizeof(char));
		scene.sparticle[i].theta = temp_float;

		if (i == scene.Ns - 1 && (is.eof() || is.tellg() < 0))
		{
			scene.health_status = false;
			is.seekg(0,ios_base::end);
			return is;
		}
	}

	if (scene.health_status)
		VisualChain::chain_length = scene.chain_length;
}

std::ostream& operator<<(std::ostream& os, Scene& scene)
{
	os.write((char*) &(scene.t), sizeof(double) / sizeof(char));
	os.write((char*) &(scene.L.x), sizeof(double) / sizeof(char));
	os.write((char*) &(scene.L.y), sizeof(double) / sizeof(char));
	os.write((char*) &(scene.chain_length), sizeof(int) / sizeof(char));
	os.write((char*) &(scene.Ns), sizeof(int) / sizeof(char));
	os.write((char*) &(scene.Nm), sizeof(int) / sizeof(char));
	for (int i = 0; i < scene.Nm; i++)
		scene.mparticle[i].Write(os);
	for (int i = 0; i < scene.Ns; i++)
		scene.sparticle[i].Write(os);
}

class SceneSet{
public:
	vector<Scene> scene;
	C2DVector L, L_min;
	stringstream address;
	string info;
	ifstream input_file;
	ofstream output_file;
	Field* field;

	SceneSet(string input_address);
	~SceneSet();
	bool Read(int skip = 0);
	void Write(int, int); // write from a time to the end of the file
	void Write(int start, int end, int limit); // write from start to the end if end-start > limit
	void Write_Every(int interval); // write every interval
	void Save_Theta_Deviation(int, int, int, string);
	void Plot_Fields(int, int, string);
	void Plot_Averaged_Fields(int smaller_grid_dim, string name);
	void Plot_Averaged_Fields_Section(int smaller_grid_dim, int y, string info);
	void Plot_Density_Contour(int smaller_grid_dim, double rho, string info);
	void Accumulate_Theta(int smaller_grid_dim, const int num_bins, const double& p_c, const double& dp, const string& info);
};

SceneSet::SceneSet(string input_address)
{
	L.Null();
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
	L.Null();

	input_file.open(address.str().c_str());
	if (input_file.is_open())
	{
		input_file.seekg(0,ios_base::end);
		int end_of_file = input_file.tellg();
		input_file.seekg(0,ios_base::beg);

		if (skip > 0)
			temp_scene.Skip_File(input_file, skip);
		temp_scene.Reset();

		while (input_file.tellg() < end_of_file && input_file.tellg() >= 0)
		{
			scene.push_back(temp_scene);
			int index = scene.size()-1;
			input_file >> scene[index];

			for (int i = 0; i < scene[index].Nm; i++)
			{
				if (abs(scene[index].mparticle[i].r.x) > L.x)
					L.x = abs(scene[index].mparticle[i].r.x);
				if (abs(scene[index].mparticle[i].r.y) > L.y)
					L.y = abs(scene[index].mparticle[i].r.y);
			}
			for (int i = 0; i < scene[index].Ns; i++)
			{
				if (abs(scene[index].sparticle[i].r.x) > L.x)
					L.x = abs(scene[index].sparticle[i].r.x);
				if (abs(scene[index].sparticle[i].r.y) > L.y)
					L.y = abs(scene[index].sparticle[i].r.y);
			}
			if (0.000001*counter*( scene[index].Ns*sizeof(VisualChain) + scene[index].Nm*sizeof(VisualMembrane) + sizeof(scene[index]) ) > (600))
			{
				cout << "File is too big" << endl;
				return(false);
			}
			if (counter % 1 != 0)
				scene.pop_back();
			counter++;
		}
	}
	else
		return (false);
	input_file.close();

	if (!scene[scene.size()-1].health_status)
		scene.pop_back();

	L.x = round(L.x+0.1);
	L.y = round(L.y+0.1);
	L_min = L;
	L.x += 0.5;
	L.y += 0.5;

	L = scene[0].L;
	L_min = L;

	return (true);
}

void SceneSet::Write(int start, int limit)
{
	if ((scene.size() - start) > limit)
	{
		output_file.open(address.str().c_str());
		for (int i = start; i < scene.size(); i++)
			if (scene[i].health_status)
				output_file << scene[i];
		output_file.close();
	}
	else
		cout << "I will not cut the file because it is short enough!" << endl;
}

void SceneSet::Write(int start, int end, int limit)
{
	if ((end - start) > limit)
	{
		output_file.open(address.str().c_str());
		for (int i = start; i < end; i++)
			if (scene[i].health_status)
				output_file << scene[i];
		output_file.close();
	}
	else
		cout << "I will not cut the file because it is short enough!" << endl;
}

void SceneSet::Write_Every(int interval)
{
	output_file.open(address.str().c_str());
	for (int i = 0; i < scene.size(); i += interval)
	{
		if (scene[i].health_status)
		{
			output_file << scene[i];
		}
	}
	output_file.close();
}

void SceneSet::Save_Theta_Deviation(int smaller_grid_dim, int start_t, int end_t, string info)
{
	Field averaged_field(smaller_grid_dim, L);
	for (int i = start_t; i < scene.size() && i < end_t; i++)
	{
		Field f(smaller_grid_dim, L);
		f.Compute(scene[i].sparticle, scene[i].Ns);
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

void SceneSet::Plot_Fields(int smaller_grid_dim, int t, string name)
{
	field = new Field(smaller_grid_dim,L);
	field->Compute((BasicParticle0*) scene[t].sparticle, scene[t].Ns);
	field->Draw(name);
}

void SceneSet::Plot_Averaged_Fields(int smaller_grid_dim, string info)
{
	Field averaged_field(smaller_grid_dim, L);
	for (int i = 0; i < scene.size(); i++)
	{
		Field f(smaller_grid_dim, L);
		f.Compute(scene[i].sparticle, scene[i].Ns);
		averaged_field.Add(&f);
	}
	averaged_field.Average();
	averaged_field.Draw(info);
	averaged_field.Reset();
}

void SceneSet::Plot_Averaged_Fields_Section(int smaller_grid_dim, int y, string info)
{
	Field averaged_field(smaller_grid_dim, L);
	for (int i = 0; i < scene.size(); i++)
	{
		Field f(smaller_grid_dim, L);
		f.Compute(scene[i].sparticle, scene[i].Ns);
		averaged_field.Add(&f);
	}
	averaged_field.Average();
	averaged_field.Draw_Section(info, y);
	averaged_field.Reset();
}

void SceneSet::Plot_Density_Contour(int smaller_grid_dim, double rho, string info)
{
	Field averaged_field(smaller_grid_dim, L);
	Field f(smaller_grid_dim, L);
	for (int i = 0; i < scene.size(); i++)
	{
		f.Compute(scene[i].sparticle, scene[i].Ns);
		averaged_field.Add(&f);
		f.Reset();
	}
	averaged_field.Average();
	averaged_field.Draw_Density_Contour(info, rho);
}

void SceneSet::Accumulate_Theta(int smaller_grid_dim, const int num_bins, const double& p_c, const double& dp, const string& info)
{
	Field f(smaller_grid_dim, L);
	Stat<double> dtheta;
	int step = (scene.size()/500);
	if (step == 0)
		step = 1;
	for (int i = 0; i < scene.size(); i+=step)
	{
		f.Compute(scene[i].sparticle, scene[i].Ns);
		f.Add_Theta(dtheta,p_c,dp);
		f.Reset();
	}
	stringstream address;
	address.str("");
	address << info << "-p-" << p_c << "-dp-" << dp << ".dat";
	dtheta.Periodic_Transform(M_PI);
	dtheta.Shift_Average();
	dtheta.Periodic_Transform(M_PI);
	dtheta.Shift_Average();
	dtheta.Compute();
	dtheta.Histogram(num_bins,address.str().c_str());
	cout << "Number of date: " << dtheta.data.size() << endl;
	dtheta.Reset();
}

#endif
