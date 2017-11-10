#ifndef _READ_
#define _READ_

#include "../shared/c2dvector.h"
#include "visualparticle.h"
#include "field.h"
#include <boost/algorithm/string.hpp>

const float mem_max = 1000;

class Scene{
public:
	VisualChain* sparticle;
	VisualMembrane* mparticle;
	Real t;
	Real omega_s;
	bool health_status;
	static Real density;
	static Real noise;
	static int Ns;
	static int Nm;
	static int chain_length;
	SavingVector L;
	SavingVector r_cm;
	Scene();
	Scene(const Scene&);
	~Scene();
	void Init(int, int);
	void Delete();
	void Periodic_Transform();
	void Translate_Rcm();
	void Draw();
	void Draw_Translated_Scaled();
	void Magnify(SavingVector r0, float d0, SavingVector r1, float d1);
	void Auto_Correlation();
	void Skip_File(std::istream& is, int n);
	friend std::istream& operator>>(std::istream& is, Scene& scene);
	friend std::ostream& operator<<(std::ostream& os, Scene& scene);
};

Real Scene::density;
Real Scene::noise;
int Scene::chain_length;
int Scene::Ns;
int Scene::Nm;

Scene::Scene()
{
	if (Ns > 0)
		sparticle = new VisualChain[Ns];
	if (Nm > 0)
		mparticle = new VisualMembrane[Nm];
}

Scene::Scene(const Scene& s)
{
	health_status = s.health_status;
	t = s.t;
	chain_length = s.chain_length;
	if (Ns > 0)
		sparticle = new VisualChain[Ns];
	if (Nm > 0)
		mparticle = new VisualMembrane[Nm];
	for (int i = 0; i < Ns; i++)
		sparticle[i] = s.sparticle[i];
	for (int i = 0; i < Nm; i++)
		mparticle[i] = s.mparticle[i];
	L = s.L;
}

Scene::~Scene()
{
	Delete();
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

void Scene::Delete()
{
	if (Ns != 0)
		delete [] sparticle;
	if (Nm != 0)
		delete [] mparticle;
}

void Scene::Periodic_Transform()
{
	for (int i = 0; i < Nm; i++)
	{
		mparticle[i].r.x -= Lx*2*((int) floor(mparticle[i].r.x / (Lx*2) + 0.5));
		mparticle[i].r.y -= Ly*2*((int) floor(mparticle[i].r.y / (Ly*2) + 0.5));
	}
	for (int i = 0; i < Ns; i++)
	{
		sparticle[i].r.x -= Lx*2*((int) floor(sparticle[i].r.x / (Lx*2) + 0.5));
		sparticle[i].r.y -= Ly*2*((int) floor(sparticle[i].r.y / (Ly*2) + 0.5));
	}
}

void Scene::Translate_Rcm()
{
		for (int i = 0; i < Nm; i++)
			mparticle[i].r -= r_cm;
		for (int i = 0; i < Ns; i++)
			sparticle[i].r -= r_cm;
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
}

void Scene::Draw_Translated_Scaled()
{
	for (int i = 0; i < Ns; i++)
		sparticle[i].Draw(r_cm.x, r_cm.y);

	for (int i = 0; i < Nm; i++)
		mparticle[i].Draw(r_cm.x, r_cm.y);

	if (Nm > 0)
		mparticle[0].Draw(r_cm.x, r_cm.y, 1,true);
}

void Scene::Magnify(SavingVector r0, float d0, SavingVector r1, float d1)
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
		mparticle[i].Draw_Magnified(r0,d0,r1,d1, i == 0);
	for (int i = 0; i < Ns; i++)
		sparticle[i].Draw_Magnified(r0,d0,r1,d1);
}
#endif

void Scene::Skip_File(std::istream& in, int n)
{
	for (int i = 0; i < n; i++)
	{
		double Lx, Ly;
		in.read((char*) &t, sizeof(double) / sizeof(char));
		in.read((char*) &Lx, sizeof(double) / sizeof(char));
		in.read((char*) &Ly, sizeof(double) / sizeof(char));
		in.read((char*) &chain_length, sizeof(int) / sizeof(char));
		in.read((char*) &Nm, sizeof(int) / sizeof(char));
		in.read((char*) &Ns, sizeof(int) / sizeof(char));

		L.x = Lx;
		L.y = Ly;

		SavingVector temp_r;
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
	int temp_int;
	double Lx, Ly;
	is.read((char*) &(scene.t), sizeof(double) / sizeof(char));
	is.read((char*) &(Lx), sizeof(double) / sizeof(char));
	is.read((char*) &(Ly), sizeof(double) / sizeof(char));
	is.read((char*) &(scene.chain_length), sizeof(int) / sizeof(char));
	is.read((char*) &(temp_int), sizeof(int) / sizeof(char));
	is.read((char*) &(temp_int), sizeof(int) / sizeof(char));

	scene.L.x = Lx;
	scene.L.y = Ly;

	if (is.eof() || is.tellg() < 0)
	{
		scene.health_status = false;
		is.seekg(0,ios_base::end);
		return is;
	}

	scene.r_cm.Null();
	for (int i = 0; i < scene.Nm; i++)
	{
		is >> scene.mparticle[i].r;
		scene.r_cm += scene.mparticle[i].r;
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
		scene.r_cm += scene.sparticle[i].r;
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

	scene.r_cm /= scene.Nm + scene.Ns;

	if (scene.t < 0 || scene.t > 1000000)
		scene.health_status = false;
	if (temp_int < 0 || temp_int > 1000000)
		scene.health_status = false;

	if (scene.health_status)
		VisualChain::chain_length = scene.chain_length;
}

std::ostream& operator<<(std::ostream& os, Scene& scene)
{
	double Lx = scene.L.x;
	double Ly = scene.L.y;
	os.write((char*) &(scene.t), sizeof(double) / sizeof(char));
	os.write((char*) &(Lx), sizeof(double) / sizeof(char));
	os.write((char*) &(Ly), sizeof(double) / sizeof(char));
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
	Scene* scene;
	int Nf;
	SavingVector L, L_min;
	stringstream address;
	string info;
	ifstream input_file;
	ofstream output_file;
	Field* field;

	SceneSet(string input_address);
	~SceneSet();
	void Reset();
	int Count_Frames();
	bool Read(int skip = 0);
	void Write(int start, int end); // write from frame start to the frame end
	void Write(float, float); // write from an specefic time to another time
	void Write_Every(int interval); // write every interval
	void Save_Theta_Deviation(int, int, int, string);
	void Periodic_Transform();
	void Translate_Rcm();
	void Draw_Path(int frame);
	void Plot_Fields(int, int, string);
	void Plot_Averaged_Fields(int smaller_grid_dim, string name);
	void Plot_Averaged_Fields_Section(int smaller_grid_dim, int y, string info);
	void Plot_Density_Contour(int smaller_grid_dim, double rho, string info);
	void Accumulate_Theta(int smaller_grid_dim, const int num_bins, const double& p_c, const double& dp, const string& info);
	void Compute_Omega();
};

SceneSet::SceneSet(string input_address)
{
	Nf = 0;
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
	Reset();
}

void SceneSet::Reset()
{
	if (Nf > 0)
		delete [] scene;
	Nf = 0;
	Scene::Ns = 0;
	Scene::Nm = 0;
}

int SceneSet::Count_Frames()
{
	int counter = 0;

	input_file.open(address.str().c_str());
	if (input_file.is_open())
	{
		input_file.seekg(0,ios_base::end);
		int end_of_file = input_file.tellg();
		input_file.seekg(0,ios_base::beg);

		double temp_double;
		int temp_int;
		input_file.read((char*) &(temp_double), sizeof(double) / sizeof(char));
		input_file.read((char*) &(temp_double), sizeof(double) / sizeof(char));
		input_file.read((char*) &(temp_double), sizeof(double) / sizeof(char));
		input_file.read((char*) &(temp_int), sizeof(int) / sizeof(char));
		input_file.read((char*) &(Scene::Ns), sizeof(int) / sizeof(char));
		input_file.read((char*) &(Scene::Nm), sizeof(int) / sizeof(char));
		input_file.seekg(0,ios_base::beg);

		Scene temp_scene;

		while (input_file.tellg() < end_of_file && input_file.tellg() >= 0)
		{
			input_file >> temp_scene;
			if (temp_scene.health_status)
				counter++;
		}

		float mem_usage = 0.000001*counter*( Scene::Ns*sizeof(VisualChain) + Scene::Nm*sizeof(VisualMembrane) + sizeof(temp_scene) );
		cout << "Memory usage:\t" << mem_usage << endl;
		if (mem_usage > (mem_max))
		{
				cout << "File is too big" << endl;
				return(-2);
		}
	}
	else
		return (-1);
	input_file.close();

	return (counter);
}

bool SceneSet::Read(int skip)
{
	int index = 0;

	Nf = Count_Frames() - skip;
	scene = new Scene[Nf];

	Scene temp_scene;
	L.Null();

	input_file.open(address.str().c_str());
	if (input_file.is_open())
	{
		input_file.seekg(0,ios_base::end);
		int end_of_file = input_file.tellg();
		input_file.seekg(0,ios_base::beg);

		if (skip > 0)
			temp_scene.Skip_File(input_file, skip);
	
		while (input_file.tellg() < end_of_file && input_file.tellg() >= 0 && index < Nf)
		{
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
			float mem_usage = 0.000001*index*( Scene::Ns*sizeof(VisualChain) + Scene::Nm*sizeof(VisualMembrane) + sizeof(temp_scene) );
			if (mem_usage > (mem_max))
			{
				cout << "File is too big" << endl;
				return(false);
			}
			index++;
		}
	}
	else
		return (false);
	input_file.close();

	Nf--;
	if (!scene[Nf-1].health_status)
		Nf--;

	L.x = round(L.x+0.1);
	L.y = round(L.y+0.1);
	L.x = max(L.x, L.y);
	L.y = L.x;
	L_min = L;
	L.x += 0.5;
	L.y += 0.5;

	#ifdef Periodic_Show
		L = scene[0].L;
		L_min = L;
	#endif

	Compute_Omega();

	#ifdef Periodic_Show
		Periodic_Transform();
	#endif
	#ifdef Translated_Scaled
		Translate_Rcm();
		L.x = L.y = 0;
		for (int index = 0; index < Nf; index++)
			for (int i = 0; i < scene[index].Nm; i++)
			{
				if (abs(scene[index].mparticle[i].r.x) > L.x)
					L.x = abs(scene[index].mparticle[i].r.x);
				if (abs(scene[index].mparticle[i].r.y) > L.y)
					L.y = abs(scene[index].mparticle[i].r.y);
			}
	#endif

	return (true);
}

void SceneSet::Write(int start, int end)
{
	output_file.open(address.str().c_str());
	for (int i = start; i < end; i++)
		if (scene[i].health_status)
			output_file << scene[i];
	output_file.close();
}

void SceneSet::Write(float start, float end)
{
	output_file.open(address.str().c_str());
	for (int i = 0; i < Nf; i++)
		if (scene[i].health_status && (scene[i].t <= end) && (scene[i].t >= start))
			output_file << scene[i];
	output_file.close();
}

void SceneSet::Write_Every(int interval)
{
	output_file.open(address.str().c_str());
	for (int i = 0; i < Nf; i += interval)
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
	for (int i = start_t; i < Nf && i < end_t; i++)
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

void SceneSet::Periodic_Transform()
{
	for (int i = 0; i < Nf; i++)
		scene[i].Periodic_Transform();
}

void SceneSet::Translate_Rcm()
{
	for (int i = 0; i < Nf; i++)
		scene[i].Translate_Rcm();
}

void SceneSet::Draw_Path(int frame)
{
	glLineWidth(1);
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	glColor4f(0, 0, 0,1);
	glBegin(GL_LINES);
	for (int i = 0; i < frame; i++)
		glVertex2f(scene[i].r_cm.x, scene[i].r_cm.y);
	glEnd();
	glLineWidth(1);
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
	for (int i = 0; i < Nf; i++)
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
	for (int i = 0; i < Nf; i++)
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
	for (int i = 0; i < Nf; i++)
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
	int step = (Nf/500);
	if (step == 0)
		step = 1;
	for (int i = 0; i < Nf; i+=step)
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

void SceneSet::Compute_Omega()
{
	SavingVector v, r;
	Real Ls;
	Real Sxx, Syy, Sxy;
	for (int step = 1; step < Nf; step++)
	{
		Ls = Sxx = Sxy = Syy = 0;
		for (int i = 0; i < scene[step].Ns; i++)
		{
			v.x = (scene[step].sparticle[i].r.x - scene[step - 1].sparticle[i].r.x) / (scene[step].t - scene[step-1].t);
			v.y = (scene[step].sparticle[i].r.y - scene[step - 1].sparticle[i].r.y) / (scene[step].t - scene[step-1].t);
			r = scene[step].sparticle[i].r - scene[step].r_cm;
			Ls += r.x * v.y - r.y * v.x;
			Sxx += r.x*r.x;
			Sxy += r.x*r.y;
			Syy += r.y*r.y;
		}
		scene[step].omega_s = Ls / (Sxx + Syy);
	}
}

#endif
