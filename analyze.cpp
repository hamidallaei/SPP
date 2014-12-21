#include<iostream>
#include<cstdlib>
#include<vector>

#include"analyze.h"
#include"shapes.h"

using namespace std;

vector<Analytic_Scene> box_set;

int t = 0;
float box_dim = 0;

float Readfile(ifstream& input_file)
{
	int counter = 0;
	static Scene temp_scene;
	float L = 0;

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
		if (counter % 100 == 0)
			cout << counter << endl;
	}
	input_file.close();

	scene.pop_back();

	for (int i = 0; i < (counter-1); i++)
	{
		for (int j = 0; j < scene[i].number_of_particles; j++)
			scene[i].particle[j].v = (scene[i+1].particle[j].r - scene[i].particle[j].r);
	}
	L = round(1000*L)/1000;
	return (L+1);
}

void Init()
{
	glClearColor (1.0, 1.0, 1.0, 0.0);
	glShadeModel (GL_FLAT);
	Init_Circle();
	float sigma = 0.5;
	Particle::radius = sigma / 2;
	Particle::tail = 1*Particle::radius;
	Particle::color.red = 0.8;
	Particle::color.green = 0.0;
	Particle::color.blue = 0.0;
}

void Display()
{
	glClear(GL_COLOR_BUFFER_BIT);
	//glPushMatrix();
//	glPopMatrix();

	scene[t].Draw();

	if (save && (t < (scene.size()-1)))
	{
		static cv::Mat img(window_width,window_height,CV_8UC3);
		glPixelStorei(GL_PACK_ALIGNMENT, (img.step & 3) ? 1 : 4);
		glPixelStorei(GL_PACK_ROW_LENGTH, img.step/img.elemSize());
		glReadPixels(0, 0, img.cols, img.rows, GL_BGR, GL_UNSIGNED_BYTE, img.data);
//		cv::flip(img, flipped, 0);
		writer << img;
	}

	glutSwapBuffers();
}


void Next_Frame(void)
{
	t++;
	if (t < 0)
		t = 0;
	if (t >= scene.size())
		t = scene.size() - 1;
	glutPostRedisplay();
}

void SpecialInput(int key, int x, int y)
{
	switch(key)
	{
		case GLUT_KEY_PAGE_UP:
			t += 99;
		break;
		case GLUT_KEY_PAGE_DOWN:
			t -= 101;
		break;
		case GLUT_KEY_UP:
			t += 9;
		break;	
		case GLUT_KEY_DOWN:
			t -= 11;
		break;
		case GLUT_KEY_LEFT:
			t -= 2;
		break;
		case GLUT_KEY_RIGHT:
			t += 0;
		break;
		default:
		break;
	}
	Next_Frame();
}

void KeyboardInput(unsigned char key, int x, int y)
{
	static bool stop = true;
	if ((key == 32) || (key = 112))
		stop = !stop;
	if (stop)
		glutIdleFunc(NULL);
	else
		glutIdleFunc(Next_Frame);
}


void Reshape(int w, int h)
{
	window_width = w;
	window_height = h;

	glViewport (0, 0, (GLsizei) w, (GLsizei) h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(-box_dim, box_dim,  -box_dim, box_dim, -1.0, 1.0);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}


int main(int argc, char** argv)
{
	int glargc = 0;
	char** glargv;

	if (argc > 2)
		save = true;

	stringstream address("");
	address << argv[argc-1];
	ifstream data_file;
	data_file.open(address.str().c_str());
	box_dim = Readfile(data_file);
	data_file.close();

	string name = address.str().c_str();
	string::size_type position_of_txt = name.find("-r-v", 0);
	name.erase(position_of_txt);
	address.str("");
	address << name << ".mpg";
//	address << name << ".mpg";
	if (save)
	{
//		writer.open(address.str().c_str(), CV_FOURCC('F','L','V','1'), 20, cv::Size(window_width,window_height), true);
		writer.open(address.str().c_str(), CV_FOURCC('P','I','M','1'), 20, cv::Size(window_width,window_height), true);
	}

	glutInit(&glargc, glargv);
	glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB);
	glutInitWindowSize (window_width, window_height);
	glutInitWindowPosition (100, 100);
	glutCreateWindow ("Created By Hamid!");
	Init ();
	glutDisplayFunc(Display);
	glutReshapeFunc(Reshape);
	glutSpecialFunc(SpecialInput);
	glutKeyboardFunc(KeyboardInput);
	glutMainLoop();

	return 0;
}

