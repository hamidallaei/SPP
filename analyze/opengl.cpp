#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#include <opencv/cv.h>
#include <opencv/highgui.h>

#include <cstdlib>
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include "read.h"
#include "../shared/c2dvector.h"

unsigned int window_width = 650;
unsigned int window_height = 650;

cv::VideoWriter writer;
SceneSet* sceneset;
string global_address;

int t = 0;
bool save = false;
float box_dim = 0;

void Init()
{
	glClearColor (1.0, 1.0, 1.0, 0.0);
	glShadeModel (GL_FLAT);
	Init_Circle();
	float sigma = 0.5;
	VisualParticle::radius = sigma / 2;
	VisualParticle::tail = 1*VisualParticle::radius;
	VisualParticle::color.red = 0.8;
	VisualParticle::color.green = 0.0;
	VisualParticle::color.blue = 0.0;
}

void Save_Movie()
{
		static cv::Mat img(window_width,window_height,CV_8UC3);
		glPixelStorei(GL_PACK_ALIGNMENT, (img.step & 3) ? 1 : 4);
		glPixelStorei(GL_PACK_ROW_LENGTH, img.step/img.elemSize());
		glReadPixels(0, 0, img.cols, img.rows, GL_BGR, GL_UNSIGNED_BYTE, img.data);
//		cv::flip(img, flipped, 0);
		writer << img;
}

void Save_Image(string address)
{
		static cv::Mat img(window_width,window_height,CV_8UC3);
		glPixelStorei(GL_PACK_ALIGNMENT, (img.step & 3) ? 1 : 4);
		glPixelStorei(GL_PACK_ROW_LENGTH, img.step/img.elemSize());
		glReadPixels(0, 0, img.cols, img.rows, GL_BGR, GL_UNSIGNED_BYTE, img.data);
		cv::imwrite(address, img); //this has a bug in current version
}

void Display()
{
	glClear(GL_COLOR_BUFFER_BIT);
	//glPushMatrix();
	//glPopMatrix();

	sceneset->scene[t].Draw();

	if (save && (t < (sceneset->scene.size()-1)))
	{
		Save_Movie();
	}

	glutSwapBuffers();
}


void Next_Frame(void)
{
	t++;
	if (t < 0)
		t = 0;
	if (t >= sceneset->scene.size())
		t = sceneset->scene.size() - 1;
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
	if ((key == 83) || (key == 115))
	{
		sceneset->Plot_Fields(41, t, sceneset->info);
		Save_Image(global_address);
	}
	if ((key == 32) || (key == 112))
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
	C2DVector::Init_Rand(321);

	int glargc = 0;
	char** glargv;

	if (argc > 2)
		save = true;

	sceneset = new SceneSet(argv[argc-1]);
	bool read_state = sceneset->Read();

	if (read_state)
	{
		box_dim = sceneset->L;

		string name = argv[argc-1];
		string::size_type position_of_txt = name.find("-r-v", 0);
		name.erase(position_of_txt);

		stringstream address("");
		address << name << ".mpg";

		if (save)
		{
//			writer.open(address.str().c_str(), CV_FOURCC('F','L','V','1'), 20, cv::Size(window_width,window_height), true);
			writer.open(address.str().c_str(), CV_FOURCC('P','I','M','1'), 20, cv::Size(window_width,window_height), true);
		}
		else
			cv::VideoCapture cap(0);

		address.str("");
		address << "screen-shot-" << name << ".png";
		global_address = address.str().c_str();

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
	}
	else
		cout << "Can not open the file" << endl;

	return 0;
}

