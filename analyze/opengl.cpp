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
#include "read-visual.h"
#include "../shared/c2dvector.h"

unsigned int window_width = 660;
unsigned int window_height = 660;

cv::VideoWriter writer;
SceneSet* sceneset;
string global_address;

int t = 0;
bool save = false;
bool magnify = false;
float box_dim = 0;
C2DVector r_lense;
C2DVector r_image;
float d0, d1;

void Init()
{
	glClearColor (1.0, 1.0, 1.0, 0.0);
	glEnable (GL_BLEND);
	glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glShadeModel (GL_FLAT);
	Init_Circle();
	float sigma = 1;
	VisualParticle::radius = sigma / 2;
	VisualParticle::tail = 1;
	VisualParticle::thickness = 1.5;
	VisualParticle::color.red = 0.8;
	VisualParticle::color.green = 0.0;
	VisualParticle::color.blue = 0.0;

//	r_lense.x = 25;
//	r_lense.y = 25;
//	r_image.x = -20;
//	r_image.y = -20;
	r_lense.x = -38;
	r_lense.y = -40;
	r_image.x = 20;
	r_image.y = 20;
	d0 = 5;
	d1 = 5;
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

void Magnify(C2DVector r0, float d0, C2DVector r1, float d1)
{
	sceneset->scene[t].Magnify(r0,d0,r1,d1);
	cout << (d1 / d0)*VisualParticle::thickness << endl;
}

void Draw_Color_Wheel(C2DVector r, float R0, float R1)
{
	float width = R1 - R0;

	int circle_points_num = 500;

		RGB color;
		glBegin( GL_TRIANGLE_FAN );
		glVertex2f(r.x, r.y);
		for (int i = 0; i < circle_points_num; i++)
		{
			float theta = (2*i*M_PI)/circle_points_num;
			HSV_To_RGB(theta, 1, 1, color);
			glColor4f(color.red, color.green, color.blue,0.5);
			glVertex2f(r.x+R1*cos(theta),r.y+R1*sin(theta));
		}
		glEnd();

		circle_points_num = 200;

		glBegin( GL_TRIANGLE_FAN );
		glVertex2f(r.x, r.y);
		glColor4f(1, 1, 1,1);
		for (int i = 0; i < circle_points_num; i++)
		{
			float theta = (2*i*M_PI)/circle_points_num;
			glVertex2f(r.x+R0*cos(theta),r.y+R0*sin(theta));
		}
		glEnd();

//	glScalef(R0,R0,1);

	glDrawElements(GL_POLYGON,circle_points_num,GL_UNSIGNED_INT,indices);

	glDisableClientState( GL_VERTEX_ARRAY );
	glLoadIdentity();
}

void Display()
{
	glClear(GL_COLOR_BUFFER_BIT);
	//glPushMatrix();
	//glPopMatrix();
	C2DVector r0,r1;


	sceneset->scene[t].Draw();
	if (magnify)
		Magnify(r_lense,d0,r_image,d1);
	r0.x = -49;
	r0.y = 49;
	Draw_Color_Wheel(r0,3,10);

	glLineWidth(2);
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	glColor4f(0, 0, 0,1);
	glBegin(GL_POLYGON);
	float l = sceneset->L - 0.2;
	glVertex2f(-l, -l);
	glVertex2f(l, -l);
	glVertex2f(l, l);
	glVertex2f(-l, l);
	glEnd();
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	glLineWidth(1);

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
			t -= 1;
		break;
	}
	Next_Frame();
}

void KeyboardInput(unsigned char key, int x, int y)
{
	static bool stop = true;
	if (key == 109)
	{
		magnify = !magnify;
		Display();
	}
	if ((key == 83) || (key == 115))
	{
		sceneset->Plot_Fields(41, t, sceneset->info);
		Save_Image(global_address);
	}
	float d_jump = 0.5;
	if (key == 97)
	{
		d0 += d_jump;
		cout << d0 << "\t" << d1 << endl;
		Display();
	}
	if (key == 122)
	{
		d0 -= d_jump;
		cout << d0 << "\t" << d1 << endl;
		Display();
	}
	if (key == 65)
	{
		d1 += d_jump;
		cout << d0 << "\t" << d1 << endl;
		Display();
	}
	if (key == 90)
	{
		d1 -= d_jump;
		cout << d0 << "\t" << d1 << endl;
		Display();
	}
	if ((key == 32) || (key == 112))
		stop = !stop;
	if (stop)
		glutIdleFunc(NULL);
	else
		glutIdleFunc(Next_Frame);
}

void MouseInput(int button, int state, int x, int y)
{
	if (state == GLUT_DOWN)
	{
		if (button == GLUT_LEFT_BUTTON)
		{
			r_lense.x = (2*x*sceneset->L) / window_width - sceneset->L;
			r_lense.y = sceneset->L - (2*y*sceneset->L) / window_width;
			cout << x << "\t" << y << "\t" << window_width << endl;
		}
		if (button == GLUT_RIGHT_BUTTON)
		{
			r_image.x = (2*x*sceneset->L) / window_width - sceneset->L;
			r_image.y = sceneset->L - (2*y*sceneset->L) / window_width;
		}
		Display();
	}
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
		glutMouseFunc(MouseInput);
		glutMainLoop();
	}
	else
		cout << "Can not open the file" << endl;

	return 0;
}

