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

#define VISUAL
#define Periodic_Show
#include "read.h"
#include "../shared/c2dvector.h"

unsigned int window_width = 680;
unsigned int window_height = 680;
const unsigned int max_height = 680;
const unsigned int max_width = 1300;

cv::VideoWriter writer;
SceneSet* sceneset;
string global_address;
string global_name;

int t = 0;
bool save = false;
bool frame_maker = false;
bool magnify = false;
SavingVector box_dim;
SavingVector r_lense;
SavingVector r_image;
float d0, d1;

void Init()
{
	glClearColor (1.0, 1.0, 1.0, 0.0);
	glEnable (GL_BLEND);
	glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glShadeModel (GL_FLAT);
	Init_Circle();
	float sigma = 1.0;
	VisualParticle::radius = sigma / 2;
	VisualParticle::tail = 1;
	VisualParticle::thickness = 1.5;

//	r_lense.x = 25;
//	r_lense.y = 25;
//	r_image.x = -20;
//	r_image.y = -20;
	r_lense.x = -38;
	r_lense.y = -40;
	r_image.x = 20;
	r_image.y = 20;
	d0 = 5;
	d1 = 25;
}

void Save_Movie()
{
		static cv::Mat img(window_height,window_width,CV_8UC3);
		glPixelStorei(GL_PACK_ALIGNMENT, (img.step & 3) ? 1 : 4);
		glPixelStorei(GL_PACK_ROW_LENGTH, img.step/img.elemSize());
		glReadPixels(0, 0, img.cols, img.rows, GL_BGR, GL_UNSIGNED_BYTE, img.data);
//		cv::flip(img, flipped, 0);
		writer << img;
}

void Save_Image(string address)
{
		static cv::Mat img(window_height,window_width,CV_8UC3);
		glPixelStorei(GL_PACK_ALIGNMENT, (img.step & 3) ? 1 : 4);
		glPixelStorei(GL_PACK_ROW_LENGTH, img.step/img.elemSize());
		glReadPixels(0, 0, img.cols, img.rows, GL_BGR, GL_UNSIGNED_BYTE, img.data);
		cv::flip(img, img, 0);
		cv::imwrite(address, img); //this has a bug in current version
}

void Save_Frame()
{
		static int num = 0;
		stringstream address;
		if (num < 10)
			address << "img000" << num << ".png";
		else
			if (num < 100)
				address << "img00" << num << ".png";
			else
				if (num < 1000)
					address << "img0" << num << ".png";
				else
					address << "img" << num << ".png";
		Save_Image(address.str().c_str());
		num++;
}

void Magnify(SavingVector r0, float d0, SavingVector r1, float d1)
{
	sceneset->scene[t].Magnify(r0,d0,r1,d1);
	cout << (d1 / d0)*VisualParticle::thickness << endl;
}

void Draw_Color_Wheel(SavingVector r, float R0, float R1)
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
	SavingVector r0,r1;

	sceneset->scene[t].Draw();
	if (magnify)
		Magnify(r_lense,d0,r_image,d1);

	SavingVector l(sceneset->L);
//	l.x -= 0.2;
//	l.y -= 0.2;

	float small_radius, big_radius;
	float l_small = min(l.x,l.y);

	small_radius = l_small / 40;
	big_radius = l_small / 8;

	r0.y = l.y - 1.2*big_radius;
	r0.x = l.x - 1.2*big_radius;
//	r0.x = -r0.x;

	Draw_Color_Wheel(r0,small_radius,big_radius);

	glLineWidth(2);
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	glColor4f(0, 0, 0,1);
	glBegin(GL_POLYGON);
	glVertex2f(-l.x, -l.y);
	glVertex2f(l.x, -l.y);
	glVertex2f(l.x, l.y);
	glVertex2f(-l.x, l.y);
	glEnd();
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	glLineWidth(1);

	if (save && (t < (sceneset->Nf-1)))
	{
		Save_Movie();
	}

	if (frame_maker && (t < (sceneset->Nf-1)))
	{
		Save_Frame();
	}

	glutSwapBuffers();
}


void Next_Frame(void)
{
	t++;
	if (t < 0)
		t = 0;
	if (t >= sceneset->Nf)
		t = sceneset->Nf - 1;
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
		cout << "Saving a snapshot" << endl;
//		sceneset->Plot_Fields(16, t, sceneset->info);
		stringstream address("");
		address << "screen-shot-t=" << sceneset->scene[t].t << "-" << global_name << ".png";
		global_address = address.str().c_str();
		Save_Image(global_address);
	}
	if ((key == 72) || (key == 104))
		frame_maker = !frame_maker;
	if ((key == 73) || (key == 105))
	{
		Save_Frame();
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
			r_lense.x = (2*x*sceneset->L.x) / window_width - sceneset->L.x;
			r_lense.y = sceneset->L.y - (2*y*sceneset->L.y) / window_height;
			cout << x << "\t" << y << "\t" << window_width << endl;
		}
		if (button == GLUT_RIGHT_BUTTON)
		{
			r_image.x = (2*x*sceneset->L.x) / window_width - sceneset->L.x;
			r_image.y = sceneset->L.y - (2*y*sceneset->L.y) / window_height;
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
	glOrtho(-box_dim.x, box_dim.x,  -box_dim.y, box_dim.y, -1.0, 1.0);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}


void Welcome()
{
	cout << "Welcome to visualization program. The following keys and their function help you to better use this program:\n" << endl;
	cout << "To stop movie: space or p" << endl;
	cout << "Next frame: right arrow" << endl;
	cout << "Previous frame: left arrow" << endl;
	cout << "Next 10 frame: up arrow" << endl;
	cout << "Previous 10 frame: down arrow" << endl;
	cout << "Next 100 frame: page up" << endl;
	cout << "Previous 100 frame: page down" << endl;
	cout << "To save current snapshot and field values: s or S" << endl;
	cout << "To activate or diactivate frame saving: h or H" << endl;
	cout << "To save current frame: i or I" << endl;
	cout << "To activate or diactivate magnifier: m" << endl;
	cout << "To make magnifier window bigger: a" << endl;
	cout << "To make magnifier window smaller: z" << endl;
	cout << "To make magnified window bigger: A" << endl;
	cout << "To make magnified window bigger: Z" << endl;
}

int main(int argc, char** argv)
{
	SavingVector::Init_Rand(321);

	int glargc = 0;
	char** glargv;

	if (argc > 2)
		save = true;

	Welcome();	
	sceneset = new SceneSet(argv[argc-1]);
	bool read_state = sceneset->Read();
//	sceneset->L = 60;

	if (read_state)
	{
		box_dim = sceneset->L;

		window_height = max_height;
		window_width = (unsigned int) round(box_dim.x*window_height / box_dim.y);
		if (window_width > max_width)
		{
			window_width = max_width;
			window_height = (unsigned int) round(box_dim.y*window_width / box_dim.x);
		}

		global_name = argv[argc-1];
		string::size_type position_of_txt = global_name.find("-r-v", 0);
		if (position_of_txt < 20000)
			global_name.erase(position_of_txt);

		stringstream address("");
		address << global_name << ".mpg";

		if (save)
		{
//			writer.open(address.str().c_str(), CV_FOURCC('F','L','V','1'), 20, cv::Size(window_height,window_width), true);
			writer.open(address.str().c_str(), CV_FOURCC('P','I','M','1'), 20, cv::Size(window_height,window_width), true);
		}
		else
			cv::VideoCapture cap(0);

		glutInit(&glargc, glargv);
		glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB);
		glutInitWindowSize (window_width, window_height);
		glutInitWindowPosition (100, 100);
		glutCreateWindow ("SPP Visualizer");
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

