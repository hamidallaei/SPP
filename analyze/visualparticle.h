#ifndef _VISUALPARTICLE_
#define _VISUALPARTICLE_

#ifdef VISUAL
	const int circle_points_num = 100;
	const double alpha_value = 0.3;
	const double center_line_illum = 0.4;

	GLfloat unitcircle[2*circle_points_num];
	GLuint indices[circle_points_num];


	void Init_Circle()
	{
		for (int i = 0; i < circle_points_num; i++)
		{
			unitcircle[2*i] = cos(2*i*PI/circle_points_num);
			unitcircle[2*i+1] = sin(2*i*PI/circle_points_num);
			indices[i] = i;
		}
	}

	struct RGB{
		float red,green,blue;
	};

	void Draw_Circle(SavingVector r, float theta, float radius, float thickness, RGB color, GLenum mode)
	{
		glEnableClientState (GL_VERTEX_ARRAY);

		glColor4f(color.red, color.green, color.blue, alpha_value);
		glLineWidth(thickness);

		glPolygonMode(GL_FRONT_AND_BACK, mode);
		glPolygonMode(GL_FRONT_AND_BACK, mode);

		glVertexPointer (2, GL_FLOAT, 0, unitcircle);

		glTranslatef(r.x,r.y,0);
		glScalef(radius,radius,1);

		glDrawElements(GL_POLYGON,circle_points_num,GL_UNSIGNED_INT,indices);

		glBegin(GL_LINES);
			glVertex2f(0, 0);
			glVertex2f(-cos(theta),-sin(theta));
		glEnd();
			glColor4f(color.red, color.green, color.blue, 1.0);
			glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
			glDrawElements(GL_POLYGON,circle_points_num,GL_UNSIGNED_INT,indices);
	
		glDisableClientState( GL_VERTEX_ARRAY );
		glLoadIdentity();
		glLineWidth(1);
	}

	void Draw_Circle(SavingVector r, float radius, float thickness, RGB color, GLenum mode)
	{
		glEnableClientState (GL_VERTEX_ARRAY);

		glColor4f(color.red, color.green, color.blue, alpha_value);
		glLineWidth(thickness);

		glPolygonMode(GL_FRONT_AND_BACK, mode);
		glPolygonMode(GL_FRONT_AND_BACK, mode);

		glVertexPointer (2, GL_FLOAT, 0, unitcircle);

		glTranslatef(r.x,r.y,0);
		glScalef(radius,radius,1);

		glDrawElements(GL_POLYGON,circle_points_num,GL_UNSIGNED_INT,indices);

		glColor4f(color.red, color.green, color.blue, 1.0);
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		glDrawElements(GL_POLYGON,circle_points_num,GL_UNSIGNED_INT,indices);
	
		glDisableClientState( GL_VERTEX_ARRAY );
		glLoadIdentity();
		glLineWidth(1);
	}

	void HSV_To_RGB(float h, float s, float v,RGB& color)
	{
		int index;
		float f, p, q, t;
		if( s == 0 ) {
			// achromatic (grey)
			color.red = color.green = color.blue = v;
			return;
		}

		h /= 2*M_PI;			// sector 0 to 5
		h *= 6;
		index = floor( h );
		f = h - index;			// factorial part of h
		p = v * ( 1 - s );
		q = v * ( 1 - s * f );
		t = v * ( 1 - s * ( 1 - f ) );

		switch( index ) {
			case 0:
				color.red = v;
				color.green = t;
				color.blue = p;
				break;
			case 1:
				color.red = q;
				color.green = v;
				color.blue = p;
				break;
			case 2:
				color.red = p;
				color.green = v;
				color.blue = t;
				break;
			case 3:
				color.red = p;
				color.green = q;
				color.blue = v;
				break;
			case 4:
				color.red = t;
				color.green = p;
				color.blue = v;
				break;
			default:		// case 5:
				color.red = v;
				color.green = p;
				color.blue = q;
				break;
		}
	}
#endif

class BasicParticle00{
public:
	SavingVector r;
	void Write(std::ostream&);
};

void BasicParticle00::Write(std::ostream& os)
{
	r.write(os);
}

class BasicParticle0{
public:
	SavingVector r;
	Saving_Real theta;
	void Write(std::ostream&);
};

void BasicParticle0::Write(std::ostream& os)
{
	r.write(os);
	Saving_Real temp_float = (Saving_Real) theta;
	os.write((char*) &temp_float,sizeof(Saving_Real) / sizeof(char));
}

class VisualMembrane : public BasicParticle00 {
public:
	#ifdef VISUAL
		static float radius;
		static float thickness;
		static GLenum mode; // GL_FILL (filled) or GL_LINE (empty) or GL_POINT

//		virtual void Draw(float, float, float);
		void Draw(float, float, float, bool);
		void Draw_Translated_Scaled(float, float, float, bool, SavingVector);
		void Draw_Magnified(SavingVector r0, float d0, SavingVector r1, float d1);
		void Draw_Magnified(SavingVector r0, float d0, SavingVector r1, float d1, bool flag);
	#endif
};

#ifdef VISUAL
	float VisualMembrane::radius = 0.5;
	float VisualMembrane::thickness = 1;
	GLenum VisualMembrane::mode = GL_FILL;

	void VisualMembrane::Draw(float x0 = 0, float y0 = 0, float scale = 1, bool flag = false)
	{
		RGB color;
		color.red = 0.5;
		color.green = 0.5;
		color.blue = 0.5;
		if (flag)
		{
			color.red = 0.0;
			color.green = 0.0;
			color.blue = 0.0;
		}

		SavingVector r_temp(r);
		r_temp.x -= x0;
		r_temp.y -= y0;

		Draw_Circle(r_temp, radius*scale, thickness, color, mode);
	}

/*	void VisualMembrane::Draw_Translated_Scaled(float x0 = 0, float y0 = 0, float scale = 1, bool flag = false, SavingVector rcm)*/
/*	{*/
/*		RGB color;*/
/*		color.red = 0.5;*/
/*		color.green = 0.5;*/
/*		color.blue = 0.5;*/
/*		if (flag)*/
/*		{*/
/*			color.red = 0.0;*/
/*			color.green = 0.0;*/
/*			color.blue = 0.0;*/
/*		}*/

/*		Draw_Circle(r - r_cm, radius*scale, thickness, color, mode);*/
/*	}*/

	void VisualMembrane::Draw_Magnified(SavingVector r0, float d0, SavingVector r1, float d1)
	{
		float scale = d1 / d0;
		SavingVector p;
		p.x = r.x - r0.x;
		p.y = r.y - r0.y;
		if (fabs(p.x) < (d0-radius) && fabs(p.y) < (d0-radius))
		{
			p *= scale;
			p += r1;
			const int thickness_factor = 3; // 2000

			SavingVector null(0);
			null.Null();

			RGB color;
			color.red = 0.5;
			color.green = 0.5;
			color.blue = 0.5;

			Draw_Circle(p, radius*scale, scale*thickness / thickness_factor, color, mode);
		}
	}

	void VisualMembrane::Draw_Magnified(SavingVector r0, float d0, SavingVector r1, float d1, bool flag)
	{
		float scale = d1 / d0;
		SavingVector p;
		p.x = r.x - r0.x;
		p.y = r.y - r0.y;
		if (fabs(p.x) < (d0-radius) && fabs(p.y) < (d0-radius))
		{
			p *= scale;
			p += r1;
			const int thickness_factor = 3; // 2000

			SavingVector null(0);
			null.Null();

			RGB color;
			color.red = 0.5;
			color.green = 0.5;
			color.blue = 0.5;
			if (flag)
			{
				color.red = 0.0;
				color.green = 0.0;
				color.blue = 0.0;
			}


			Draw_Circle(p, radius*scale, scale*thickness / thickness_factor, color, mode);
		}
	}
#endif


class VisualParticle : public BasicParticle0 {
public:
	#ifdef VISUAL
		static float radius;
		static float tail;
		static float thickness;
		static GLenum mode; // GL_FILL (filled) or GL_LINE (empty) or GL_POINT

		void Find_Color(RGB&);
		void Draw(float, float, float);
		void Draw_Translated_Scaled(float, float, float, SavingVector r_cm);
		void Draw_Magnified(SavingVector r0, float d0, SavingVector r1, float d1);
	#endif
};

#ifdef VISUAL
	float VisualParticle::radius;
	float VisualParticle::tail;
	float VisualParticle::thickness;
	GLenum VisualParticle::mode = GL_FILL;

	void VisualParticle::Find_Color(RGB& input_color)
	{
		float f, p, q, t, s = 1, value = 1;

		theta = theta - floor(theta/(2*M_PI))*2*M_PI;
		HSV_To_RGB(theta,1,1, input_color);
	}

	void VisualParticle::Draw(float x0 = 0, float y0 = 0, float scale = 1)
	{
		RGB color;
		Find_Color(color);

		SavingVector r_temp(r);
		r_temp.x -= x0;
		r_temp.y -= y0;

		Draw_Circle(r_temp, theta, radius*scale, thickness, color, mode);
	}

/*	void VisualParticle::Draw_Translated_Scaled(float x0 = 0, float y0 = 0, float scale = 1, SavingVector r_cm)*/
/*	{*/
/*		RGB color;*/
/*		Find_Color(color);*/
/*		Draw_Circle(r - r_cm, theta, radius*scale, thickness, color, mode);*/
/*	}*/

	void VisualParticle::Draw_Magnified(SavingVector r0, float d0, SavingVector r1, float d1)
	{
		float scale = d1 / d0;
		SavingVector p;
		p.x = r.x - r0.x;
		p.y = r.y - r0.y;
		if (fabs(p.x) < (d0-radius) && fabs(p.y) < (d0-radius))
		{
			p *= scale;
			p += r1;
			const int thickness_factor = 3; // 2000

			RGB color;
			Find_Color(color);
			Draw_Circle(p, theta, radius*scale, scale*thickness / thickness_factor, color, mode);
		}
	}
#endif


class VisualChain : public VisualParticle {
public:
	static int chain_length;

	#ifdef VISUAL
		void Draw(float, float, float);
		void Draw(float, float, float, SavingVector);
		void Draw_Magnified(SavingVector r0, float d0, SavingVector r1, float d1);
	#endif
};

int VisualChain::chain_length = 2;

#ifdef VISUAL
	void VisualChain::Draw(float x0 = 0, float y0 = 0, float scale = 1)
	{
		SavingVector v_temp, r_temp;
		v_temp.x = cos(theta);
		v_temp.y = sin(theta);

		RGB color;
		Find_Color(color);

		for (int i = 0; i < chain_length; i++)
		{
			SavingVector s_i = v_temp*(((1-chain_length)/2.0 + i)/(chain_length-1));

			r_temp = r + s_i;
			r_temp.x -= x0;
			r_temp.y -= y0;

			Draw_Circle(r_temp, theta, radius*scale, thickness, color, mode);
		}

		glLineWidth(2);
		glColor4f(center_line_illum*color.red, center_line_illum*color.green, center_line_illum*color.blue, 1.0);
		glTranslatef(r.x,r.y,0);
		glBegin(GL_LINES);
			glVertex2f(-0.7*cos(theta),-0.7*sin(theta));
			glVertex2f(0.7*cos(theta),0.7*sin(theta));
		glEnd();
	
		glDisableClientState( GL_VERTEX_ARRAY );
		glLoadIdentity();
		glLineWidth(1);
	}

/*	void VisualChain::Draw_Translated_Scaled(float x0 = 0, float y0 = 0, float scale = 1, SavingVector r_cm)*/
/*	{*/
/*		SavingVector v_temp, r_temp;*/
/*		v_temp.x = cos(theta);*/
/*		v_temp.y = sin(theta);*/

/*		RGB color;*/
/*		Find_Color(color);*/

/*		for (int i = 0; i < chain_length; i++)*/
/*		{*/
/*			SavingVector s_i = v_temp*(((1-chain_length)/2.0 + i)/(chain_length-1));*/

/*			r_temp = r + s_i - r_cm;*/

/*			Draw_Circle(r_temp, theta, radius*scale, thickness, color, mode);*/
/*		}*/

/*		glLineWidth(2);*/
/*		glColor4f(center_line_illum*color.red, center_line_illum*color.green, center_line_illum*color.blue, 1.0);*/
/*		glTranslatef(r.x - r_cm.x,r.y - r_cm.y,0);*/
/*		glBegin(GL_LINES);*/
/*			glVertex2f(-0.7*cos(theta),-0.7*sin(theta));*/
/*			glVertex2f(0.7*cos(theta),0.7*sin(theta));*/
/*		glEnd();*/
/*	*/
/*		glDisableClientState( GL_VERTEX_ARRAY );*/
/*		glLoadIdentity();*/
/*		glLineWidth(1);*/
/*	}*/

	void VisualChain::Draw_Magnified(SavingVector r0, float d0, SavingVector r1, float d1)
	{
		float scale = d1 / d0;
		SavingVector v_temp, p;
		SavingVector r_new = (r - r0)*scale + r1;
		v_temp.x = cos(theta);
		v_temp.y = sin(theta);

		RGB color;
		Find_Color(color);

		for (int i = 0; i < chain_length; i++)
		{
			SavingVector s_i = v_temp*(((1-chain_length)/2.0 + i)/(chain_length-1));
			SavingVector p;
			p = r - r0 + s_i;

			if (fabs(p.x) < (d0-radius) && fabs(p.y) < (d0-radius))
			{
				p *= scale;
				p += r1;
				const int thickness_factor = 3; // 2000

				Draw_Circle(p, theta, radius*scale, scale*thickness / thickness_factor, color, mode);
			}
		}

		if (fabs(r.x - r0.x) < (d0-radius) && fabs(r.y - r0.y) < (d0-radius))
		{
			float lw = round(scale*2);
			glLineWidth(lw);
			glColor4f(center_line_illum*color.red, center_line_illum*color.green, center_line_illum*color.blue, 1.0);
			glTranslatef(r_new.x,r_new.y,0);
			glBegin(GL_LINES);
				glVertex2f(-0.7*scale*cos(theta),-0.7*scale*sin(theta));
				glVertex2f(0.7*scale*cos(theta),0.7*scale*sin(theta));
			glEnd();

			glDisableClientState( GL_VERTEX_ARRAY );
			glLoadIdentity();
			glLineWidth(1);
		}

	}
#endif

#endif
