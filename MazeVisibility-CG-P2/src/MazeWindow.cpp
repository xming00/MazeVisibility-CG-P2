/************************************************************************
     File:        MazeWindow.cpp

     Author:     
                  Stephen Chenney, schenney@cs.wisc.edu
     Modifier
                  Yu-Chi Lai, yu-chi@cs.wisc.edu

     Comment:    
						(c) 2001-2002 Stephen Chenney, University of Wisconsin at Madison

						Class header file for the MazeWindow class. The MazeWindow is
						the window in which the viewer's view of the maze is displayed.
		

     Platform:    Visio Studio.Net 2003 (converted to 2005)

*************************************************************************/

#include "MazeWindow.h"
#include <Fl/math.h>
#include <Fl/gl.h>
#include <GL/glu.h>
#include <stdio.h>


//*************************************************************************
//
// * Constructor
//=========================================================================
MazeWindow::
MazeWindow(int x, int y, int width, int height, const char *label,Maze *m)
	: Fl_Gl_Window(x, y, width, height, label)
//=========================================================================
{
	maze = m;

	// The mouse button isn't down and there is no key pressed.
	down = false;
	z_key = 0;
}


//*************************************************************************
//
// * Set the maze. Also causes a redraw.
//=========================================================================
void MazeWindow::
Set_Maze(Maze *m)
//=========================================================================
{
	// Change the maze
	maze = m;

	// Force a redraw
	redraw();
}

//*************************************************************************
//
//  
//=========================================================================
void MazeWindow::NormalizeVector(float* v3)
{
//=========================================================================
	float v3_Len = sqrt(pow((v3[0]), 2) + pow((v3[1]), 2)+ pow((v3[2]), 2));
	v3[0] /= v3_Len;
	v3[1] /= v3_Len;
	v3[2] /= v3_Len;
}

//*************************************************************************
//
//  ¥~¿n
//=========================================================================
void MazeWindow::ComputeNormalOfPlane(float* result, float* a, float* b)
{
//=========================================================================
	result[0] = a[1] * b[2] - a[2] * b[1];
	result[1] = a[2] * b[0] - a[0] * b[2];
	result[2] = a[0] * b[1] - a[1] * b[0];
}

//*************************************************************************
//
//  perspective.
//=========================================================================
void MazeWindow::perspectivef2(float fovyInDegrees, float aspectRatio, float znear, float zfar)
{
	float ymax, xmax;
    float temp, temp2, temp3, temp4;
	float matrix[16];
	ymax = znear * tanf(fovyInDegrees * M_PI / 360.0);
    // ymin = -ymax;
    // xmin = -ymax * aspectRatio;
    xmax = ymax * aspectRatio;
	glMatrixMode(GL_PROJECTION);
	glGetFloatv(GL_PROJECTION_MATRIX, matrix);
	frustumf2(matrix, -xmax, xmax, -ymax, ymax, znear, zfar);

}

void MazeWindow::frustumf2(float* matrix, float left, float right, float bottom, float top,	float znear, float zfar) 
{
	float temp, temp2, temp3, temp4;
	temp = 2.0 * znear;
	temp2 = right - left;
	temp3 = top - bottom;
	temp4 = zfar - znear;
	matrix[0] = temp / temp2;
	matrix[1] = 0.0;
	matrix[2] = 0.0;
	matrix[3] = 0.0;
	matrix[4] = 0.0;
	matrix[5] = temp / temp3;
	matrix[6] = 0.0;
	matrix[7] = 0.0;
	matrix[8] = (right + left) / temp2;
	matrix[9] = (top + bottom) / temp3;
	matrix[10] = (-zfar - znear) / temp4;
	matrix[11] = -1.0;
	matrix[12] = 0.0;
	matrix[13] = 0.0;
	matrix[14] = (-temp * zfar) / temp4;
	matrix[15] = 0.0;
	glLoadMatrixf(matrix);
}


//*************************************************************************
//
//  lookat.
//=========================================================================
void MazeWindow::lookat( float* eyePosition3D, float* center3D, float* upVector3D)
{
	float forward[3], side[3], up[3];
	float matrix2[16], resultMatrix[16];
	// --------------------
	forward[0] = center3D[0] - eyePosition3D[0];
	forward[1] = center3D[1] - eyePosition3D[1];
	forward[2] = center3D[2] - eyePosition3D[2];
	NormalizeVector(forward);
	// --------------------
	// Side = forward x up
	ComputeNormalOfPlane(side, forward, upVector3D);
	NormalizeVector(side);
	//--------------------
	// Recompute up as: up = side x forward
	ComputeNormalOfPlane(up, side, forward);
	// --------------------
	matrix2[0] = side[0];
	matrix2[4] = side[1];
	matrix2[8] = side[2];
	matrix2[12] = 0.0;
	// --------------------
	matrix2[1] = up[0];
	matrix2[5] = up[1];
	matrix2[9] = up[2];
	matrix2[13] = 0.0;
	// --------------------
	matrix2[2] = -forward[0];
	matrix2[6] = -forward[1];
	matrix2[10] = -forward[2];
	matrix2[14] = 0.0;
	// --------------------
	matrix2[3] = matrix2[7] = matrix2[11] = 0.0;
	matrix2[15] = 1.0;
	// --------------------
	glMatrixMode(GL_MODELVIEW);
	
	glMultMatrixf(matrix2);
	glTranslatef(-eyePosition3D[0], -eyePosition3D[1], -eyePosition3D[2]);
	// --------------------

}
//*************************************************************************
//
//  lookat.
//=========================================================================
void MazeWindow::seeMat(float* mat)
{
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			std::cout << mat[j * 4 + i] << " ";
		}
		std::cout << std::endl;
	}

}

//*************************************************************************
//
// * draw() method invoked whenever the view changes or the window
//   otherwise needs to be redrawn.
//=========================================================================
void MazeWindow::
draw(void)
//=========================================================================
{
	float   focal_length;
	float mat[16] = {0};
	for (int i = 0; i < 16; i++)
	{
		mat[i] = i;
	}
	if ( ! valid() ) {
		// The OpenGL context may have been changed
		// Set up the viewport to fill the window.
		glViewport(0, 0, w(), h());

		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		// We are using orthogonal viewing for 2D. This puts 0,0 in the
		// middle of the screen, and makes the image size in view space
		// the same size as the window.
		gluOrtho2D(-w() * 0.5, w() * 0.5, -h() * 0.5, h() * 0.5);

		// Sets the clear color to black.
		glClearColor(0.0, 0.0, 0.0, 1.0);
	}

	// Clear the screen.
	glClear(GL_COLOR_BUFFER_BIT);

	
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();


	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	glBegin(GL_QUADS);
		// Draw the "floor". It is an infinite plane perpendicular to
		// vertical, so we know it projects to cover the entire bottom
		// half of the screen. Walls of the maze will be drawn over the top
		// of it.
		glColor3f(0.2f, 0.2f, 0.2f);
		glVertex2f(-w() * 0.5f, -h() * 0.5f);
		glVertex2f( w() * 0.5f, -h() * 0.5f);
		glVertex2f( w() * 0.5f, 0.0       );
		glVertex2f(-w() * 0.5f, 0.0       );

		// Draw the ceiling. It will project to the entire top half
		// of the window.
		glColor3f(0.4f, 0.4f, 0.4f);
		glVertex2f( w() * 0.5f,  h() * 0.5f);
		glVertex2f(-w() * 0.5f,  h() * 0.5f);
		glVertex2f(-w() * 0.5f, 0.0       );
		glVertex2f( w() * 0.5f, 0.0       );
	glEnd();


	if ( maze ) {
		// Set the focal length. We can do this because we know the
		// field of view and the size of the image in view space. Note
		// the static member function of the Maze class for converting
		// radians to degrees. There is also one defined for going backwards.
		focal_length = w()
						 / (float)(2.0*tan(Maze::To_Radians(maze->viewer_fov)*0.5));
		
		// Draw the 3D view of the maze (the visible walls.) You write this.
		// Note that all the information that is required to do the
		// transformations and projection is contained in the Maze class,
		// plus the focal length.
		glClear(GL_DEPTH_BUFFER_BIT);

		float aspect = (float)w() / h();
		//gluPerspective(maze->viewer_fov, aspect, 0.01, (double)200);
		perspectivef2(maze->viewer_fov, aspect, 0.01, (double)200);
		float viewer_pos[3] = { maze->viewer_posn[Maze::Y], 0.0f ,  maze->viewer_posn[Maze::X] }; // xy -> zx
		float eye[3] = { viewer_pos[Maze::X], viewer_pos[Maze::Y], viewer_pos[Maze::Z] };

		float center[3] = { viewer_pos[Maze::X] + sin(Maze::To_Radians(maze->viewer_dir)),
							viewer_pos[Maze::Y],
							viewer_pos[Maze::Z] + cos(Maze::To_Radians(maze->viewer_dir)) };
		float up[3] = { 0.0f, 1.0f, 0.0f };
	/*	gluLookAt(maze->viewer_posn[Maze::Y], 0.0f, maze->viewer_posn[Maze::X],
			viewer_pos[Maze::X] + sin(Maze::To_Radians(maze->viewer_dir)),
			viewer_pos[Maze::Y],
			viewer_pos[Maze::Z] + cos(Maze::To_Radians(maze->viewer_dir)),
			0.0f, 1.0f, 0.0f);*/
		lookat(eye, center, up);
		
		glGetFloatv(GL_PROJECTION_MATRIX, mat);
		seeMat(mat);
		maze->Draw_View(focal_length);
	}
}


//*************************************************************************
//
// *
//=========================================================================
bool MazeWindow::
Drag(float dt)
//=========================================================================
{
	float   x_move, y_move, z_move;

	if ( down ) {
		int	dx = x_down - x_last;
		int   dy = y_down - y_last;
		float dist;

		// Set the viewing direction based on horizontal mouse motion.
		maze->Set_View_Dir(d_down + 360.0f * dx / (float)w());

		// Set the viewer's linear motion based on a speed (derived from
		// vertical mouse motion), the elapsed time and the viewing direction.
		dist = 10.0f * dt * dy / (float)h();
		x_move = dist * (float)cos(Maze::To_Radians(maze->viewer_dir));
		y_move = dist * (float)sin(Maze::To_Radians(maze->viewer_dir));
	}
	else {
		x_move = 0.0;
		y_move = 0.0;
	}

	// Update the z posn
	z_move = z_key * 0.1f;
	z_key = 0;

	// Tell the maze how much the view has moved. It may restrict the motion
	// if it tries to go through walls.
	maze->Move_View_Posn(x_move, y_move, z_move);

	return true;
}


//*************************************************************************
//
// *
//=========================================================================
bool MazeWindow::
Update(float dt)
//=========================================================================
{
	// Update the view

	if ( down || z_key ) // Only do anything if the mouse button is down.
		return Drag(dt);

	// Nothing changed, so no need for a redraw.
	return false;
}


//*************************************************************************
//
// *
//=========================================================================
int MazeWindow::
handle(int event)
//=========================================================================
{
	if (!maze)
		return Fl_Gl_Window::handle(event);

	// Event handling routine.
	switch ( event ) {
		case FL_PUSH:
			down = true;
			x_last = x_down = Fl::event_x();
			y_last = y_down = Fl::event_y();
			d_down = maze->viewer_dir;
			return 1;
		case FL_DRAG:
			x_last = Fl::event_x();
			y_last = Fl::event_y();
			return 1;
			case FL_RELEASE:
			down = false;
			return 1;
		case FL_KEYBOARD:
			/*
			if ( Fl::event_key() == FL_Up )	{
				z_key = 1;
				return 1;
			}
			if ( Fl::event_key() == FL_Down ){
				z_key = -1;
				return 1;
			}
			*/
			return Fl_Gl_Window::handle(event);
		case FL_FOCUS:
		case FL_UNFOCUS:
			return 1;
	}

	// Pass any other event types on the superclass.
	return Fl_Gl_Window::handle(event);
}


