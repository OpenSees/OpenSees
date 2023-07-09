/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */
                                                                        
// $Revision: 1.1.1.1 $
// $Date: 2000-09-15 08:23:27 $
// $Source: /usr/local/cvs/OpenSees/SRC/renderer/simple.c,v $
                                                                        
                                                                        
/*
 *  simple.c
 *
 *  WGL/OpenGL program
 *
 *  Copyright (C) 1997 by Nate Robins (ndr@pobox.com)
 *
 *  This program is freely distributable without licensing fees and is
 *  provided without guarantee or warrantee expressed or implied. This
 *  program is not in the public domain.
 */


/* includes */
#include <windows.h>		
#include <gl/gl.h>			
#include <stdio.h>



/* main()
 *  Entry point
 */
int drawIT(void) 
{
	int i;
    HDC       hDC1, hDC2;			/* device context */
    HGLRC     hRC1, hRC2;			/* opengl context */
    HWND      hWnd1, hWnd2;			/* window */

    /* create the windows */
    hWnd1 = oglCreateWindow("ONE", 0, 0, 200, 200, &hRC1, &hDC1);
    if (hWnd1 == NULL)
      exit(1);
	    hWnd2 = oglCreateWindow("TWO", 200, 200, 200, 200, &hRC2, &hDC2);
    if (hWnd2 == NULL)
      exit(1);

    /* now we can start changing state & rendering */
    for (i=0; i<100; i++) {
        /* rotate a triangle around */
		wglMakeCurrent(hDC1, hRC1);
        glClear(GL_COLOR_BUFFER_BIT);
		glRotatef(1.0, 0.0, 0.0, 1.0);
		glBegin(GL_TRIANGLES);
		glColor3f(1.0, 0.0, 0.0);
		glVertex2i( 0,  1);
		glColor3f(0.0, 1.0, 0.0);
		glVertex2i(-1, -1);
		glColor3f(0.0, 0.0, 1.0);
		glVertex2i( 1, -1);
		glEnd();
		glFlush();
		SwapBuffers(hDC1);		/* nop if singlebuffered */

		wglMakeCurrent(hDC2, hRC2);
        glClear(GL_COLOR_BUFFER_BIT);
		glRotatef(1.0, 0.0, 0.0, 1.0);
		glBegin(GL_TRIANGLES);
		glColor3f(0.0, 0.0, 0.0);
		glVertex2i( 0,  1);
		glColor3f(1.0, 1.0, 0.0);
		glVertex2i(-1, -1);
		glColor3f(0.0, 0.0, 1.0);
		glVertex2i( 1, -1);
		glEnd();
		glFlush();
		SwapBuffers(hDC2);		/* nop if singlebuffered */
    }
	getch();
    /* clean up */
	
	oglDeleteWindow("ONE",hWnd1, hRC1, hDC1);
	oglDeleteWindow("TWO",hWnd2, hRC2, hDC2);

	return 0;
	
}

int main(int argv, char **argc) {
	int i;
	for (i=0; i<10; i++)
			drawIT();
	return 0;
}

