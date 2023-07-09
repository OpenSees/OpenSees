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
                                                                        
// $Revision: 1.3 $
// $Date: 2008-02-15 23:43:31 $
// $Source: /usr/local/cvs/OpenSees/SRC/renderer/main.cpp,v $
                                                                        
                                                                        
#include <stdlib.h>
#include <iostream>
using std::cin;

#include <OPS_Globals.h>
#include <StandardStream.h>
#include <OPS_Stream.h>

StandardStream sserr;
OPS_Stream *opserrPtr  = &sserr;

double        ops_Dt = 0;
Domain       *ops_TheActiveDomain = 0;
Element      *ops_TheActiveElement = 0;

#include <OpenGLRenderer.h>
#include <PlainMap.h>

#ifdef _AGL

#include <OpenGL/glu.h>
#include <OpenGL/gl.h>

#elif _GLX

#include <GL/glu.h>

#endif

static void init(void)
{
  
   glEnable (GL_BLEND);
   glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
   /*
   glShadeModel (GL_FLAT);
   glClearColor (0.0, 0.0, 0.0, 0.0);
   */
}

static void drawLeftTriangle(void)
{
   /* draw yellow triangle on LHS of screen */

   glBegin (GL_TRIANGLES);
      glColor4f(1.0, 1.0, 0.0, 0.75);
      glVertex3f(0.1, 0.9, 0.0); 
      glVertex3f(0.1, 0.1, 0.0); 
      glVertex3f(0.7, 0.5, 0.0); 
   glEnd();
}

static void drawRightTriangle(void)
{
   /* draw cyan triangle on RHS of screen */

   glBegin (GL_TRIANGLES);
      glColor4f(0.0, 1.0, 1.0, 0.75);
      glVertex3f(0.9, 0.9, 0.0); 
      glVertex3f(0.3, 0.5, 0.0); 
      glVertex3f(0.9, 0.1, 0.0); 
   glEnd();
}


int doImage(Renderer &theRenderer) {

  theRenderer.startImage();

  init();
  drawLeftTriangle();
  drawRightTriangle();

  theRenderer.doneImage();
  return 0;
}

int main(int argc, char **argv)
{
  PlainMap theMap;
  OpenGLRenderer theRenderer1("ONE", 10 , 10 , 100, 200, theMap);
  OpenGLRenderer theRenderer2("TWO", 100 , 100 , 300, 400, theMap);


  theRenderer2.setFillMode("fill");
  doImage(theRenderer1);
  doImage(theRenderer2);

  char d;
  opserr << "Enter a char: " ;
  cin >> d;

  theRenderer2.setFillMode("wire");
  theRenderer1.setFillMode("fill");
  doImage(theRenderer1);
  doImage(theRenderer2);

  opserr << "Enter a char: " ;
  cin >> d;

  exit(0);
  return 0;
}

