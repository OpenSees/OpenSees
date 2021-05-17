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
                                                                        
// $Revision: 1.20 $
// $Date: 2008-04-15 18:01:55 $
// $Source: /usr/local/cvs/OpenSees/SRC/renderer/OpenGlRenderer.cpp,v $
                                                                        
                                                                        
// Written: fmk 
// Revision: A
//
// Description: This file contains the class definition for OpenGLRenderer.
// OpenGLRenderer is an class which diplays using X11 or openGL.
//
// What: "@(#) OpenGLRenderer.h, revA"

#include <OpenGLRenderer.h>
#include <ColorMap.h>
#include <stdio.h>
#include <stdlib.h>
#include <Matrix.h>
#include <iomanip>
using std::ios;

#include <string.h>

#ifdef _WGL
#include <windows.h>
//#include <GL/glaux.h>
#include <GL/glu.h>
#include <GL/gl.h>

#elif _GLX

#include <GL/gl.h>
#include <GL/glx.h>
#include <GL/glu.h>

#elif _AGL

#include <OpenGL/glu.h>
#include <OpenGL/gl.h>

#endif

#define PARALLEL_MODE 0
#define PERSPECTIVE_MODE 1

#define WIRE_MODE 0
#define FILL_MODE 1

//#include <db.H>
#include <Vector.h>

OpenGLRenderer::OpenGLRenderer(const char *_title, int _xLoc, int _yLoc, 
			       int _width, int _height, 
			       ColorMap &_theMap)
  :Renderer(_title, _theMap),  
  windowTitle(0), height(_height), width(_width), xLoc(_xLoc), yLoc(_yLoc),
  count(-1), theOutputFileName(0), 
  theDevice(0),
  vrp(3), vuv(3), vpn(3), cop(3), ViewMat(4,4), 
  projectionMode(0), vpWindow(4), ProjMat(4,4),
  portWindow(4)
{

  // set the WindowDevices title, height, wdth, xLoc and yLoc
  windowTitle = new char [strlen(_title)+1];
  strcpy(windowTitle, _title);

  fillMode = WIRE_MODE;
  projectionMode = PARALLEL_MODE;
  portWindow(0) = -1.0; portWindow(1) = 1.0;
  portWindow(2) = -1.0;  portWindow(3) = 1.0;

  theDevice = new OpenGlDevice();
  theDevice->WINOPEN(_title, _xLoc, _yLoc, _width, _height);

  theDevice->CLEAR();

  // glEnable (GL_BLEND);
  // glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  glClearColor(1.0f,1.0f,1.0f,1.0f);

  glEnable(GL_DEPTH_TEST);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  glEnd();
  glFlush();
  theDevice->ENDIMAGE();
}




OpenGLRenderer::OpenGLRenderer(const char *_title, int _xLoc, int _yLoc, 
			       int _width, int _height,	
			       ColorMap &_theMap, 
			       const char *outputFileName, 
			       const char *bitmapFileName)
  :Renderer(_title, _theMap),  
  windowTitle(0), height(_height), width(_width), xLoc(_xLoc), yLoc(_yLoc),
  count(-1), theOutputFileName(0), 
  theDevice(0),
  vrp(3), vuv(3), vpn(3), cop(3), ViewMat(4,4), 
  projectionMode(0), vpWindow(4), ProjMat(4,4),
  portWindow(4)
{
  // set the WindowDevices title, height, wdth, xLoc and yLoc
  windowTitle = new char [strlen(_title)+1];
  strcpy(windowTitle, _title);

  fillMode = WIRE_MODE;
  projectionMode = PARALLEL_MODE;
  portWindow(0) = -1.0; portWindow(1) = 1.0;
  portWindow(2) = -1.0;  portWindow(3) = 1.0;

  theDevice = new OpenGlDevice();
  if (bitmapFileName != 0) {
		opserr << "OpenGLRenderer:;OpenGlRenderer - feature to save image only to BMP removed\n";
  }
  theDevice->WINOPEN(_title, _xLoc, _yLoc, _width, _height);;
  
  theDevice->CLEAR();


  // open the file for  making the movie
  if (outputFileName != 0) {
    windowTitle = new char [strlen(_title)+1];
    theOutputFileName = new char [strlen(outputFileName)+1];

    strcpy(windowTitle, _title);
    strcpy(theOutputFileName, outputFileName);
    theFile.open(theOutputFileName, ios::out);
    if (theFile.bad()) {
      opserr << "WARNING - OpenGLRenderer::OpenGLRenderer() - could not open file: " << outputFileName << endln;
      theOutputFileName = 0;
    } else {
      theFile << windowTitle << endln;
      theFile << 1.0*xLoc << " " << 1.0*yLoc << " " << 1.0*width << " " << 1.0*height << endln;
    }
  }

  theDevice->CLEAR();
  glClearColor(1.0f,1.0f,1.0f,1.0f);

  glEnable(GL_DEPTH_TEST);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  glFlush();
  theDevice->ENDIMAGE();
}

OpenGLRenderer::~OpenGLRenderer()
{
  if (theDevice != 0)
    delete theDevice;

  if (windowTitle != 0)
    delete [] windowTitle;

  if (theOutputFileName != 0) {
    theFile.close();
    delete [] theOutputFileName;
  }
}


int 
OpenGLRenderer::clearImage(void)
{
  theDevice->CLEAR();
  glClearColor(1.0f,1.0f,1.0f,1.0f);

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glFlush();
 
#ifdef _UNIX
  theDevice->ENDIMAGE();
#endif

  return 0;
}

int 
OpenGLRenderer::saveImage(const char *imageName)
{
  return theDevice->saveImage(imageName, 0);
}


int 
OpenGLRenderer::startImage(void)
{

  theMap->startImage();

  theDevice->STARTIMAGE();

  if (fillMode == WIRE_MODE) {
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
  } else {
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  }

  // glEnable(GL_BLEND);
  // glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  // glEnable(GL_LINE_SMOOTH);
  // glEnable(GL_POINT_SMOOTH);
  // glEnable(GL_POLYGON_SMOOTH);

    //
    // set up the viewing transformation
    //
    static Vector u(3), v(3), n(3);

    // vpnCopy used to store vpn, if no vpn yet specified vpn = VRP-PRP
    static Vector vpnCopy(3);

    if (vpn(0) == 0.0 && vpn(1) == 0.0 && vpn(2) == 0.0) {
      vpnCopy(0) = cop(0) - vrp(0);
      vpnCopy(1) = cop(1) - vrp(1);
      vpnCopy(2) = cop(2) - vrp(2);
    } else {
      vpnCopy(0) = vpn(0);
      vpnCopy(1) = vpn(1);
      vpnCopy(2) = vpn(2);
    }

    for (int i=0; i<3; i++) {
	n(i) = vpnCopy(i);
	v(i) = vuv(i);
    }

    if (n.Normalize() != 0) {
	opserr << "View::update() - VPN cannot have zero length\n";
	return -1;
    }

    // u = v % n;
    u(0) = v(1)*n(2) - v(2)*n(1) ;
    u(1) = v(2)*n(0) - v(0)*n(2) ;
    u(2) = v(0)*n(1) - v(1)*n(0) ;

    if (u.Normalize() != 0) {
	opserr << "View::update() - VUV X VPN cannot have zero length\n";
	return -1;
    }

    //v = n % u;
    v(0) = n(1)*u(2) - n(2)*u(1) ;
    v(1) = n(2)*u(0) - n(0)*u(2) ;
    v(2) = n(0)*u(1) - n(1)*u(0) ;

    v.Normalize();
    
    ViewMat(0,0) = u(0); ViewMat(0,1) = u(1); ViewMat(0,2) = u(2); ViewMat(0,3) = -(vrp^u);
    ViewMat(1,0) = v(0); ViewMat(1,1) = v(1); ViewMat(1,2) = v(2); ViewMat(1,3) = -(vrp^v);
    ViewMat(2,0) = n(0); ViewMat(2,1) = n(1); ViewMat(2,2) = n(2); ViewMat(2,3) = -(vrp^n);
    ViewMat(3,0) =  0.0; ViewMat(3,1) =  0.0; ViewMat(3,2) = 0.0;  ViewMat(3,3) =  1.0;

    for (int j=0; j<4; j++)
	for (int k=0; k<4; k++)
	    viewData[j+k*4] = ViewMat(j,k);
    
    glMatrixMode(GL_MODELVIEW);
    glLoadMatrixf(viewData);

   //
   // set up the projection transformation
   //

    float midU, midV, dopU, dopV, dopN, shU, shV, F, B;
    midU = (vpWindow(0) + vpWindow(1))/2;
    midV = (vpWindow(2) + vpWindow(3))/2;

    // PRP (COP) in viewing system
    float PRPu = u(0)*cop(0) + u(1)*cop(1) + u(2)*cop(2) + ViewMat(0,3);
    float PRPv = v(0)*cop(0) + v(1)*cop(1) + v(2)*cop(2) + ViewMat(1,3);
    float PRPn = n(0)*cop(0) + n(1)*cop(1) + n(2)*cop(2) + ViewMat(2,3);

    float diffU, diffV;// prpN;
	
    diffU = vpWindow(1)-vpWindow(0);
    diffV = vpWindow(3)-vpWindow(2);      

    dopU = midU - PRPu;
    dopV = midV - PRPv;
    dopN = -PRPn;

    shU = dopU/dopN;
    shV = dopV/dopN;
    F = clippingPlanes[0];
    B = clippingPlanes[1];

    /******* the equiv of the viewData transformation using glu **************
    glLoadIdentity();
    gluLookAt(cop[0),cop[1),cop[2),vrp[0),vrp[1),vrp[2),vuv[0),vuv[1),vuv[2));
    glTranslatef(0.0, 0.0, PRPn);  // note the PRPn transformation here
    **************************************************************************/

    if (projectionMode == PARALLEL_MODE) {

      ProjMat(0,0) = 2.0/diffU; ProjMat(0,1) = 0.0; ProjMat(0,2) = 2.0*shU/diffU; 
      ProjMat(0,3) = -2*midU/diffU,

	ProjMat(1,0) = 0.0; ProjMat(1,1) = 2/diffV; ProjMat(1,2) = 2*shV/diffV; 
	ProjMat(1,3) = -2*midV/diffV;

	ProjMat(2,0) =  0.0; ProjMat(2,1) =  0.0; ProjMat(2,2) = 1.0;  ProjMat(2,3) =  0.0;
	ProjMat(3,0) =  0.0; ProjMat(3,1) =  0.0; ProjMat(3,2) = 0.0;  ProjMat(3,3) =  1.0;
      
	for (int jj=0; jj<4; jj++)
	  for (int kk=0; kk<4; kk++)
	    projData[jj+kk*4] = ProjMat(jj,kk);

	glMatrixMode(GL_PROJECTION);
	glLoadMatrixf(projData);
	glOrtho(-1.0, 1.0, -1.0, 1.0, -F, -B);

    } else { // perspective projection

	double VRPn = -PRPn;  // VRP after T(-PRP)
	float near2 = PRPn-F;
	float far2  = PRPn-B;

	/**************** Having trouble with the following transformation ******
	float zMin = near2/far2;
	float a,b,c,e,f,g,h;

	a = 2.0*PRPn/(diffU * far);
	b = 2.0*PRPn/(diffV * far);
	c = -1.0/far;

	e = 1/(1-zMin);
	f = -zMin * e;
	
	g = PRPu - shU * PRPn;
	h = PRPv - shV * PRPn;

	ProjMat.Set(     a,   0.0,         0.0,    0.0,
		       0.0,     b,         0.0,    0.0,
		     a*shU, b*shV,         e*c,     -c,
		       -a*g,  -b*h, e*c*PRPn+f, -c*PRPn);  

	for (int jj=0; jj<4; jj++)
	  for (int kk=0; kk<4; kk++)
	    projData[jj*4+kk) = ProjMat.m(jj)[kk);

	glMatrixMode(GL_PROJECTION);
	glLoadMatrixf(projData);
	******* so in the meantime use the following - NO SHEARING  **************/

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(cop(0),cop(1),cop(2),vrp(0),vrp(1),vrp(2),vuv(0),vuv(1),vuv(2));

	glMatrixMode(GL_PROJECTION);
	float left = near2*diffU/VRPn/2;
	float bottom = near2*diffV/VRPn/2;
	float right = -left;
	float top = -bottom;
	glLoadIdentity();
	glFrustum(left,right,bottom,top,near2,far2);	

    }

    int xMin = floor((portWindow(0)+1) * width/2.0);
    int xDiff = floor((portWindow(1)-portWindow(0)) * width/2.0);
    int yMin = floor((portWindow(2)+1) * height/2.0);
    int yDiff = floor((portWindow(3)-portWindow(2)) * height/2.0);

    glViewport((GLsizei) xMin, (GLsizei) yMin, (GLsizei) xDiff, (GLsizei) yDiff);

    if (theOutputFileName != 0) {
	theFile << "StartImage\n";
	theFile << "VRP " << vrp(0) << " " << vrp(1) << " " 
	    << vrp(2) << " " << endln;
	theFile << "VPN " << vpn(0) << " " << vpn(1) << " " 
	    << vpn(2) << " " << endln;
	theFile << "VUV " << vuv(0) << " " << vuv(1) << " " 
	    << vuv(2) << " " << endln;
	theFile << "COP " << cop(0) << " " << cop(1) << " " 
	    << cop(2) << " " << endln;
	
	theFile << "PROJECTIONMODE " << 1.0*projectionMode << endln;
	theFile << "VPWINDOW " << vpWindow(0) << " " << vpWindow(1) << " "
	    << vpWindow(2) << " " << vpWindow(3) << " " << endln;
	theFile << "PLANES " << clippingPlanes[0] << " " << clippingPlanes[1] << "\n";
	theFile << "PORTWINDOW " << portWindow(0) << " " << portWindow(1) << " "
	    << portWindow(2) << " " << portWindow(3) << " " << endln;
    } 
 
    // done
    return 0;
}


int 
OpenGLRenderer::doneImage(void)
{
  if (theOutputFileName != 0) {
    theFile << "DoneImage\n";
  }
    
  if (count != -1) {
    count++;
  }

  theDevice->ENDIMAGE();    
  return 0;
}

int 
OpenGLRenderer::drawPoint(const Vector &pos1, float V1, int tag, int mode, int numPixels)
{
    glPointSize(numPixels);

    glBegin(GL_POINTS);
    float r, g, b;

    theMap->getRGB(V1, r, g, b);
    glColor3f(r,g,b);
    glVertex3f(pos1(0),pos1(1),pos1(2));
    
    if (theOutputFileName != 0) {
	theFile << "Point\n" << pos1(0) << " " << pos1(1) << " " << pos1(2) 
	    << " " << r << " " << g << " " << b << " " << endln;
    }

    glEnd();

    return 0;  
}


int 
OpenGLRenderer::drawPoint(const Vector &pos1, const Vector &rgb, int tag, int mode, int numPixels)
{
    glPointSize(numPixels);

    glBegin(GL_POINTS);
    float r, g, b;
    r = rgb(0);
    g = rgb(1);
    b = rgb(2);

    glColor3f(r,g,b);
    glVertex3f(pos1(0),pos1(1),pos1(2));
    
    if (theOutputFileName != 0) {
	theFile << "Point\n" << pos1(0) << " " << pos1(1) << " " << pos1(2) 
	    << " " << r << " " << g << " " << b << " " << endln;
    }

    glEnd();

    return 0;  
}



int 
OpenGLRenderer::drawLine(const Vector &pos1, const Vector &pos2, 
			 float V1, float V2, int tag, int mode, int width, int style)
{
    // open gl does a divide by zero error if points are the same - so check
    if (pos1(0) == pos2(0) && pos1(1) == pos2(1) && pos1(2) == pos2(2))
      return 0;

    width = 2; // S. Mazzoni
    glLineWidth(width);

    glBegin(GL_LINES);
    float r, g, b;

    theMap->getRGB(V1, r, g, b);
    glColor3f(r,g,b);

    glVertex3f(pos1(0),pos1(1),pos1(2));
    
    if (theOutputFileName != 0) {
	theFile << "Line\n" << pos1(0) << " " << pos1(1) << " " << pos1(2) 
	    << " " << r << " " << g << " " << b << " " << endln;
    }

    theMap->getRGB(V2, r, g, b);
    glColor3f(r,g,b);

    glVertex3f(pos2(0),pos2(1),pos2(2));
    glEnd();

    if (theOutputFileName != 0) {
	theFile << pos2(0) << " " << pos2(1) << " " << pos2(2) << " " << r 
	    << " " << g << " " << b << " " << endln;
    }

    return 0;  
}



int 
OpenGLRenderer::drawLine(const Vector &end1, const Vector &end2, 
			 const Vector &rgb1, const Vector &rgb2,
			 int tag, int mode, int width, int style)
{
    // open gl does a divide by zero error if points are the same
    if (end1(0) == end2(0) && end1(1) == end2(1) && end1(2) == end2(2))
      return 0;

    width = 2; // S. Mazzoni
    glLineWidth(width);

    glBegin(GL_LINES);
    float r, g, b;
    r = rgb1(0);
    g = rgb1(1);
    b = rgb1(2);
  
    if (theOutputFileName != 0) {
	theFile << "Line\n" << end1(0) << " " << end1(1) << " " << end1(2) 
	    << " " << r << " " << g << " " << b << " " << endln;
    }
    glColor3f(r,g,b);
    glVertex3f(end1(0),end1(1),end1(2));

    r = rgb2(0);
    g = rgb2(1);
    b = rgb2(2);
  
    if (theOutputFileName != 0) {
	theFile << end2(0) << " " << end2(1) << " " << end2(2) << " " << r 
	    << " " << g << " " << b << " " << endln;
    }
    glColor3f(r,g,b);

    glVertex3f(end2(0),end2(1),end2(2));
    glEnd();

    return 0;
}



int 
OpenGLRenderer::drawPolygon(const Matrix &pos, const Vector &data, int tag, int mode)
{
#ifdef _G3DEBUG
  if (pos.noCols() != 3) {
    opserr << "OpenGLRenderer::drawPolygon - matrix needs 3 cols\n";
    return -1;
  }
  if (pos.noRows() != data.Size()) {
    opserr << "OpenGLRenderer::drawPolygon - matrix & vector incompatable\n";
    return -1;
  }
#endif

  double posX, posY, posZ, value;
  float r,g,b;

    glBegin(GL_POLYGON);
    int numRows = pos.noRows();
    for (int i=0; i<numRows; i++) {
      posX = pos(i,0);
      posY = pos(i,1);
      posZ = pos(i,2);
      value = data(i);
      theMap->getRGB(value, r, g, b);

    if (theOutputFileName != 0) {
	theFile << posX << " " << posY << " " << posZ << " " << r 
	    << " " << g << " " << b << " " << endln;
    }
      glColor3f(r,g,b);
      glVertex3f(posX, posY, posZ);
    }

    glEnd();

    return 0;
}


int 
OpenGLRenderer::drawPolygon(const Matrix &pos, const Matrix &rgbData, int tag, int mode)

{
#ifdef _G3DEBUG
  if (pos.noCols() != 3 || rgbData.noCols() != 3) {
    opserr << "OpenGLRenderer::drawPolygon - matrix needs 3 cols\n";
    return -1;
  }
  if (pos.noRows() != rgbData.noRows()) {
    opserr << "OpenGLRenderer::drawPolygon - matrix & vector incompatable\n";
    return -1;
  }
#endif

  double posX, posY, posZ;
  float r,g,b;

    glBegin(GL_POLYGON);
    int numRows = pos.noRows();
    for (int i=0; i<numRows; i++) {
      posX = pos(i,0);
      posY = pos(i,1);
      posZ = pos(i,2);
      r = rgbData(i,0);
      g = rgbData(i,1);
      b = rgbData(i,2);

    if (theOutputFileName != 0) {
	theFile << posX << " " << posY << " " << posZ << " " << r 
	    << " " << g << " " << b << " " << endln;
    }
      glColor3f(r,g,b);
      glVertex3f(posX, posY, posZ);
    }

    glEnd();

    return 0;
}



int 
OpenGLRenderer::drawText(const Vector &pos, char *text, int length,
			 char horizontalJustify, char verticalJustify)
{
    // add POINTs to the FACE  
    int size = pos.Size();
    float x,y,z;
    if (size == 1) {
	x = pos(0);
	y = 0;
	z = 0;
    } else if (size == 2) {
	x = pos(0);
	y = pos(1);
	z = 0;
    } else {
	x = pos(0);
	y = pos(1);
	z = pos(2);
    }  

    theDevice->drawText(x,y,z, text, length, horizontalJustify, verticalJustify);

    return 0;
}

int 
OpenGLRenderer::setVRP(float x, float y, float z)
{
  vrp(0) = x;
  vrp(1) = y;
  vrp(2) = z;

  return 0;
}

int 
OpenGLRenderer::setVPN(float x, float y, float z)
{
  vpn(0) = x;
  vpn(1) = y;
  vpn(2) = z;

  return 0;
}

int 
OpenGLRenderer::setVUP(float x, float y, float z)
{
  vuv(0) = x;
  vuv(1) = y;
  vuv(2) = z;

  return 0;
}

int 
OpenGLRenderer::setViewWindow(float umin, float umax, float vmin, float vmax)
{
  if (umin > umax || vmin > vmax) {
      opserr << "OpenGLRenderer::setViewWindow() - invalid window ";
      opserr << umin << " "<< umax << " "<< vmin << " "<< vmax << endln;
      return -1;
  }

  vpWindow(0) = umin;
  vpWindow(1) = umax;
  vpWindow(2) = vmin;
  vpWindow(3) = vmax;

  return 0;
}

int 
OpenGLRenderer::setPlaneDist(float anear, float afar) 
{
  if ((anear < afar)) {
      opserr << "OpenGLRenderer::setClippingPlanes() - invalid planes";
      opserr << anear << " " << afar << endln;
      return -1;
  }

  clippingPlanes[0] = anear;
  clippingPlanes[1] = afar;

  return 0;
}

int 
OpenGLRenderer::setProjectionMode(const char *newMode)
{
  if ((strcmp(newMode, "parallel") == 0) || (strcmp(newMode, "Parallel") == 0))
    projectionMode = PARALLEL_MODE;
  else if ((strcmp(newMode, "perspective") == 0) || (strcmp(newMode, "Perspective") == 0))
    projectionMode = PERSPECTIVE_MODE;
  return 0;
}

int 
OpenGLRenderer::setFillMode(const char *newMode)
{
  if ((strcmp(newMode, "wire") == 0) || (strcmp(newMode, "Wire") == 0))
    fillMode = WIRE_MODE;
  else if ((strcmp(newMode, "fill") == 0) || (strcmp(newMode, "Fill") == 0))
    fillMode = FILL_MODE;

  return 0;
}

// eye location
int 
OpenGLRenderer::setPRP(float u, float v, float n){
  cop(0) = u;
  cop(1) = v;
  cop(2) = n;

  return 0;
}
    
int 
OpenGLRenderer::setPortWindow(float left, float right, 
			     float bottom, float top)
{
  if (left < -1 || right > 1 || bottom < -1 || top > 1
      || left > right || bottom > top) {
      
      opserr << "OpenGLRenderer::setPortWindow() - bounds invalid ";
      opserr << left << " "<< right << " "<< bottom << " "<< top << endln;
      return -1;
  }

  portWindow(0) = left;
  portWindow(1) = right;
  portWindow(2) = bottom;
  portWindow(3) = top;

  return 0;
}

