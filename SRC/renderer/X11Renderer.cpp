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
                                                                        
// $Revision: 1.2 $
// $Date: 2000-12-12 07:11:52 $
// $Source: /usr/local/cvs/OpenSees/SRC/renderer/X11Renderer.cpp,v $
                                                                        
                                                                        
// File: ~/renderer/X11Renderer.C
//
// Written: fmk 
// Created: 10/98
// Revision: A
//
// Description: This file contains the class definition for X11Renderer.
// X11Renderer is an class which diplays using X11 or openGL.
//
// What: "@(#) X11Renderer.h, revA"

#include <X11Renderer.h>
#include <ColorMap.h>
#include <string.h>

#include <db.H>

#include <Matrix.h>
#include <Vector.h>
#include <View.h>
#include <Projection.h>
#include <Viewport.h>
#include <Clipping.h>
#include <WindowDevice.h>
#include <Scan.h>

X11Renderer::X11Renderer(char *title, int xLoc, int yLoc, int width, int height,
			 ColorMap &_theMap)
  :Renderer(_theMap)
{
    theView = new View;
    theProjection = new Projection;
    theClipping = new Clipping;
    theViewport = new Viewport;
    theScan = new ScanLineConverter;
    theDevice = new WindowDevice;
    
    theScan->setDevice(*theDevice);
    theScan->setProjection(*theProjection);  
    theViewport->setDevice(*theDevice);
  
    theDevice->WINOPEN(title, xLoc, yLoc, width, height);
    theScan->setFillMode(0);    

    aFile = 0;
}


X11Renderer::X11Renderer(char *title, int xLoc, int yLoc, int width, int height,
			 ColorMap &_theMap,
			 char *fileName)
  :Renderer(_theMap)
{
    theView = new View;
    theProjection = new Projection;
    theClipping = new Clipping;
    theViewport = new Viewport;
    theScan = new ScanLineConverter;
    theDevice = new WindowDevice;
    
    theScan->setDevice(*theDevice);
    theScan->setProjection(*theProjection);  
    theViewport->setDevice(*theDevice);
  
    theDevice->WINOPEN(title, xLoc, yLoc, width, height);
    theScan->setFillMode(0);    

    // if file is specified, copy name and open the file
    if (fileName != 0) {
      if (strlen(fileName) > MAX_FILENAMELENGTH) 
	g3ErrorHandler->warning("WARNING - X11Renderer::X11Renderer() - fileName %s too long, max %d\n",
			      fileName, MAX_FILENAMELENGTH);
      else {
	strcpy(theFileName, fileName);    
	theFile.open(fileName, ios::out);
	if (!theFile) {
	  g3ErrorHandler->warning("WARNING - X11Renderer::X11Renderer() - could not open file %s\n",fileName);	  
	  aFile = 0;
	} else {
	  aFile = 1;
	  theFile << title << endl;
	  theFile << xLoc << " " << yLoc << " " << width << " " << height << endl;
	}
      } 
    }
}

X11Renderer::~X11Renderer()
{
    delete    theView;
    delete    theProjection;
    delete    theClipping;
    delete    theViewport;
    delete    theScan;
    delete    theDevice;
    
    if (aFile == 1)
      theFile.close();
}


int 
X11Renderer::clearImage(void)
{
  theDevice->CLEAR();
  return 0;
}

int 
X11Renderer::startImage(void)
{
  theView->update();
  theProjection->update();
  theClipping->update();
  theViewport->update();
  theScan->update();
  theDevice->STARTIMAGE();

  if (aFile == 1) {
    theFile << "StartImage\n";
    theFile << "VRP " << theView->vrp[0] << " " << theView->vrp[1] << " " << theView->vrp[2] << " " << endl;
    theFile << "VPN " << theView->vpn[0] << " " << theView->vpn[1] << " " << theView->vpn[2] << " " << endl;
    theFile << "VUV " << theView->vuv[0] << " " << theView->vuv[1] << " " << theView->vuv[2] << " " << endl;
    theFile << "COP " << theProjection->cop[0] << " " << theProjection->cop[1] << " " << theProjection->cop[2] << " " << endl;
    theFile << "PROJECTIONMODE " << theProjection->projection_mode << endl;
    theFile << "VPWINDOW " << theProjection->vpwindow[0] << " " << theProjection->vpwindow[1] << " "
	    << theProjection->vpwindow[2] << " " << theProjection->vpwindow[3] << " " << endl;
    theFile << "PLANES " << theProjection->planedist[0] << " " << theProjection->planedist[2] << " " << endl;
    theFile << "PORTWINDOW " << theViewport->portwindow[0] << " " << theViewport->portwindow[1] << " "
	    << theViewport->portwindow[2] << " " << theViewport->portwindow[3] << " " << endl;

  } 

  return 0;
}

int 
X11Renderer::doneImage(void)
{
  theDevice->ENDIMAGE();

  if (aFile == 1) {
    theFile << "DoneImage\n";
  }
  return 0;
}

int 
X11Renderer::drawLine(const Vector &pos1, const Vector &pos2, 
		       float V1, float V2)
{
  FACE *theFace = new FACE();	
  
  int size;
  float x,y,z, r, g, b;
  MYPOINT *point;

  // add POINTs to the FACE  
  size = pos1.Size();
  if (size == 1) {
    x = pos1(0);
    y = 0;
    z = 0;
  } else if (size == 2) {
    x = pos1(0);
    y = pos1(1);
    z = 0;
  } else {
    x = pos1(0);
    y = pos1(1);
    z = pos1(2);
  }  
  point = new MYPOINT(1,x,y,z);
  r = theMap->getRed(V1);
  g = theMap->getGreen(V1);
  b = theMap->getBlue(V1);
  point->r = r;
  point->g = g;
  point->b = b;
    
  theFace->AddPoint(*point);

  if (aFile == 1) {
    theFile << "Line\n" << x << " " << y << " " << z << " " << r << " " << g << " " << b << " " << endl;
  }

  size = pos2.Size();
  if (size == 1) {
    x = pos2(0);
    y = 0;
    z = 0;
  } else if (size == 2) {
    x = pos2(0);
    y = pos2(1);
    z = 0;
  } else {
    x = pos2(0);
    y = pos2(1);
    z = pos2(2);
  }  
  point = new MYPOINT(2,x,y,z);
  r = theMap->getRed(V2);
  g = theMap->getGreen(V2);
  b = theMap->getBlue(V2);
  point->r = r;
  point->g = g;
  point->b = b;
    
  theFace->AddPoint(*point);

  if (aFile == 1) {
    theFile << x << " " << y << " " << z << " " << r << " " << g << " " << b << " " << endl;
  }

  FACE &res1 = theView->transform(*theFace);
  // cerr << "X11Renderer: face after view " << theFace;      
  FACE &res2 = theProjection->transform(res1);
  // cerr << "X11Renderer: face after projection " << theFace;        
  FACE &res3 = theClipping->transform(res2);
  // cerr << "X11Renderer: face after clipping " << theFace;          
  FACE &res4 = theViewport->transform(res3);
  // cerr << "X11Renderer: face after viewport " << theFace;            
  theScan->scanLine(res4);
  return 0;  
}

int 
X11Renderer::drawLine(const Vector &pos1, const Vector &pos2, 
		      const Vector &rgb1, const Vector &rgb2)
{
  FACE *theFace = new FACE();	
  
  int size;
  float x,y,z, r, g, b;
  MYPOINT *point;

  // add POINTs to the FACE  
  size = pos1.Size();
  if (size == 1) {
    x = pos1(0);
    y = 0;
    z = 0;
  } else if (size == 2) {
    x = pos1(0);
    y = pos1(1);
    z = 0;
  } else {
    x = pos1(0);
    y = pos1(1);
    z = pos1(2);
  }  
  point = new MYPOINT(1,x,y,z);
  r = rgb1(0);
  g = rgb1(1);
  b = rgb1(2);
  point->r = r;
  point->g = g;
  point->b = b;
    
  theFace->AddPoint(*point);

  if (aFile == 1) {
    theFile << "Line\n" << x << " " << y << " " << z << " " << r << " " << g << " " << b << " " << endl;
  }

  size = pos2.Size();
  if (size == 1) {
    x = pos2(0);
    y = 0;
    z = 0;
  } else if (size == 2) {
    x = pos2(0);
    y = pos2(1);
    z = 0;
  } else {
    x = pos2(0);
    y = pos2(1);
    z = pos2(2);
  }  
  point = new MYPOINT(2,x,y,z);
  r = rgb2(0);
  g = rgb2(1);
  b = rgb2(2);
  point->r = r;
  point->g = g;
  point->b = b;
    
  theFace->AddPoint(*point);

  if (aFile == 1) {
    theFile << x << " " << y << " " << z << " " << r << " " << g << " " << b << " " << endl;
  }

  FACE &res1 = theView->transform(*theFace);
  // cerr << "X11Renderer: face after view " << theFace;      
  FACE &res2 = theProjection->transform(res1);
  // cerr << "X11Renderer: face after projection " << theFace;        
  FACE &res3 = theClipping->transform(res2);
  // cerr << "X11Renderer: face after clipping " << theFace;          
  FACE &res4 = theViewport->transform(res3);
  // cerr << "X11Renderer: face after viewport " << theFace;            
  theScan->scanLine(res4);
  return 0;  
}





int 
X11Renderer::drawPolygon(const Matrix &pos, const Vector &data)

{
#ifdef _G3DEBUG
  if (pos.noCols() != 3) {
    g3ErrorHandler->warning("OpenGLRenderer::drawPolygon - matrix needs 3 cols\n");
    return -1;
  }
  if (pos.noRows() != data.Size()) {
    g3ErrorHandler->warning("OpenGLRenderer::drawPolygon - matrix & vector incompatable\n");
    return -1;
  }
#endif

  FACE *theFace = new FACE();	

  int size;
  float posX,posY,posZ, r, g, b;
  double value;
  MYPOINT *point;

  // add POINTs to the FACE  
  int numRows = pos.noRows();
  for (int i=0; i<numRows; i++) {
    posX = pos(i,0);
    posY = pos(i,1);
    posZ = pos(i,2);
    value = data(i);
    r = theMap->getRed(value);
    g = theMap->getGreen(value);
    b = theMap->getBlue(value);      

    if (aFile == 1) {
      theFile << posX << " " << posY << " " << posZ << " " << r 
	      << " " << g << " " << b << " " << endl;
    }	

    point = new MYPOINT(1,posX,posY,posZ);
    point->r = r;
    point->g = g;
    point->b = b;
    
    theFace->AddPoint(*point);
  }

    
  theFace->AddPoint(*point);

  // display the face
  FACE &res1 = theView->transform(*theFace);
  // cerr << "X11Renderer: face after view " << theFace;      
  FACE &res2 = theProjection->transform(res1);
  // cerr << "X11Renderer: face after projection " << theFace;        
  FACE &res3 = theClipping->transform(res2);
  // cerr << "X11Renderer: face after clipping " << theFace;          
  FACE &res4 = theViewport->transform(res3);
  // cerr << "X11Renderer: face after viewport " << theFace;            
  theScan->scanPolygon(res4);
  return 0;  
}



int 
X11Renderer::drawGText(const Vector &pos, char *text, int length)
{
  MYPOINT *point;

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
  point = new MYPOINT(1,x,y,z);

  MYPOINT *res =0;
  res = theView->transformP(point);
  if (res != 0) {
      res = theProjection->transformP(res);
  }
  if (res != 0) {
      res = theClipping->transformP(res);
  }

  if (res != 0) {
      res = theViewport->transformP(res);
  }
  if (res != 0) {
      float x = res->p[0];
      float y = res->p[1];
      theDevice->drawText(x,y, text, length);
  }
  
  if (res != 0)
      delete res;
  
  return 0;
}

int 
X11Renderer::drawLText(const Vector &pos, char *text, int length)
{
//  MYPOINT *point;

  // add POINTs to the FACE  
  int size = pos.Size();
  float x,y;
  if (size == 1) {
    x = pos(0);
    y = 0;
  } else if (size == 2) {
    x = pos(0);
    y = pos(1);
  } else {
    x = pos(0);
    y = pos(1);
  }  

  theDevice->drawText(x,y, text, length);
  
  return 0;
}




int
X11Renderer::displayFace(FACE &theFace)
{
  // cerr << "X11Renderer: input face " << theFace;    
  FACE &res1 = theView->transform(theFace);
  // cerr << "X11Renderer: face after view " << theFace;      
  FACE &res2 = theProjection->transform(res1);
  // cerr << "X11Renderer: face after projection " << theFace;        
  FACE &res3 = theClipping->transform(res2);
  // cerr << "X11Renderer: face after clipping " << theFace;          
  FACE &res4 = theViewport->transform(res3);
  // cerr << "X11Renderer: face after viewport " << theFace;            
  theScan->scanPolygon(res4);
  return 0;
}




int 
X11Renderer::setVRP(float x, float y, float z)
{
  theView->vrp[0] = x;
  theView->vrp[1] = y;
  theView->vrp[2] = z;

  return 0;
}

int 
X11Renderer::setVPN(float x, float y, float z)
{
  theView->vpn[0] = x;
  theView->vpn[1] = y;
  theView->vpn[2] = z;

  return 0;
}

int 
X11Renderer::setVUP(float x, float y, float z)
{
  theView->vuv[0] = x;
  theView->vuv[1] = y;
  theView->vuv[2] = z;

  return 0;
}

int 
X11Renderer::setViewWindow(float umin, float umax, float vmin, float vmax)
{
  if (umin > umax || vmin > vmax) {
      cerr << "X11Renderer::setViewWindow() - invalid window ";
      cerr << umin << " "<< umax << " "<< vmin << " "<< vmax << endl;
      return -1;
  }

  theProjection->vpwindow[0] = umin;
  theProjection->vpwindow[1] = umax;
  theProjection->vpwindow[2] = vmin;
  theProjection->vpwindow[3] = vmax;

  return 0;
}

int 
X11Renderer::setPlaneDist(float anear, float afar) 
{

   if ((anear < 0.0) || (afar < 0.0)) {
      cerr << "X11Renderer::setPlaneDist() - invalid planes";
      cerr << anear << " " << afar << endl;
      return -1;
  }
  
  theProjection->planedist[0] = anear;
  theProjection->planedist[2] = afar;

  return 0;
}

int 
X11Renderer::setProjectionMode(int newMode)
{
  theProjection->projection_mode = newMode;
  return 0;
}

int 
X11Renderer::setFillMode(int newMode)
{
  theScan->setFillMode(newMode);
  return 0;
}

// eye location
int 
X11Renderer::setPRP(float u, float v, float n){
  theProjection->cop[0] = u;
  theProjection->cop[1] = v;
  theProjection->cop[2] = n;

  return 0;
}
    
int 
X11Renderer::setPortWindow(float left, float right, 
			     float bottom, float top)
{
  if (left < -1 || right > 1 || bottom < -1 || top > 1
      || left > right || bottom > top) {
      
      cerr << "X11Renderer::setPortWindow() - bounds invalid ";
      cerr << left << " "<< right << " "<< bottom << " "<< top << endl;
      return -1;
  }

  theViewport->portwindow[0] = left;
  theViewport->portwindow[1] = right;
  theViewport->portwindow[2] = bottom;
  theViewport->portwindow[3] = top;

  return 0;
}

