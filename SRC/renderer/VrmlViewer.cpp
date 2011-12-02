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
                                                                        
// $Revision: 1.4 $
// $Date: 2003-02-19 15:43:25 $
// $Source: /usr/local/cvs/OpenSees/SRC/renderer/VrmlViewer.cpp,v $
                                                                        
                                                                        
// File: ~/graphics/VrmlViewer.C
//
// Written: fmk 
// Created: 1/99
// Revision: A
//
// Description: This file contains the class definition for VrmlViewer.
// VrmlViewer is an class which diplays using X11 or openGL.
//
// What: "@(#) VrmlViewer.h, revA"

#include <VrmlViewer.h>
#include <ColorMap.h>

#include <OPS_Globals.h>
#include <iomanip>
using std::ios;

#include <string.h>

#include <stdlib.h>
#include <stdio.h>

#include <Vector.h>

VrmlViewer::VrmlViewer(char *fileName,
		       Domain &_theDomain, ColorMap &_theMap)
  :Renderer(_theMap)
{
    strcpy(vrmlFileName, fileName);
    vrmlFile = new fstream(vrmlFileName, ios::out);
    if (vrmlFile == 0) {
	opserr << "FATAL - VrmlViewer::VrmlViewer() - could not open file ";
	opserr << fileName << endln;
	exit(-1);
    }
  (*vrmlFile) << "#VRML V2.0 utf8 \n";        
}

VrmlViewer::~VrmlViewer()
{
    vrmlFile->close();
    delete vrmlFile;
}

int 
VrmlViewer::clearImage(void)
{
    // open the file again
    vrmlFile = new fstream(vrmlFileName, ios::out);
    if (vrmlFile == 0) {
	opserr << "FATAL - VrmlViewer::clearImage() - could not open file ";
	opserr << vrmlFileName << endln;
	return -1;
    }
  (*vrmlFile) << "#VRML V2.0 utf8 \n";    
    return 0;
}


int 
VrmlViewer::doneImage(void)
{
    vrmlFile->close();
    return 0;
}

int 
VrmlViewer::drawLine(const Vector &pos1, const Vector &pos2, 
		       float V1, float V2)
{
  float x,y,z, r, g, b;


    (*vrmlFile) << "Shape { geometry IndexedLineSet ";
  (*vrmlFile) << "{ coord Coordinate { \n\t\t point [\n";
  
  int size = pos1.Size();
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
  (*vrmlFile) << "\t\t\t " << x << "  " << y << "  " << z << ",\n";
  


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

  (*vrmlFile) << "\t\t\t " << x << "  " << y << "  " << z << " ] }\n";  
  (*vrmlFile) << "          coordIndex [ 0 1 ]\n";  
  (*vrmlFile) << "          colorPerVertex TRUE\n ";    
  (*vrmlFile) << "          color Color { \n\t color [\n";  
  
  r = theMap->getRed(V1);
  g = theMap->getGreen(V1);
  b = theMap->getBlue(V1);

  (*vrmlFile) << "\t\t\t " << r << "  " << g << "  " << b << ",\n";
  
  r = theMap->getRed(V2);
  g = theMap->getGreen(V2);
  b = theMap->getBlue(V2);

  (*vrmlFile) << "\t\t\t " << r << "  " << g << "  " << b << " ] \n }}}\n";
  
  return 0;
}





int 
VrmlViewer::drawTriangle(const Vector &pos1, const Vector &pos2,
			   const Vector &pos3,
			   float V1, float V2, float V3)
{
  int size;
  float x,y,z, r, g, b;

  (*vrmlFile) << "Shape { geometry IndexedFaceSet ";
  (*vrmlFile) << "{ coord Coordinate { \n\t\t point [\n";  
  
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

  (*vrmlFile) << "\t\t\t " << x << "  " << y << "  " << z << ",\n";

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

  (*vrmlFile) << "\t\t\t " << x << "  " << y << "  " << z << ",\n";
  
  size = pos3.Size();
  if (size == 1) {
    x = pos3(0);
    y = 0;
    z = 0;
  } else if (size == 2) {
    x = pos3(0);
    y = pos3(1);
    z = 0;
  } else {
    x = pos3(0);
    y = pos3(1);
    z = pos3(2);
  }  

  (*vrmlFile) << "\t\t\t " << x << "  " << y << "  " << z << " ] }\n";  
  (*vrmlFile) << "          coordIndex [ 0 1 2 ]\n";  
  (*vrmlFile) << "          colorPerVertex TRUE\n ";    
  (*vrmlFile) << "          color Color { \n\t color [\n";  
  
  r = theMap->getRed(V1);
  g = theMap->getGreen(V1);
  b = theMap->getBlue(V1);
  
  (*vrmlFile) << "\t\t\t " << r << "  " << g << "  " << b << ",\n";
      
  r = theMap->getRed(V2);
  g = theMap->getGreen(V2);
  b = theMap->getBlue(V2);
  
  (*vrmlFile) << "\t\t\t " << r << "  " << g << "  " << b << ",\n";  

  r = theMap->getRed(V3);
  g = theMap->getGreen(V3);
  b = theMap->getBlue(V3);

  (*vrmlFile) << "\t\t\t " << r << "  " << g << "  " << b << " ] \n }}}\n";
  
  return 0;
}


int 
VrmlViewer::setVRP(float x, float y, float z)
{
  return 0;
}


int 
VrmlViewer::setVPN(float x, float y, float z) 
{
  return 0;
}

int 
VrmlViewer::setVUP(float x, float y, float z) 
{
  return 0;
}

int 
VrmlViewer::setViewWindow(float, float, float, float) 
{
  return 0;
}

int 
VrmlViewer::setPlaneDist(float, float) 
{
  return 0;
}

int 
VrmlViewer::setProjectionMode(int) 
{
  return 0;
}

int 
VrmlViewer::setFillMode(int)    
{
  return 0;
}

int 
VrmlViewer::setPRP(float u, float v, float n) 
{
  return 0;
}

int 
VrmlViewer::setPortWindow(float, float, float, float)
{
  return 0;
}

