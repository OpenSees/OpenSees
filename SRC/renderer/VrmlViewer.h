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
// $Source: /usr/local/cvs/OpenSees/SRC/renderer/VrmlViewer.h,v $
                                                                        
                                                                        
// File: ~/graphics/VrmlViewer.h
//
// Written: fmk 
// Created: 1/99
// Revision: A
//
// Description: This file contains the class definition for VrmlViewer.
// A VrmlViewer object is used to create an image of the domain, this image is 
// stored in a file. The language used to store the image is VRML.
//
// What: "@(#) VrmlViewer.h, revA"

#ifndef VrmlViewer_h
#define VrmlViewer_h

#include <Renderer.h>

class ColorMap;
class Domain;
class fstream;

class VrmlViewer : public Renderer
{
 public:
    VrmlViewer(char *fileName, Domain &theDomain, ColorMap &theMap);    
    
    ~VrmlViewer();

    int clearImage(void);    
    int doneImage(void);
    
    int drawLine(const Vector &, const Vector &, 
		 float V1, float V2);
    
    int drawTriangle(const Vector &, const Vector &, const Vector &,
		     float V1, float V2, float V3);


    int setVRP(float x, float y, float z); // point on view plane    
    int setVPN(float x, float y, float z); // view plane normal
    int setVUP(float x, float y, float z); // view-up vector
	
    // the following are in view coordinates	
    int setViewWindow(float, float, float, float); // view bounds
                               // umin, umax, vmin, vmax

    int setPlaneDist(float, float); // location of
                               // near and far clipping planes

    int setProjectionMode(int); // 
    int setFillMode(int);    // 1 = wire, otherwise fill
    
    int setPRP(float u, float v, float n); // eye location if 
	                         // perspective, dirn to +ViewPlane if parallel

    // the following are in normalized coordinates
    int setPortWindow(float, float, float, float); // view port
                              // left, right, bottom, top [-1,1,-1,1]


 protected:
    
 private:
    fstream *vrmlFile;
    char vrmlFileName[50];    
    ColorMap *theMap;
};

#endif

