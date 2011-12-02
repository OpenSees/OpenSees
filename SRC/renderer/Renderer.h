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
// $Source: /usr/local/cvs/OpenSees/SRC/renderer/Renderer.h,v $
                                                                        
                                                                        
// File: ~/renderer/Renderer.h
//
// Written: fmk 
// Created: 10/98
// Revision: A
//
// Description: This file contains the class definition for Renderer.
// Renderer is an abstract base class. An Renderer object is used
// to process 3d graphics.
//
// What: "@(#) Renderer.h, revA"

#ifndef Renderer_h
#define Renderer_h

class Vector;
class ColorMap;

class Renderer
{
  public:
    Renderer(ColorMap &theMap);    
    virtual ~Renderer();

    // method to set the color map
    void setColorMap(ColorMap &theMap);

    // method to clear the current image
    virtual int clearImage(void) =0;    

    // methods to be invoked when image processing is to start or is finished
    virtual int startImage(void) =0;    
    virtual int doneImage(void) =0;

    // methods invoked by the objects to display themselves    
    virtual int drawLine(const Vector &, const Vector &, 
			 float V1, float V2) =0;

    virtual int drawLine(const Vector &end1, const Vector &end2, 
			 const Vector &rgb1, const Vector &rgb2) =0;
    
    virtual int drawTriangle(const Vector &, const Vector &, const Vector &,
			     float V1, float V2, float V3) =0;

    virtual int drawVector(const Vector &position, const Vector &value);
    
    virtual int drawGText(const Vector &posGlobal, char *string, int length);    
    virtual int drawLText(const Vector &posLocal, char *string, int length); 


    // 
    // the following are for setting up the vieing system
    //

    // the following are in world coordinates & define view coord system
    virtual int setVRP(float x, float y, float z) =0; // point on view plane    
    virtual int setVPN(float x, float y, float z) =0; // view plane normal
    virtual int setVUP(float x, float y, float z) =0; // view-up vector
	
    // the following are in view coordinates	
    virtual int setViewWindow(float, float, float, float) =0; // view bounds
                               // umin, umax, vmin, vmax

    virtual int setPlaneDist(float, float) =0; // location of
                               // near and far clipping planes

    virtual int setProjectionMode(int) =0; // 
    virtual int setFillMode(int) =0;    // 1 = wire, otherwise fill
    
    virtual int setPRP(float u, float v, float n) =0; // eye location if 
	                         // perspective, dirn to +ViewPlane if parallel

    // the following are in normalized coordinates
    virtual int setPortWindow(float, float, float, float) =0; // view port
                              // left, right, bottom, top [-1,1,-1,1]


  protected:
    ColorMap *theMap;
    
  private:
};


#endif

