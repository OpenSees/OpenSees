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
// $Date: 2000-09-15 08:23:25 $
// $Source: /usr/local/cvs/OpenSees/SRC/renderer/Clipping.cpp,v $
                                                                        
                                                                        
#include "Clipping.h"
#include "Projection.h"

#ifndef TRUE
#define TRUE 0
#endif

#ifndef FALSE
#define FALSE 1
#endif

Clipping::Clipping()
{

}

Clipping::~Clipping()
{

}

int
Clipping::update(void){
  X = 1.0;
  Y = 1.0;
  Zfar = -1.0;
  Znear = 0.0;
  return 0;
}

FACE &CLIP_X(FACE &, float, float, float, float); // creates new, deletes old
FACE &CLIP_Y(FACE &, float, float, float, float); // creates new, deletes old
FACE &CLIP_Z(FACE &, float, float, float, float); // creates new, deletes old

FACE &
Clipping::transform(FACE &input)
{
  // cerr << "Clipping:1\n" << input;
  FACE &clippedX = CLIP_X(input,X,Y,Zfar,Znear);
  // cerr << "Clipping:2\n" << clippedX;
  FACE &clippedY = CLIP_Y(clippedX,X,Y,Zfar,Znear);
  // cerr << "Clipping:3\n" << clippedY;
  FACE &clippedZ = CLIP_Z(clippedY,X,Y,Zfar,Znear);
  // cerr << "Clipping:4\n" << clippedZ;
  return clippedZ;
}


MYPOINT *
Clipping::transformP(MYPOINT *input)
{
    if (input != 0) {
	float x = input->p[0];
	float y = input->p[1];
	float z = input->p[2];
	
	if (x < -X || x > X || y < -Y || Y > Y ||
	    z > Znear || z < Zfar) {
	    
	    delete input;
	    return 0;
	}
    }
    return input;
}


FACE &
CLIP_X(FACE &input, float X, float Y, float Zfar, float Znear)
{
// this first procedure is used to clip against pos and
  FACE  *output;
  MYPOINT *Point1, *Point2, *P1, *P2;
  float x1,y1,z1,x2,y2,z2;
  float sx, sy, sz, t;
  float svwx,svwy,svwz, vwx1,vwx2,vwy1,vwy2,vwz1,vwz2; // point viewing stuff
  int out1, out1a, out2;
 // float snx,sny,snz, nx1,nx2,ny1,ny2,nz1,nz2;          // the point normals
  float r1,r2, g1,g2, b1, b2, sr, sg, sb;

  // set up the new face
  output = new FACE();

  // check for quick return
  MYPOINT *point = input.pointList.GetHead();
  if (point == 0) {  
      delete &input;
      return *output;      
  }
  
  // first clip against the pos X axis, remember I have a cube so it's easy
  FOR_EACH_PAIR(Point1,Point2,input.pointList,{
      x1 = Point1->X();    x2 = Point2->X();
      y1 = Point1->Y();    y2 = Point2->Y();
      z1 = Point1->Z();    z2 = Point2->Z();

      r1 = Point1->r;     r2 = Point2->r;
      g1 = Point1->g;     g2 = Point2->g;
      b1 = Point1->b;     b2 = Point2->b;

      vwx1 = Point1->Viewp[0]; vwx2 = Point2->Viewp[0];
      vwy1 = Point1->Viewp[1]; vwy2 = Point2->Viewp[1];
      vwz1 = Point1->Viewp[2]; vwz2 = Point2->Viewp[2];

       
      // first determine in or out
      if (x1 > X) 
	  out1 = TRUE;
      else
	  out1 = FALSE;
      if (x2 > X)
	  out2 = TRUE;
      else
	  out2 = FALSE;
      
      if (out1 == FALSE || out2 == FALSE) {
	  // clip if either of points is in against pos X
	      
	  sx = x2-x1;
	  sy = y2-y1;
	  sz = z2-z1;
      
	  sr = r2 - r1;  sg = g2 - g1;  sb = b2 - b1;
	  svwx = vwx2 - vwx1; svwy = vwy2 - vwy1; svwz = vwz2 - vwz1;

	  if (out1 == TRUE) {
	      t = (X-x1)/sx;
	      x1 = x1 +t*sx;
	      y1 = y1 +t*sy;
	      z1 = z1 +t*sz;
	
	      r1 = r1 + t*sr;
	      g1 = g1 + t*sg;
	      b1 = b1 + t*sb;

	      vwx1 = vwx1 + t*svwx;
	      vwy1 = vwy1 + t*svwy;
	      vwz1 = vwz1 + t*svwz;
	      
	  }
	  if (out2 == TRUE) {
	      t = (X-x1)/sx;
	      x2 = x1 +t*sx;
	      y2 = y1 +t*sy;
	      z2 = z1 +t*sz;
		  
	      r2 = r1 + t*sr;
	      g2 = g1 + t*sg;
	      b2 = b1 + t*sb;

	      vwx2 = vwx1 + t*svwx;
	      vwy2 = vwy1 + t*svwy;
	      vwz2 = vwz1 + t*svwz;
	  }          

	  // now clip against -X
	  // first determine in or out
	  out1a = out1;
	  if (x1 < -X) 
	      out1 = TRUE;
	  else
	      out1 = FALSE;
	  if (x2 < -X)
	      out2 = TRUE;
	  else
	      out2 = FALSE;
      
      // clip if either in
      if (out1 == FALSE || out2 == FALSE) {
  
	sx = x2-x1;
	sy = y2-y1;
	sz = z2-z1;

	sr = r2 - r1;  sg = g2 - g1;  sb = b2 - b1;
	svwx = vwx2 - vwx1; svwy = vwy2 - vwy1; svwz = vwz2 - vwz1;
	
	if (out1 == TRUE) {
	  t = (-X-x1)/sx;
	  x1 = x1 +t*sx;
	  y1 = y1 +t*sy;
	  z1 = z1 +t*sz;
	  
	  r1 = r1 + t*sr;
	  g1 = g1 + t*sg;
	  b1 = b1 + t*sb;

	  vwx1 = vwx1 + t*svwx;
	  vwy1 = vwy1 + t*svwy;
	  vwz1 = vwz1 + t*svwz;
	  
	}
	if (out2 == TRUE) {
	  t = (-X-x1)/sx;
	  x2 = x1 +t*sx;
	  y2 = y1 +t*sy;
	  z2 = z1 +t*sz;
	  
	  r2 = r1 + t*sr;
	  g2 = g1 + t*sg;
	  b2 = b1 + t*sb;

	  vwx2 = vwx1 + t*svwx;
	  vwy2 = vwy1 + t*svwy;
	  vwz2 = vwz1 + t*svwz;

	}          
	
	// only add new P1 if P1 not inside both
	if ((out1 == TRUE && out1a == FALSE) ||
	    (out1 == FALSE && out1a == TRUE))  {  
	  P1 = new MYPOINT(0,x1,y1,z1,0,0,0,vwx1,vwy1,vwz1);
	  P1->r = r1; P1->b = b1; P1->g = g1;
	  output->AddPoint(*P1);
	}
	P2 = new MYPOINT(0,x2,y2,z2,0,0,0,vwx2,vwy2,vwz2);
	P2->r = r2; P2->b = b2; P2->g = g2;
	output->AddPoint(*P2);
      }  // else do nothing
    } // else do nothing again
  });

  delete &input;
  return *output;
}




FACE &
CLIP_Y(FACE &input, float X, float Y, float Zfar, float Znear)
{
  // this first procedure is used to clip against pos and
  FACE  *output;
  MYPOINT *Point1, *Point2, *P1, *P2;
  float x1,y1,z1,x2,y2,z2;
  float sx, sy, sz, t;
  int out1, out1a, out2;

  float svwx,svwy,svwz, vwx1,vwx2,vwy1,vwy2,vwz1,vwz2; // point viewing stuff
  float r1,r2, g1,g2, b1, b2, sr, sg, sb;

  output = new FACE();

  // check for quick return
  MYPOINT *point = input.pointList.GetHead();
  if (point == 0) {  
      delete &input;
      return *output;      
  }
  
  
  
  // first clip against the pos X axis, remember I have a cube so it's easy
  FOR_EACH_PAIR(Point1,Point2,input.pointList,{
    x1 = Point1->X();
    x2 = Point2->X();
    y1 = Point1->Y();
    y2 = Point2->Y();
    z1 = Point1->Z();
    z2 = Point2->Z();

    vwx1 = Point1->Viewp[0]; vwx2 = Point2->Viewp[0];
    vwy1 = Point1->Viewp[1]; vwy2 = Point2->Viewp[1];
    vwz1 = Point1->Viewp[2]; vwz2 = Point2->Viewp[2];

    r1 = Point1->r;     r2 = Point2->r;
    g1 = Point1->g;     g2 = Point2->g;
    b1 = Point1->b;     b2 = Point2->b;


    // first determine in or out
    if (y1 > Y) 
      out1 = TRUE;
    else
      out1 = FALSE;
    if (y2 > Y)
      out2 = TRUE;
    else
      out2 = FALSE;

    if (out1 == FALSE || out2 == FALSE) {
      // clip if either of points is in against pos X

      sx = x2-x1;
      sy = y2-y1;
      sz = z2-z1;

      sr = r2 - r1;  sg = g2 - g1;  sb = b2 - b1;
      svwx = vwx2 - vwx1; svwy = vwy2 - vwy1; svwz = vwz2 - vwz1;

      if (out1 == TRUE) {
	t = (Y-y1)/sy;
	x1 = x1 +t*sx;
	y1 = y1 +t*sy;
	z1 = z1 +t*sz;

	r1 = r1 + t*sr;
	g1 = g1 + t*sg;
	b1 = b1 + t*sb;

	vwx1 = vwx1 + t*svwx;
	vwy1 = vwy1 + t*svwy;
	vwz1 = vwz1 + t*svwz;

      }
      if (out2 == TRUE) {
	t = (Y-y1)/sy;
	x2 = x1 +t*sx;
	y2 = y1 +t*sy;
	z2 = z1 +t*sz;
	
	r2 = r1 + t*sr;
	g2 = g1 + t*sg;
	b2 = b1 + t*sb;

	vwx2 = vwx1 + t*svwx;
	vwy2 = vwy1 + t*svwy;
	vwz2 = vwz1 + t*svwz;

      }          

      // now clip against -Y

      // first determine in or out
      out1a = out1;
      if (y1 < -Y) 
	out1 = TRUE;
      else
	out1 = FALSE;
      if (y2 < -Y)
	out2 = TRUE;
      else
	out2 = FALSE;

      // clip if either in
      if (out1 == FALSE || out2 == FALSE) {
	
	sx = x2-x1;
	sy = y2-y1;
	sz = z2-z1;

	sr = r2 - r1;  sg = g2 - g1;  sb = b2 - b1;	
	svwx = vwx2 - vwx1; svwy = vwy2 - vwy1; svwz = vwz2 - vwz1;
	
	if (out1 == TRUE) {
	  t = (-Y-y1)/sy;
	  x1 = x1 +t*sx;
	  y1 = y1 +t*sy;
	  z1 = z1 +t*sz;

	  r1 = r1 + t*sr;
	  g1 = g1 + t*sg;
	  b1 = b1 + t*sb;

	  vwx1 = vwx1 + t*svwx;
	  vwy1 = vwy1 + t*svwy;
	  vwz1 = vwz1 + t*svwz;

	}
	if (out2 == TRUE) {
	  t = (-Y-y1)/sy;
	  x2 = x1 +t*sx;
	  y2 = y1 +t*sy;
	  z2 = z1 +t*sz;
	  
	  r2 = r1 + t*sr;
	  g2 = g1 + t*sg;
	  b2 = b1 + t*sb;

	  vwx2 = vwx1 + t*svwx;
	  vwy2 = vwy1 + t*svwy;
	  vwz2 = vwz1 + t*svwz;
	  
	}          

	// only add new P1 if P1 not inside both
	if ((out1 == TRUE && out1a == FALSE) ||
	    (out1 == FALSE && out1a == TRUE))  {  
	  P1 = new MYPOINT(0,x1,y1,z1,0,0,0,vwx1,vwy1,vwz1);
	  P1->r = r1; P1->b = b1; P1->g = g1;
	  output->AddPoint(*P1);
	}
	P2 = new MYPOINT(0,x2,y2,z2,0,0,0,vwx2,vwy2,vwz2);
	P2->r = r2; P2->b = b2; P2->g = g2;
	output->AddPoint(*P2);
      }  // else do nothing
    } // else do nothing again
  });

  delete &input;
  return *output;
}


FACE &
CLIP_Z(FACE &input, float X, float Y, float Zfar, float Znear)
{
  // this first procedure is used to clip against pos and
  FACE  *output;
  MYPOINT *Point1, *Point2, *P1, *P2;
  float x1,y1,z1,x2,y2,z2;
  float sx, sy, sz, t;
  int out1, out1a, out2;

  float svwx,svwy,svwz, vwx1,vwx2,vwy1,vwy2,vwz1,vwz2; 
  float sr, sg, sb, r1,r2, g1,g2, b1, b2;// r,g,b;

  // set up the new face
  output = new FACE();

  // check for quick return
  MYPOINT *point = input.pointList.GetHead();
  if (point == 0) {  
      delete &input;
      return *output;      
  }  
  
  // first clip against the pos X axis, remember I have a cube so it's easy
  FOR_EACH_PAIR(Point1,Point2,input.pointList,{
    x1 = Point1->X();
    x2 = Point2->X();
    y1 = Point1->Y();
    y2 = Point2->Y();
    z1 = Point1->Z();
    z2 = Point2->Z();

    r1 = Point1->r;     r2 = Point2->r;
    g1 = Point1->g;     g2 = Point2->g;
    b1 = Point1->b;     b2 = Point2->b;
    

    vwx1 = Point1->Viewp[0]; vwx2 = Point2->Viewp[0];
    vwy1 = Point1->Viewp[1]; vwy2 = Point2->Viewp[1];
    vwz1 = Point1->Viewp[2]; vwz2 = Point2->Viewp[2];

    r1 = Point1->r;     r2 = Point2->r;
    g1 = Point1->g;     g2 = Point2->g;
    b1 = Point1->b;     b2 = Point2->b;

    // first determine in or out
    if (z1 < Zfar) 
      out1 = TRUE;
    else
      out1 = FALSE;
    if (z2 < Zfar)
      out2 = TRUE;
    else
      out2 = FALSE;
    
    if (out1 == FALSE || out2 == FALSE) {
      // clip if either of points is in against Zfar
      
      sx = x2-x1;
      sy = y2-y1;
      sz = z2-z1;

      sr = r2 - r1;  sg = g2 - g1;  sb = b2 - b1;	      
      svwx = vwx2 - vwx1; svwy = vwy2 - vwy1; svwz = vwz2 - vwz1;
      
      if (out1 == TRUE) {
	t = (Zfar-z1)/sz;
	x1 = x1 +t*sx;
	y1 = y1 +t*sy;
	z1 = z1 +t*sz;

	r1 = r1 + t*sr;
	g1 = g1 + t*sg;
	b1 = b1 + t*sb;
	
	vwx1 = vwx1 + t*svwx;
	vwy1 = vwy1 + t*svwy;
	vwz1 = vwz1 + t*svwz;

      }
      if (out2 == TRUE) {
	t = (Zfar-z1)/sz;
	x2 = x1 +t*sx;
	y2 = y1 +t*sy;
	z2 = z1 +t*sz;

	r2 = r1 + t*sr;
	g2 = g1 + t*sg;
	b2 = b1 + t*sb;

	vwx2 = vwx1 + t*svwx;
	vwy2 = vwy1 + t*svwy;
	vwz2 = vwz1 + t*svwz;
	
      }          

      // now clip against Znear
      // first determine in or out
      out1a = out1;
      if (z1 > Znear) 
	out1 = TRUE;
      else
	out1 = FALSE;
      if (z2 > Znear)
	out2 = TRUE;
      else
	out2 = FALSE;


      // clip if either in
      if (out1 == FALSE || out2 == FALSE) {
  
	sx = x2-x1;
	sy = y2-y1;
	sz = z2-z1;
	
	sr = r2 - r1;  sg = g2 - g1;  sb = b2 - b1;	      
	svwx = vwx2 - vwx1; svwy = vwy2 - vwy1; svwz = vwz2 - vwz1;
	
	if (out1 == TRUE) {
	  t = (Znear-z1)/sz;
	  x1 = x1 +t*sx;
	  y1 = y1 +t*sy;
	  z1 = z1 +t*sz;

	  r1 = r1 + t*sr;
	  g1 = g1 + t*sg;
	  b1 = b1 + t*sb;	  

	  vwx1 = vwx1 + t*svwx;
	  vwy1 = vwy1 + t*svwy;
	  vwz1 = vwz1 + t*svwz;

	}
	if (out2 == TRUE) {
	  t = (Znear-z1)/sz;
	  x2 = x1 +t*sx;
	  y2 = y1 +t*sy;
	  z2 = z1 +t*sz;
	  
	  r2 = r1 + t*sr;
	  g2 = g1 + t*sg;
	  b2 = b1 + t*sb;

	  vwx2 = vwx1 + t*svwx;
	  vwy2 = vwy1 + t*svwy;
	  vwz2 = vwz1 + t*svwz;
	  
	}          
	
	if ((out1 == TRUE && out1a == FALSE) ||
	    (out1 == FALSE && out1a == TRUE))  {  


	  // only add new P1 if P1 not inside both
	  P1 = new MYPOINT(0,x1,y1,z1,0,0,0,vwx1,vwy1,vwz1);
	  P1->r = r1; P1->b = b1; P1->g = g1;
	  output->AddPoint(*P1);
	}
	P2 = new MYPOINT(0,x2,y2,z2,0,0,0,vwx2,vwy2,vwz2);
	P2->r = r1; P2->b = b1; P2->g = g1;
	output->AddPoint(*P2);
      }  // else do nothing
    } // else do nothing again
  });

  delete &input;
  return *output;
}

