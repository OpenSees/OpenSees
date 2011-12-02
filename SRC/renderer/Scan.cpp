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
// $Date: 2003-02-14 23:01:58 $
// $Source: /usr/local/cvs/OpenSees/SRC/renderer/Scan.cpp,v $
                                                                        
                                                                        
#include <math.h>
#include <limits.h>
#include <stdlib.h>

#include "Scan.h"
#include "Device.h"
#include "Projection.h"

// #define DEBUG_SCAN

/*** NOTE:: device.ENDIMAGE(); *****/

ScanLineConverter::ScanLineConverter()
:zBuffer(0), sizeBuffer(0)
{
  scanMode = WINDING;
  fillMode = 1;
}

ScanLineConverter::~ScanLineConverter()
{
  if (zBuffer != 0)
      delete [] zBuffer;
}

void 
ScanLineConverter::setScanMode(int mode)
 {
  scanMode = mode;
}

void 
ScanLineConverter::setFillMode(int mode)
{
  fillMode = mode;
}


int
ScanLineConverter::update(void)
{
     // first set the zbuffer
     int i;// j;
     
     length = theDevice->GetHeight(); 
     width = theDevice->GetWidth();

     int size = length*width;
     if (size != sizeBuffer) {
	 if (zBuffer != 0)
	     delete [] zBuffer;
	 zBuffer = new float[size];
       if (zBuffer == 0) {
	 opserr << "ScanLineConverter::update() - out of memory\n";
	 return -1;
       }
       sizeBuffer = size;
     }

     for (i=0; i<length*width; i++)
       zBuffer[i] = -2; // we have clipped out all from 0 -2

    return 0;
}


void
ScanLineConverter::scanLine(FACE &face)
{
#ifdef _UNIX
    MYPOINT *point1 = face.pointList.GetHead();
    if (point1 != 0) {
	if (fillMode == 1) {
	    DrawWire(face);
	} else {
	    MYPOINT *point1 = face.pointList.SetToHead();	
      	MYPOINT *point2 = face.pointList.Forward();	

	    if (point2 != 0) {
		int x1 = (int)point1->X();
		int y1 = (int)point1->Y();    
		float zCur = point1->Z();    	    
		float rCur = point1->r;
		float gCur = point1->g;
		float bCur = point1->b;

		int x2 = (int)point2->X();
		int y2 = (int)point2->Y();    	
		float z2 = point2->Z();    		    
		float r2 = point2->r;
		float g2 = point2->g;
		float b2 = point2->b;	    
		/*
		rCur = 1.0; gCur = 0; bCur = 0;
		r2 = 0; g2 = 0; b2 = 1;
		*/
		int dx = x2 - x1;
		int dy = y2 - y1;
		float dz = z2 - zCur;
		float dr = r2 -rCur;
		float dg = g2 -gCur;
		float db = b2 -bCur;
		
		if (dx == 0 && dy == 0) {
		    ;
		} else if (dx == 0) {
		    int L = dy;
		    if (L < 0) L = -L;
		    dz = dz/L; dr = dr/L; dg = dg/L; db = db/L;
	            if (y1 < y2) { 
			for (int y=y1; y<=y2; y++) {
			    if (zBuffer[x1+width*y] < zCur) {
				theDevice->C3F(rCur, gCur, bCur);
				zBuffer[x1+width*y] = zCur;
				theDevice->V2F(x1,y);
			    }
			    zCur += dz;
			    rCur += dr;
			    gCur += dg;
			    bCur += db;
			}
		    } else {
			for (int y=y1; y>=y2; y--) {
			    if (zBuffer[x1+width*y] < zCur) {
				theDevice->C3F(rCur, gCur, bCur);
				zBuffer[x1+width*y] = zCur;
				theDevice->V2F(x1,y);
			    }
			    zCur += dz;
			    rCur += dr;
			    gCur += dg;
			    bCur += db;
			}
		    }
			
		} else if (dy == 0) {
		    int L = dx;
		    if (L < 0) L = -L;
		    dz = dz/L; dr = dr/L; dg = dg/L; db = db/L;
		    if (x1 < x2) {
			for (int x=x1; x<=x2; x++) {
			    if (zBuffer[x+width*y1] < zCur) {
				theDevice->C3F(rCur, gCur, bCur);
				zBuffer[x+width*y1] = zCur;
				theDevice->V2F(x,y1);
			    }
			    zCur += dz;
			    rCur += dr;
			    gCur += dg;
			    bCur += db;
			} 
		    } else {
			for (int x=x1; x>=x2; x--) {
			    if (zBuffer[x+width*y1] < zCur) {
				theDevice->C3F(rCur, gCur, bCur);
				zBuffer[x+width*y1] = zCur;
				theDevice->V2F(x,y1);
			    }
			    zCur += dz;
			    rCur += dr;
			    gCur += dg;
			    bCur += db;
			} 
		    } 			
		} else {
		    float yCur = point1->Y();
		    float xCur = point1->X(); 
		    float dx = (point2->X() -xCur);
		    float dy = (point2->Y() -yCur);
		    float length = sqrt(dx*dx + dy*dy);
		    dx /= length;
		    dy /= length;
		    int L = (int)length;
		    dz = dz/length; dr = dr/length; dg = dg/length; db = db/length;
		    for (int l=0; l<L; l++) {
			int x = (int)xCur;
			int y = (int)yCur;			    
			if (zBuffer[x+width*y] < zCur) {
			    theDevice->C3F(rCur, gCur, bCur);
			    zBuffer[x+width*y] = zCur;
			    theDevice->V2F(x,y);
			}
			xCur += dx;
			yCur += dy;			    
			zCur += dz;
			rCur += dr;
			gCur += dg;
			bCur += db;
		    }
		}
 	    }
	}
    }
#else
 MYPOINT *point1 = face.pointList.GetHead();
    if (point1 != 0) {
	if (fillMode == 1) {
	    DrawWire(face);
	} else {
	    MYPOINT *point1 = face.pointList.SetToHead();	
      	MYPOINT *point2 = face.pointList.Forward();	

	    if (point2 != 0) {
		int x1 = (int)point1->X();
		int y1 = (int)point1->Y();    
		float zCur = point1->Z();    	    
		float rCur = point1->r;
		float gCur = point1->g;
		float bCur = point1->b;

		int x2 = (int)point2->X();
		int y2 = (int)point2->Y();    	
		float z2 = point2->Z();    		    
		float r2 = point2->r;
		float g2 = point2->g;
		float b2 = point2->b;
				
 
  
  theDevice->BGNCLOSEDLINE();
  theDevice->C3F(rCur,gCur,bCur);
  theDevice->V2F(x1,y1);
  theDevice->C3F(r2,g2,b2);
  theDevice->V2F(x2,y2);
  theDevice->ENDCLOSEDLINE();
		}
	}
	}
  

#endif

    // now we delete the face
    FACE *theFace = &face;
    delete theFace;
}

void
ScanLineConverter::scanPolygon(FACE &face)
{
  MYPOINT *point = face.pointList.GetHead();
  if (point != 0) {
    if (fillMode == 1) {
      DrawWire(face);
    } else 
      DrawFill(face);
  }    
  
  // now we delete the face
  FACE *theFace = &face;
  delete theFace;
}

void
ScanLineConverter::DrawWire(FACE &face)
{
  MYPOINT *point;

  // set the color of wire to rgb of first point
  point = face.pointList.GetHead();
  theDevice->C3F(point->r, point->g, point->b);

  theDevice->BGNCLOSEDLINE();
  FOR_EACH(point, face.pointList)
    {
      theDevice->V2F(point->X(), point->Y());
    }
  theDevice->ENDCLOSEDLINE();
}
  

LIST<EDGE> FaceEdges(ORDER_DECREASING);
LIST<EDGE> ActiveEdges(ORDER_INCREASING);
  
void
ScanLineConverter::DrawFill(FACE &face)
{
//  MYPOINT *point;
  int y1, y2, cnt, slope;
  int Ymax, Ymin, Pixel, xleft, xright, dX;
  float Rcur, Gcur, Bcur, Zcur, dR_dX, dG_dX, dB_dX, dZ_dX;// Zfar;
//  float R,G,B,r1,g1,b1,r2,g2,b2,ka,kd,r,g,b,Iar,Iag,Iab,DotLN,N;
  MYPOINT *Point1, *Point2;
  EDGE *Edge, *Edge1, *Edge2;

  Ymax = 0;
  Ymin = 1000000; 

  FaceEdges.Clear();
  ActiveEdges.Clear();

  // for each pair of vertices create an edge and add to list of face edges
  FOR_EACH_PAIR(Point1,Point2,face.pointList,{
     y1 = (int)Point1->Y();
     y2 = (int)Point2->Y();
     if (y1 > y2) {
        if ((y1-1) >= y2) {  // only add if it will be drawn
	   Edge = new EDGE(Point1, Point2, 1);
           FaceEdges.AddSorted(*Edge, &Edge->Y_Top);
           if (y1 > Ymax) Ymax = y1;
           if (y2 < Ymin) Ymin = y2;
        }
     }
     else 
 	if ((y2-1) >= y1) {   // only add if it will be drawn
  	   Edge = new EDGE(Point2, Point1, -1);
     	   FaceEdges.AddSorted(*Edge, &Edge->Y_Top);
	   if (y2 > Ymax) Ymax = y2;
           if (y1 < Ymin) Ymin = y1;
     	  }
   });

   Ymax = Ymax - 1; // we start 1 scan line down          

   // now for each scanline determine which pixels should be off and 
   // which should be on
   for (int ScanLine = Ymax; ScanLine >= Ymin; ScanLine--) {

     // for each edge on active list increment the attributes & check 
     // if it should still be there
     FOR_EACH(Edge, ActiveEdges) {  
       if (Edge->edge_YBot() > ScanLine)  
	 Edge = ActiveEdges.RemoveCurrent();
       else
	 Edge->INCREMENT(); 
     }	

     // remove edges from faceedge list if it should be in active edges
     FOR_EACH(Edge, FaceEdges) {
       if (Edge->edge_YTop() > ScanLine) {
	 Edge = FaceEdges.RemoveCurrent();
	 
	 ActiveEdges.AddSorted(*Edge, &Edge->X_Cur);
       }
     }

     // resort the edges on the active list
     ActiveEdges.Resort();

     // if WINDING we must ensure if at bottom edges sorted also by slope
     if (scanMode == WINDING) {
       FOR_EACH_PAIR(Edge1, Edge2, ActiveEdges, {
	 if ((ScanLine == Edge1->edge_YBot())
	      && (ScanLine == Edge2->edge_YBot())
	      && (Edge1->edge_Xcur() == Edge2->edge_Xcur()))
	     
	     if (Edge1->edge_Slope() == 1) {
	       Edge1->edge_SlopeSet(-1);
	       Edge2->edge_SlopeSet(1);
	     }
       });
     }

     // now go across
     cnt = 0; slope = 0;
     FOR_EACH_PAIR(Edge1, Edge2, ActiveEdges, {
       xleft = Edge1->edge_Xcur();
       xright = Edge2->edge_Xcur();
       slope = slope + Edge1->edge_Slope() + Edge2->edge_Slope();

       if (xleft < xright) {

	 cnt++;
	 if (((cnt%2) == 1) || (scanMode == WINDING && slope != 0)) {

	   Zcur = Edge1->edge_Zcur();
	   Rcur = Edge1->edge_Rcur();
	   Gcur = Edge1->edge_Gcur();
	   Bcur = Edge1->edge_Bcur();
	   dX = xleft - xright;
	   dZ_dX = (Zcur - Edge2->edge_Zcur()) / dX;
	   dR_dX = (Rcur - Edge2->edge_Rcur()) / dX;
	   dB_dX = (Bcur - Edge2->edge_Bcur()) / dX;
	   dG_dX = (Gcur - Edge2->edge_Gcur()) / dX;

	   theDevice->BGNPOINT();
	   for (Pixel = xleft; Pixel < xright; Pixel++) {
	     if (zBuffer[Pixel+width*ScanLine] < Zcur) {
	       theDevice->C3F(Rcur, Gcur, Bcur);
	       zBuffer[Pixel+width*ScanLine] = Zcur;
	       theDevice->V2F(Pixel,ScanLine);
	     }
	     Zcur += dZ_dX;
	     Rcur += dR_dX;
	     Gcur += dG_dX;
	     Bcur += dB_dX;
	   }
	   theDevice->ENDPOINT();
	 }
       }	   

       else if (xleft == xright) {  // crossing edges or a vertex

	 if ((ScanLine == Edge1->edge_YTop() 
	      && ScanLine == Edge2->edge_YTop()) || 
	     (ScanLine == Edge1->edge_YBot() 
	      && ScanLine == Edge2->edge_YBot())  )

	   // meeting vertices top or bottom -- for even odd count as 2
	   cnt += 1;

	 else if ((ScanLine == Edge1->edge_YTop() 
		   && ScanLine == Edge2->edge_YBot()) || 
		  (ScanLine == Edge1->edge_YBot() 
		   && ScanLine == Edge2->edge_YTop())   )

	   // meeting vertices, one top one bottom -- for even odd count as 1
	   cnt += 2;

	 else   // crossing edges
	   cnt += 1;
       }
     });
   }
}

  
void
ScanLineConverter::setDevice(Device &device)
{
    theDevice = &device;
}

  
void
ScanLineConverter::setProjection(Projection &projection)
{
    theProjection = &projection;
}
