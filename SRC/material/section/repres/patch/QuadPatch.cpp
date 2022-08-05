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
// $Date: 2008-04-14 21:41:12 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/repres/patch/QuadPatch.cpp,v $
                                                                        
                                                                        
// File: QuadPatch.C
// Written by Remo M. de Souza
// December 1998

#include <Matrix.h>
#include <Patch.h>
#include <QuadPatch.h>
#include <QuadCell.h>
#include <elementAPI.h>

void* OPS_QuadPatch()
{
    if(OPS_GetNumRemainingInputArgs() < 11) {
	opserr<<"insufficient arguments for QuadPatch\n";
	return 0;
    }

    // get idata
    int numData = 3;
    int idata[3];
    if(OPS_GetIntInput(&numData,&idata[0]) < 0) return 0;

    // get data
    static Matrix vertexCoords(4,2);
    double data[8];
    numData = 8;
    if(OPS_GetDoubleInput(&numData,&data[0]) < 0) return 0;
    for(int i=0; i<4; i++) {
	for(int j=0; j<2; j++) {
	    vertexCoords(i,j) = data[i*2+j];
	}
    }
<<<<<<< HEAD
	
	
	bool computePositive = false;
	while (OPS_GetNumRemainingInputArgs() > 0) {
      const char* poa = OPS_GetString();
      if (strcmp(poa,"-positiveArea") == 0) {
	      computePositive = true;
      }
    }
	
	// Check if counter clockwise or not
    if (computePositive==true) { 
		double v1z, v1y, v2z, v2y, theta1, theta2, centroidz, centroidy;
		double z0, y0, z1, y1, result;
		  
		y0 = vertexCoords(0,0) ;
		z0 = vertexCoords(0,1) ;

		y1 = vertexCoords(1,0) ;
		z1 = vertexCoords(1,1) ;  

		  
		centroidx=0;
		centroidy=0;
		for(int j=0; j<4; j++) {
				centroidy=centroidy+vertexCoords(j,0);
				centroidz=centroidz+vertexCoords(j,1);
			}

		centroidy=centroidy/4;
		centroidz=centroidz/4;  
		  
		  
		v1y=y0-centroidy;
		v1z=z0-centroidz;
		  
		v2y=y1-centroidy;
		v2z=z1-centroidz;
		  
		theta1 = atan2 (v1y,v1z) * 180 / PI;
		theta2 = atan2 (v2y,v2z) * 180 / PI;

		// error of points definition 
		if (theta2==theta1 or theta2==180+theta1 or theta2==-180+theta1)  {
			return -1;
		}
		
		
		// clock-wise definition by user and needs to modify
		// Rotate coordinates so we get all positive areas

		if (theta2<theta1)  {
		  int vertexCoordtmp [4][2];
		
		  for(int i=0; i<4; i++) {
			for(int j=0; j<2; j++) {
			  vertexCoordtmp[i][j]=vertexCoords(i,j);
			}
		  }
		
		  for(int i=0; i<4; i++) {
			for(int j=0; j<2; j++) {
			  vertexCoords(3-i,j)=vertexCoordtmp[i][j];
			}
		  }
		
		}
	}
	
	
=======

>>>>>>> parent of 9d6aaf17 (Update QuadPatch.cpp)
    return new QuadPatch(idata[0],idata[1],idata[2],vertexCoords);
}

void* OPS_RectPatch()
{
    if(OPS_GetNumRemainingInputArgs() < 7) {
	opserr<<"insufficient arguments for RectPatch\n";
	return 0;
    }

    // get idata
    int numData = 3;
    int idata[3];
    if(OPS_GetIntInput(&numData,&idata[0]) < 0) return 0;

    // get data
    static Matrix vertexCoords(4,2);
    double data[4];
    numData = 4;
    if(OPS_GetDoubleInput(&numData,&data[0]) < 0) return 0;
    double dyOverdz = (data[2]-data[0])/(data[3]-data[1]);
    if (dyOverdz < 0) {
      // Swap coordinates so we get all positive areas
      double tmp = data[1];
      data[1] = data[3];
      data[3] = tmp;
    }
    vertexCoords(0,0) = data[0];
    vertexCoords(0,1) = data[1];
    vertexCoords(1,0) = data[2];
    vertexCoords(1,1) = data[1];
    vertexCoords(2,0) = data[2];
    vertexCoords(2,1) = data[3];
    vertexCoords(3,0) = data[0];
    vertexCoords(3,1) = data[3];

    return new QuadPatch(idata[0],idata[1],idata[2],vertexCoords);
}


QuadPatch::QuadPatch():
                       matID(0), nDivIJ(1), nDivJK(1), vertCoord(4,2)
{

}


QuadPatch::QuadPatch(int materialID, int numSubdivIJ, int numSubdivJK,
                     const Matrix &vertexCoords):
                       matID(materialID), 
                       nDivIJ(numSubdivIJ), nDivJK(numSubdivJK),
                       vertCoord(vertexCoords)

{

}  


QuadPatch::~QuadPatch()
{

}

void QuadPatch::setMaterialID(int materialID)
{
   matID = materialID;
}


void QuadPatch::setDiscretization(int numSubdivIJ, int numSubdivJK)
{
   nDivIJ = numSubdivIJ;
   nDivJK = numSubdivJK;
}

void QuadPatch::setVertCoords(const Matrix &vertexCoords)
{
   vertCoord = vertexCoords;
}

int QuadPatch::getMaterialID(void) const
{
   return matID;
}
 
void QuadPatch::getDiscretization(int &numSubdivIJ, int &numSubdivJK) const
{
   numSubdivIJ = nDivIJ;
   numSubdivJK = nDivJK;
}

const Matrix & QuadPatch::getVertCoords (void) const
{
   return vertCoord;
}

int QuadPatch::getNumCells (void) const
{
   return nDivIJ * nDivJK;
}

Cell **
QuadPatch::getCells (void) const
{
   double deltaXi;
   double deltaEta; 
   Matrix cellVertCoord(4,2);
   Vector N(4);
   double xi, eta;
   int i, j, k, r, s;
   int numCells;
   Cell **cells;
  
   if (nDivIJ > 0  && nDivJK > 0)
   {
      numCells  = this->getNumCells();

      cells = new Cell*  [numCells];
      
      if (!cells)
         return 0;

      deltaXi  = 2.0 / nDivIJ;
      deltaEta = 2.0 / nDivJK;

      k = 0;
      for (j = 0; j < nDivJK; j++)
         for (i = 0; i < nDivIJ; i++)
         {
            // compute natural coordinates

            cellVertCoord(0,0) = -1.0 + deltaXi  * i;
            cellVertCoord(0,1) = -1.0 + deltaEta * j;
            cellVertCoord(1,0) = -1.0 + deltaXi  * (i+1);
            cellVertCoord(1,1) = cellVertCoord(0,1);
            cellVertCoord(2,0) = cellVertCoord(1,0);
            cellVertCoord(2,1) = -1.0 + deltaEta * (j+1);
            cellVertCoord(3,0) = cellVertCoord(0,0);
            cellVertCoord(3,1) = cellVertCoord(2,1);

            // map to cartesian coordinates using bilinear
            // shape functions

            for (r = 0; r < 4; r++)
            {
               xi  = cellVertCoord(r,0);
               eta = cellVertCoord(r,1);
 
               N(0) = (1.0 - xi)*(1.0 - eta)/4.0;
               N(1) = (1.0 + xi)*(1.0 - eta)/4.0;
               N(2) = (1.0 + xi)*(1.0 + eta)/4.0;
               N(3) = (1.0 - xi)*(1.0 + eta)/4.0;

               cellVertCoord(r,0) = 0.0;
               cellVertCoord(r,1) = 0.0;

               for (s = 0; s < 4; s++)
               {
                  cellVertCoord(r,0) += N(s) * vertCoord(s,0);
                  cellVertCoord(r,1) += N(s) * vertCoord(s,1);
               }
            }  

            cells[k] = new QuadCell(cellVertCoord); 
            //opserr << "\ncreating cells Cell " << k << " :" << cells[k];
            k++; 
         }
   }
   else
      return 0;

   return cells;
}


Patch * 
QuadPatch::getCopy (void) const
{
   QuadPatch *theCopy = new QuadPatch (matID, nDivIJ, nDivJK, vertCoord);
   return theCopy;
}
 
void QuadPatch::Print(OPS_Stream &s, int flag) const
{
   s << "\nPatch Type: QuadPatch";
   s << "\nMaterial Id: " << matID;
   s << "\nNumber of subdivisions in the IJ direction: " << nDivIJ;
   s << "\nNumber of subdivisions in the JK direction: " << nDivJK;
   s << "\nVertex Coordinates: " << vertCoord;
}


OPS_Stream &operator<<(OPS_Stream &s, QuadPatch &quadPatch)
{
   quadPatch.Print(s);
   return s;
}
