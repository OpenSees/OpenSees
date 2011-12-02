#ifndef GeometricBrickDecorator_h
#define GeometricBrickDecorator_h

#include <stdlib.h>
#include <iostream>
#include <OPS_Globals.h>
#include <StandardStream.h>
#include <math.h>


#include <ArrayOfTaggedObjects.h>

// includes for the domain classes
#include <Domain.h>
#include <Node.h>
#include <Element.h>
#include <ZeroLength.h>

#include <ElasticMaterial.h>
#include <ElasticIsotropicThreeDimensional.h>

#include <SP_Constraint.h>

#include <SingleDomAllSP_Iter.h>
#include <NodeIter.h>

#include <Vertex.h>
#include <VertexIter.h>
#include <Graph.h>

#include <Vector.h>
#include <Matrix.h>

#include <map>
#include <set>

class GeometricBrickDecorator {



  public :
    
GeometricBrickDecorator();

virtual ~GeometricBrickDecorator();
    
void setBrick(Element* elem);
    
void setDomain(Domain * theDomain);
    
void clearBrick();
    
void clearDomain();
  
void clearAll();
    
bool isClear();
    
bool isEmpty();

bool isZero(double a, double b);
    
bool isPointInVolume(const Vector &pP);
    
double evalSignedVol(const Vector &pA, 
		     const Vector &pB,
		     const Vector &pC,
		     const Vector &pP);
    
double eval3OrderDet(const Vector &pA, 
		     const Vector &pB,
		     const Vector &pC);
    
void getFace(int which, ID& face, ID& faceID);
    
double getMinMaxCrds(int which, int whichtoret);
    
bool compareFaceToFace(int which, ID &faceOther);
  

bool isFaceinPlane(int which, const Vector &pP);
  

bool isFaceinVertPlane(int which, double xy, double zmin, double zmax, 
       		       int whichCrd);
    
bool isLeftBoundary(double xMin, double xMax, double yMin, double yMax, double zMin, double zMax);

bool isRightBoundary(double xMin, double xMax, double yMin, double yMax, double zMin, double zMax);
      
bool isFrontBoundary(double xMin, double xMax, double yMin, double yMax, double zMin, double zMax);
      
bool isRearBoundary(double xMin, double xMax, double yMin, double yMax, double zMin, double zMax);

bool isBottomBoundary(double xMin, double xMax, double yMin, double yMax, double zMin, double zMax);

bool isBoundaryLayerEle(double xMin, double xMax, double yMin, double yMax, double zMin, double zMax);

    private :
      
Domain * myDomain;
    
Element * myBrick;

};
    
#endif



