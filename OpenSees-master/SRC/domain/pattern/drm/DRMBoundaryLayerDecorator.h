#ifndef DRMBoundaryLayerDecorator_h
#define DRMBoundaryLayerDecorator_h

#include <Domain.h>
#include <Element.h>
#include <Node.h>
#include <Vector.h>
#include <Matrix.h>
#include <map>
#include <set>


using namespace std;

//class Element;

class DRMBoundaryLayerDecorator
{

  
  //friend class Element;  
  public :

  DRMBoundaryLayerDecorator();

  virtual ~DRMBoundaryLayerDecorator();

  void setArray(int *eNodes);

  void setID(ID* eNodes);

  void setMap(std::map<int,int> &eNodes);

  void setSet(std::set<int,std::less<int > > &eNodes);

  void setDomain(Domain *theDomain);
  
  void clearDomain();

  void setBrick(Element* brickTag);

  void clearBrick();

  void clear();
  
  void computeDRMLoad(Vector &drmLoad, 
		      const Vector &displ,
		      const Vector &veloc,
		      const Vector &accel);

  void applyDRMLoad(double cfact,
		    Vector &drmLoad,
		    const Vector &displ,
		    const Vector &veloc,
		    const Vector &accel);	     
  
  
  
  private:
  
  int cons;
  Domain * myDomain;
  Element *myBrick;
  std::map<int,int> eNodeMap;
  std::set<int,std::less<int > > eNodeSet;
  int *eNodeArray;
  ID* eNodeID;

  Vector* Peff_m;
  Vector* Peff_k;
  Vector* Peff_d;
  
  void get_E_B_Nodes(ID &e, ID &b);

  void zeroSubmatrix(Matrix &M, int e, int b);  
};


#endif


