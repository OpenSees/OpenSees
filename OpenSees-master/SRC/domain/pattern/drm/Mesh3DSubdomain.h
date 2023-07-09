#ifndef Mesh3DSubdomain_h
#define Mesh3DSubdomain_h

#include <Element.h>
#include <map>
#include <set>

class Domain;


class Mesh3DSubdomain {

  public :
    
    //general constructor
    Mesh3DSubdomain(Domain * inpDomain) ;


    Mesh3DSubdomain() ;


    //destructor
    virtual ~Mesh3DSubdomain();
    
    int lastEle();

    int lastNode();
    
    int lastStartNode();
  
    void allocate_e_Nodes(double xMin, double xMax, 
			  double yMin, double yMax,
			  double zMin, double zMax,
			  std::map<int,int>& eNodes);

    void allocateBoundaryLayerElements(double xMin, double xMax,
				       double yMin, double yMax,
				       double zMin, double zMax,
				       std::map<int,Element*>& elements,
				       std::map<int,Vector*>& storage,
				       std::map<int,int>& storage2);

    private :
     
    Domain * myDomain;
    int myLastStartNodeTag;
    int myLastEle;
    int myLastNode;
};

#endif






