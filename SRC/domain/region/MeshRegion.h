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
// $Date: 2003-02-14 23:01:02 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/region/MeshRegion.h,v $
                                                                        
                                                                        
// Written: fmk 
//
// Description: This file contains the class definition for MeshRegion.
// A Region is a part of the domain which is defined by a set of
// Elements and Nodes (all the end nodes of the elements are in the region, 
// as are all elements whose end nodes are in the region)
//
// What: "@(#) Region.h, revA"

#ifndef MeshRegion_h
#define MeshRegion_h

#include <DomainComponent.h>
#include <ID.h>

class Element;
class Node;
class ElementRecorder;
class NodeRecorder;

class MeshRegion : public DomainComponent
{
  public:
    MeshRegion(int tag);
    MeshRegion(int tag, int classTag);    
    virtual ~MeshRegion();

    // methods dealing with setting up the region
    virtual int setNodes(const ID &theNodes);
    virtual int setElements(const ID &theEles);
    virtual void setExtraEles(const ID& theEles) {xEles=theEles;}

    // methods getting the ID's of nodes & ele
    virtual const ID &getNodes(void);
    virtual const ID &getElements(void);
    virtual const ID &getExtraEles() const {return xEles;}

    // methods dealing with setting parameters in the region
    virtual int setRayleighDampingFactors(double alphaM, 
					  double betaK, 
					  double betaK0,
					  double betaKc);

    // methods to send & recv data for database/parallel applications
    virtual int sendSelf(int commitTag, Channel &theChannel);
    virtual int recvSelf(int commitTag, Channel &theChannel, 
			 FEM_ObjectBroker &theBroker);
    virtual void Print(OPS_Stream &s, int flag =0);       

  protected:
    
  private:
    double alphaM, betaK, betaK0, betaKc;

    ID *theNodes;
    ID *theElements;
    ID xEles;

    int	   currentGeoTag;
    int    lastGeoSendTag;
    int dbNod;
    int dbEle;
};


#endif

