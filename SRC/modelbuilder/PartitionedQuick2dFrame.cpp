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
// $Date: 2003-02-14 23:01:47 $
// $Source: /usr/local/cvs/OpenSees/SRC/modelbuilder/PartitionedQuick2dFrame.cpp,v $
                                                                        
                                                                        

#include <PartitionedQuick2dFrame.h>
#include <PartitionedDomain.h>
#include <Subdomain.h>
#include <iOPS_Stream.h>
#include <fstream.h>
#include <stdlib.h>

#include <beam2d02.h>
#include <beam2d03.h>
#include <beam2d04.h>

#include <NodalLoad.h>
#include <Node.h>
#include <SP_Constraint.h>
#include <MP_Constraint.h>

//  PartitionedQuick2dFrameModel();
//	constructor
PartitionedQuick2dFrame::PartitionedQuick2dFrame(PartitionedDomain &theDomain,
						 int numXX, int numYY, 
						 int eleTType)
:PartitionedModelBuilder(theDomain,
	      PartitionedModelBuilder_TAGS_PartitionedQuick2dFrameModel),
 numX(numXX), numY(numYY), eleType(eleTType),
 numEle(0),numNode(0),numSPs(0),numMPs(0),numLCs(0)
{
    // does nothing
}


PartitionedQuick2dFrame::PartitionedQuick2dFrame(Subdomain &theDomain)
:PartitionedModelBuilder(theDomain,
	     PartitionedModelBuilder_TAGS_PartitionedQuick2dFrameModel),
 numX(0), numY(0), eleType(0),
 numEle(0),numNode(0),numSPs(0),numMPs(0),numLCs(0)
{
    // does nothing
}



// ~PartitionedQuick2dFrame();    
//	destructor,

PartitionedQuick2dFrame::~PartitionedQuick2dFrame()
{
    // does nothing
}    

int
PartitionedQuick2dFrame::buildInterface(int numSubdomains) 
{
  
    Domain *theDomain = this->getDomainPtr();    
    if (theDomain == 0) {
	opserr << "FATAL:PlaneFrame::buildModel(void) -";
	opserr << " no associated domain.\n";
	exit(-1);
	return -1;
    }	

    if (numSubdomains == 0) {
	opserr << "FATAL:PlaneFrame::buildModel(void) -";
	opserr << " numSubdomains = 0\n";
	exit(-1);
	return -2;
    }	
    
    /*
     *  Read in the Nodal Data
     */
    int i,j;
    double  LL_x, LL_y, UR_x, UR_y, DIV_x, DIV_y, CURR_y, CURR_x;
    int    NO_x, NO_y;
    LL_x = 0;
    LL_y = 0;
    UR_x = 1;
    UR_y = 2;

    NO_x = numX + 1;
    NO_y = numY + 1;
    numNode = NO_x * NO_y;
    numSPs = NO_x;
    int numLoads = NO_y-1;
    numLCs = 1;
    
    DIV_x = (UR_x - LL_x) / (NO_x - 1);
    DIV_y = (UR_y - LL_y) / (NO_y - 1);
    CURR_y = LL_y;

    int numSubstructureStories = numY / numSubdomains;

    DIV_y *= numSubstructureStories;
    LL_y += DIV_y;
    int ok;
    int nodeNum = 1 + NO_x*numSubstructureStories;

    for (i=1; i<numSubdomains; i++) {
	CURR_x = LL_x;

	for (j=0; j<NO_x; j++) {

	    Node *nodePtr = new Node(nodeNum,3,CURR_x,CURR_y);
	    if (nodePtr != 0) {
		ok = theDomain->addNode(nodePtr);
		if (ok < 0) {
		    opserr << "PartitionedQuick2dFrame::PartitionedQuick2dFrame::Input ";
		    opserr << " could not add Node " << nodeNum << endln;
		    return -1;
		}
		nodeNum++;
		CURR_x += DIV_x;
	    } else {
		opserr << "PartitionedQuick2dFrame::PartitionedQuick2dFrame::Input ";
		opserr << " ran out of memory creating Node " << nodeNum << endln;
		return -1;
	    }

	}
	CURR_y += DIV_y;
	nodeNum += (numSubstructureStories-1)*NO_x;
    }


    double x_ld = 10;
    double y_ld = 0;
    double theta_ld = 0;


    //
    // create the nodal loads
    //
    numLoads = numSubdomains-1;

    Vector forces(3);
    forces(0) = x_ld; forces(1) = y_ld; forces(2) = theta_ld;
    for (i=0; i<numLoads; i++) {

	nodeNum = (i+1)*NO_x*numSubstructureStories+1;
	NodalLoad *nodeLoadPtr = new NodalLoad(nodeNum,forces);
	if (nodeLoadPtr != 0) {
	    ok =  theDomain->addNodalLoad(nodeLoadPtr);
	    if (ok < 0) {
		opserr << "PartitionedQuick2dFrame::PartitionedQuick2dFrame::Input ";
		opserr << " could not add nodal load " << nodeNum << endln;
		return -1;
	    }
	} else {
	    opserr << "PartitionedQuick2dFrame::PartitionedQuick2dFrame::Input ran out of ";
	    opserr << "memory creating Nodal Load " << nodeNum << endln;
	    return -1;
	}	    	
    }


    return 0;
}





int
PartitionedQuick2dFrame::buildSubdomain(int subdomain, 
					int numSubdomains,
					Subdomain &theSubdomain) 
{
    subdomain--; // assumes subdomains labelled 0 through NSub-1 in rest of this
    
    int numSect, numMat, numLoads;
    int i,j;

    int    x_rst,y_rst,theta_rst;
    float  x_ld, y_ld, theta_ld;    
    
    numSect = 2;
    numMat = 2;

    /*
     * Allocate Space for all the variables
     */
    int SectDat = 2; // A, I
    int MatDat = 1; // E

    Vector Sections(numSect*SectDat);
    Vector Materials(numMat*MatDat);

    /*
     *  Set the Sections Data & Material Data
     */
    

    for (i=0; i<numSect; i++) {
	Sections(i*SectDat) = 1; Sections(i*SectDat+1) = 1;
    }

    for (i=0; i<numMat; i++) {
	Materials(i*MatDat) = 1;
    }

    /*
     *  Read in the Nodal Data
     */
    double  LL_x, LL_y, UR_x, UR_y, DIV_x, DIV_y, CURR_y, CURR_x;
    int    NO_x, NO_y;
    LL_x = 0;
    LL_y = 0;
    UR_x = 1;
    UR_y = 2;

    NO_x = numX + 1;
    NO_y = numY + 1;
    numNode = NO_x * NO_y;
    numSPs = NO_x;
    numLoads = NO_y-1;
    numLCs = 1;
    
    DIV_x = (UR_x - LL_x) / (NO_x - 1);
    DIV_y = (UR_y - LL_y) / (NO_y - 1);
    CURR_y = LL_y;

    int numSubdomainStories = numY / numSubdomains;
    int nodeNum = 1 + subdomain*numSubdomainStories*NO_x;
    LL_y += DIV_y*numSubdomainStories*subdomain;
    if (subdomain == (numSubdomains-1))
      numSubdomainStories += numY - numSubdomainStories*numSubdomains;

    int ok;
    for (i=0; i<numSubdomainStories+1; i++) {
	CURR_x = LL_x;
	for (j=0; j<NO_x; j++) {
	    Node *nodePtr = new Node(nodeNum,3,CURR_x,CURR_y);
	    if (nodePtr != 0) {
	      if (i == 0 && subdomain != 0) 
		ok = theSubdomain.addExternalNode(nodePtr);
	      else if (i == numSubdomainStories && subdomain != (numSubdomains-1))
		ok = theSubdomain.addExternalNode(nodePtr);
	      else 
		ok = theSubdomain.addNode(nodePtr);

	      if (ok < 0) {
		opserr << "PartitionedQuick2dFrame::PartitionedQuick2dFrame::Input ";
		opserr << " could not add Node " << nodeNum << endln;
		return -1;
	      }
	      nodeNum++;
	      CURR_x += DIV_x;
	    } else {
		opserr << "PartitionedQuick2dFrame::PartitionedQuick2dFrame::Input ";
		opserr << " ran out of memory creating Node " << nodeNum << endln;
		return -1;
	    }
	}
	CURR_y += DIV_y;
    }

    /*
     *  do the beam stuff
     */


    // set the starting node num and element numbers 
    numSubdomainStories = numY / numSubdomains;
    int start_node = 1 + subdomain*numSubdomainStories*NO_x;
    int eleNum = 0 + NO_x * numSubdomainStories*subdomain + (NO_x-1)*numSubdomainStories*subdomain;
    if (subdomain == (numSubdomains-1))
      numSubdomainStories += numY - numSubdomainStories*numSubdomains;


    for (i=0; i<numSubdomainStories; i++) {
	int inode = start_node;
        int jnode = start_node + NO_x;
	Element *elePtr =0;	

	/* add the columns */
	double A = Sections(0);
	double E = Materials(0);
	double I = Sections(1);
	
        for (j=0; j < NO_x; j++) {
	    // create the element
	    if (eleType == 2)
		elePtr = new beam2d02(eleNum,A,E,I,inode,jnode);
	    else if (eleType == 3)
		elePtr = new beam2d03(eleNum,A,E,I,inode,jnode);
	    else
		elePtr = new beam2d04(eleNum,A,E,I,inode,jnode);
	    
	    // add it to the domain
	    if (elePtr != 0) {
		ok = theSubdomain.addElement(elePtr);

		if (ok < 0) {
		    opserr << "PartitionedQuick2dFrame::PartitionedQuick2dFrame::Input ";
		    opserr << " could not add Element " << eleNum << endln;
		    return -1;
		}
	    } else {
		opserr << "PartitionedQuick2dFrame::PartitionedQuick2dFrame::Input ";
		opserr << " ran out of memory creating Element " << eleNum << endln;
		return -1;
	    }	    
	    inode++;
	    jnode++;
	    eleNum++;
	}
	/*add the beams */
	A = Sections(2);
	E = Materials(1);
	I = Sections(3);
	jnode = inode++;	
        for (j=0; j< (NO_x-1); j++) {
	    // create the element
	    if (eleType == 2)
		elePtr = new beam2d02(eleNum,A,E,I,jnode,inode);
	    else if (eleType == 3)
		elePtr = new beam2d03(eleNum,A,E,I,jnode,inode);
	    else
		elePtr = new beam2d04(eleNum,A,E,I,jnode,inode);

	    // add it to the domain
	    if (elePtr != 0) {
		ok = theSubdomain.addElement(elePtr);
		if (ok < 0) {
		    opserr << "PartitionedQuick2dFrame::PartitionedQuick2dFrame::Input ";
		    opserr << " could not add Element " << eleNum << endln;
		    return -1;
		}
	    } else {
		opserr << "PartitionedQuick2dFrame::PartitionedQuick2dFrame::Input ";
		opserr << " ran out of memory creating Element " << eleNum << endln;
		return -1;
	    }	    
	    inode++;
	    jnode++;
	    eleNum++;
	}
        start_node += NO_x;
    }

    /*
     *  Read in the Restraint Data
     */

    x_rst = 1;
    y_rst = 1;
    theta_rst = 1;

    if (subdomain == 0) {
      nodeNum = 1;
      for (i=0; i<NO_x; i++) {
	SP_Constraint *spPtr;
	
	if (x_rst == 1) {
	  spPtr = new SP_Constraint(nodeNum,0,0.0);
	  // add it to the domain
	  if (spPtr != 0) {
		ok = theSubdomain.addSP_Constraint(spPtr);
		if (ok < 0) {
		    opserr << "PartitionedQuick2dFrame::PartitionedQuick2dFrame::Input ";
		    opserr << " could not add SP_Constraint " << nodeNum << endln;
		    return -1;
		}
	    } else {
		opserr << "PartitionedQuick2dFrame::PartitionedQuick2dFrame::Input ran out of ";
		opserr << "memory creating SP_Constraint " << nodeNum << endln;
		return -1;
	    }	    
	}
	
	if (y_rst == 1) {
	    spPtr = new SP_Constraint(nodeNum,1,0.0);
	    // add it to the domain
	    if (spPtr != 0) {
		ok = theSubdomain.addSP_Constraint(spPtr);
		if (ok < 0) {
		    opserr << "PartitionedQuick2dFrame::PartitionedQuick2dFrame::Input ";
		    opserr << " could not add SP_Constraint " << nodeNum << endln;
		    return -1;
		}
	    } else {
		opserr << "PartitionedQuick2dFrame::PartitionedQuick2dFrame::Input ran out of ";
		opserr << "memory creating SP_Constraint " << nodeNum << endln;
		return -1;
	    }	    
	}	

	if (theta_rst == 1) {
	    spPtr = new SP_Constraint(nodeNum,2,0.0);
	    // add it to the domain
	    if (spPtr != 0) {
		ok = theSubdomain.addSP_Constraint(spPtr);
		if (ok < 0) {
		    opserr << "PartitionedQuick2dFrame::PartitionedQuick2dFrame::Input ";
		    opserr << " could not add SP_Constraint " << nodeNum << endln;
		    return -1;
		}
	    } else {
		opserr << "PartitionedQuick2dFrame::PartitionedQuick2dFrame::Input ran out of ";
		opserr << "memory creating SP_Constraint " << nodeNum << endln;
		return -1;
	    }	    
	}	
	nodeNum++;
      }
    }


    numSubdomainStories = numY / numSubdomains;
    start_node = 1 + subdomain*numSubdomainStories*NO_x;
    if (subdomain == (numSubdomains-1))
      numSubdomainStories += numY - numSubdomainStories*numSubdomains;
    
    x_ld = 10;
    y_ld = 0;
    theta_ld = 0;

    //
    // create the nodal loads
    //

    Vector forces(3);
    forces(0) = x_ld; forces(1) = y_ld; forces(2) = theta_ld;
    for (i=1; i<numSubdomainStories; i++) { // don't add loads at interface nodes

	nodeNum = start_node+i*NO_x;
	NodalLoad *nodeLoadPtr = new NodalLoad(nodeNum,forces);
	if (nodeLoadPtr != 0) {
	    ok =  theSubdomain.addNodalLoad(nodeLoadPtr);
	    if (ok < 0) {
		opserr << "PartitionedQuick2dFrame::PartitionedQuick2dFrame::Input ";
		opserr << " could not add nodal load " << nodeNum << endln;
		return -1;
	    }
	} else {
	    opserr << "PartitionedQuick2dFrame::PartitionedQuick2dFrame::Input ran out of ";
	    opserr << "memory creating Nodal Load " << nodeNum << endln;
	    return -1;
	}	    	
    }

    if (subdomain == (numSubdomains-1)) {
	nodeNum = start_node+numSubdomainStories*NO_x;
	NodalLoad *nodeLoadPtr = new NodalLoad(nodeNum,forces);
	if (nodeLoadPtr != 0) {
	    ok = 	theSubdomain.addNodalLoad(nodeLoadPtr);
	    if (ok < 0) {
		opserr << "PartitionedQuick2dFrame::PartitionedQuick2dFrame::Input ";
		opserr << " could not add nodal load " << nodeNum << endln;
		return -1;
	    }
	} else {
	    opserr << "PartitionedQuick2dFrame::PartitionedQuick2dFrame::Input ran out of ";
	    opserr << "memory creating Nodal Load " << nodeNum << endln;
	    return -1;
	}	    	
    }

    return 0;
}



int 
PartitionedQuick2dFrame::sendSelf(int cTag, Channel &theChannel)
{
  int dataTag = this->getDbTag();
  ID data(3);
  data(0) = numX;
  data(1) = numY;
  data(2) = eleType;

  theChannel.sendID(dataTag, cTag, data);
  return 0;
}

int 
PartitionedQuick2dFrame::recvSelf(int cTag, Channel &theChannel, 
				  FEM_ObjectBroker &theBroker)
{
  int dataTag = this->getDbTag();
  ID data(3);
  theChannel.recvID(dataTag, cTag, data);
  numX = data(0);
  numY = data(1);
  eleType = data(2);
  
  return 0;

}
