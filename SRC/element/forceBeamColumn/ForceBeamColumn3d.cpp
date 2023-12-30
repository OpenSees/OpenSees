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
                                                                        
// $Revision: 1.33 $
// $Date: 2010-09-13 21:26:10 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/forceBeamColumn/ForceBeamColumn3d.cpp,v $

/*
 * References
 *

State Determination Algorithm
---
Neuenhofer, A. and F. C. Filippou (1997). "Evaluation of Nonlinear Frame Finite
Element Models." Journal of Structural Engineering, 123(7):958-966.

Spacone, E., V. Ciampi, and F. C. Filippou (1996). "Mixed Formulation of
Nonlinear Beam Finite Element." Computers and Structures, 58(1):71-83.


Plastic Hinge Integration
---
Scott, M. H. and G. L. Fenves (2006). "Plastic Hinge Integration Methods for
Force-Based Beam-Column Elements." Journal of Structural Engineering,
132(2):244-252.


Analytical Response Sensitivity (DDM)
---
Scott, M. H., P. Franchin, G. L. Fenves, and F. C. Filippou (2004).
"Response Sensitivity for Nonlinear Beam-Column Elements."
Journal of Structural Engineering, 130(9):1281-1288.


Software Design
---
Scott, M. H., G. L. Fenves, F. T. McKenna, and F. C. Filippou (2007).
"Software Patterns for Nonlinear Beam-Column Models."
Journal of Structural Engineering, Approved for publication, February 2007.

 *
 */

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

#include <Information.h>
#include <Parameter.h>
#include <ForceBeamColumn3d.h>
#include <MatrixUtil.h>
#include <Domain.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Renderer.h>
#include <math.h>
#include <elementAPI.h>
#include <string>
#include <fstream>
using std::ifstream;
#include <fstream>
#include <ElementResponse.h>
#include <CompositeResponse.h>
#include <ElementalLoad.h>
#include <ElementIter.h>

#define DefaultLoverGJ 1.0e-10

Matrix ForceBeamColumn3d::theMatrix(12,12);
Vector ForceBeamColumn3d::theVector(12);
double ForceBeamColumn3d::workArea[200];

Vector ForceBeamColumn3d::vsSubdivide[maxNumSections];
Matrix ForceBeamColumn3d::fsSubdivide[maxNumSections];
Vector ForceBeamColumn3d::SsrSubdivide[maxNumSections];

void* OPS_ForceBeamColumn3d()
{
    int dampingTag = 0;
    Damping* theDamping = 0;
    if (OPS_GetNumRemainingInputArgs() < 5) {
	opserr<<"insufficient arguments:eleTag,iNode,jNode,transfTag,integrationTag\n";
	return 0;
    }

    int ndm = OPS_GetNDM();
    int ndf = OPS_GetNDF();
    if(ndm != 3 || ndf != 6) {
	opserr<<"ndm must be 3 and ndf must be 6\n";
	return 0;
    }

    // inputs: 
    int iData[5];
    int numData = 5;
    if(OPS_GetIntInput(&numData,&iData[0]) < 0) {
	opserr << "WARNING invalid int inputs\n";
	return 0;
    }

    // options
    double mass = 0.0, tol=1e-12, subFac=10.0;
    int maxIter = 10, numSub = 4;
    numData = 1;
    while(OPS_GetNumRemainingInputArgs() > 0) {
	const char* type = OPS_GetString();
	if(strcmp(type,"-iter") == 0) {
	    if(OPS_GetNumRemainingInputArgs() > 1) {
		if(OPS_GetIntInput(&numData,&maxIter) < 0) {
		    opserr << "WARNING invalid maxIter\n";
		    return 0;
		}
		if(OPS_GetDoubleInput(&numData,&tol) < 0) {
		    opserr << "WARNING invalid tol\n";
		    return 0;
		}
	    }
	} else if(strcmp(type,"-subdivide") == 0) {
	    if(OPS_GetNumRemainingInputArgs() > 1) {
		if(OPS_GetIntInput(&numData,&numSub) < 0) {
		    opserr << "WARNING invalid numSubdivide\n";
		    return 0;
		}
		if(OPS_GetDoubleInput(&numData,&subFac) < 0) {
		    opserr << "WARNING invalid subdivideFactor\n";
		    return 0;
		}
	    }
	} else if(strcmp(type,"-mass") == 0) {
	    if(OPS_GetNumRemainingInputArgs() > 0) {
		if(OPS_GetDoubleInput(&numData,&mass) < 0) {
		    opserr << "WARNING invalid mass\n";
		    return 0;
		}
	    }
	}
    //Tang.S
    else if (strcmp(type, "-damp") == 0) {

        if (OPS_GetNumRemainingInputArgs() > 0) {
            if (OPS_GetIntInput(&numData, &dampingTag) < 0) return 0;
            theDamping = OPS_getDamping(dampingTag);
            if (theDamping == 0) {
                opserr << "damping not found\n";
                return 0;
            }
        }
    }
    }

    // check transf
    CrdTransf* theTransf = OPS_getCrdTransf(iData[3]);
    if(theTransf == 0) {
	opserr<<"coord transfomration not found\n";
	return 0;
    }

    // check beam integrataion
    BeamIntegrationRule* theRule = OPS_getBeamIntegrationRule(iData[4]);
    if(theRule == 0) {
	opserr<<"beam integration not found\n";
	return 0;
    }
    BeamIntegration* bi = theRule->getBeamIntegration();
    if(bi == 0) {
	opserr<<"beam integration is null\n";
	return 0;
    }

    // check sections
    const ID& secTags = theRule->getSectionTags();
    SectionForceDeformation** sections = new SectionForceDeformation *[secTags.Size()];
    for(int i=0; i<secTags.Size(); i++) {
	sections[i] = OPS_getSectionForceDeformation(secTags(i));
	if(sections[i] == 0) {
	    opserr<<"section "<<secTags(i)<<"not found\n";
	    return 0;
	}
    }

    Element *theEle =  new ForceBeamColumn3d(iData[0],iData[1],iData[2],secTags.Size(),sections,
					     *bi,*theTransf,mass,maxIter,tol,numSub,subFac,theDamping);
    delete [] sections;
    return theEle;
}

void *OPS_ForceBeamColumn3d(const ID &info) {
    // data needed
    int iData[5];
    double mass = 0.0, tol = 1e-12, subFac=10.0;
    int maxIter = 10, numSub = 4;
    int numData;
    int dampingTag = 0;
    Damping* theDamping = 0;

    int ndm = OPS_GetNDM();
    int ndf = OPS_GetNDF();
    if (ndm != 3 || ndf != 6) {
        opserr << "ndm must be 3 and ndf must be 6\n";
        return 0;
    }

    // 1. regular elements
    if (info.Size() == 0) {
        numData = 3;
        if (OPS_GetNumRemainingInputArgs() < numData) {
            opserr << "insufficient "
                      "arguments:eleTag,iNode,jNode\n";
            return 0;
        }
        if (OPS_GetIntInput(&numData, &iData[0]) < 0) {
            opserr << "WARNING invalid int inputs\n";
            return 0;
        }
    }

    // 2. regular elements or save data
    if (info.Size() == 0 || info(0) == 1) {
        numData = 2;
        if (OPS_GetNumRemainingInputArgs() < numData) {
            opserr << "insufficient "
                      "arguments:transfTag,integrationTag\n";
            return 0;
        }
        if (OPS_GetIntInput(&numData, &iData[3]) < 0) {
            opserr << "WARNING invalid int inputs\n";
            return 0;
        }

        numData = 1;
        while (OPS_GetNumRemainingInputArgs() > 0) {
            const char *type = OPS_GetString();
            if (strcmp(type, "-iter") == 0) {
                if (OPS_GetNumRemainingInputArgs() > 1) {
                    if (OPS_GetIntInput(&numData, &maxIter) < 0) {
                        opserr << "WARNING invalid maxIter\n";
                        return 0;
                    }
                    if (OPS_GetDoubleInput(&numData, &tol) < 0) {
                        opserr << "WARNING invalid tol\n";
                        return 0;
                    }
                }
            } else if(strcmp(type,"-subdivide") == 0) {
	      if(OPS_GetNumRemainingInputArgs() > 1) {
		if(OPS_GetIntInput(&numData,&numSub) < 0) {
		  opserr << "WARNING invalid numSubdivide\n";
		  return 0;
		}
		if(OPS_GetDoubleInput(&numData,&subFac) < 0) {
		  opserr << "WARNING invalid subdivideFactor\n";
		  return 0;
		}
	      }
	    }
	    else if (strcmp(type, "-mass") == 0) {
                if (OPS_GetNumRemainingInputArgs() > 0) {
                    if (OPS_GetDoubleInput(&numData, &mass) < 0) {
                        opserr << "WARNING invalid mass\n";
                        return 0;
                    }
                }
        }
        //Tang.S
        else if (strcmp(type, "-damp") == 0) {

                if (OPS_GetNumRemainingInputArgs() > 0) {
                    if (OPS_GetIntInput(&numData, &dampingTag) < 0) return 0;
                    theDamping = OPS_getDamping(dampingTag);
                    if (theDamping == 0) {
                        opserr << "damping not found\n";
                        return 0;
                    }
                }
            }
        }
    }

    // 3: save data
    static std::map<int, Vector> meshdata;
    if (info.Size() > 0 && info(0) == 1) {
        if (info.Size() < 2) {
            opserr << "WARNING: need info -- inmesh, meshtag\n";
            return 0;
        }

        // save the data for a mesh
        Vector &mdata = meshdata[info(1)];
        mdata.resize(7);
        mdata(0) = iData[3];
        mdata(1) = iData[4];
        mdata(2) = mass;
        mdata(3) = tol;
        mdata(4) = maxIter;
	mdata(5) = numSub;
	mdata(6) = subFac;	
        return &meshdata;
    }

    // 4: load data
    if (info.Size() > 0 && info(0) == 2) {
        if (info.Size() < 5) {
            opserr << "WARNING: need info -- inmesh, meshtag, "
                      "eleTag, nd1, nd2\n";
            return 0;
        }

        // get the data for a mesh
        Vector &mdata = meshdata[info(1)];
        if (mdata.Size() < 5) return 0;

        iData[0] = info(2);
        iData[1] = info(3);
        iData[2] = info(4);
        iData[3] = mdata(0);
        iData[4] = mdata(1);
        mass = mdata(2);
        tol = mdata(3);
        maxIter = mdata(4);
	numSub = mdata(5);
	subFac = mdata(6);
    }

    // 5: create element
    CrdTransf *theTransf = OPS_getCrdTransf(iData[3]);
    if (theTransf == 0) {
        opserr << "coord transfomration not found\n";
        return 0;
    }

    // check beam integrataion
    BeamIntegrationRule *theRule =
        OPS_getBeamIntegrationRule(iData[4]);
    if (theRule == 0) {
        opserr << "beam integration not found\n";
        return 0;
    }
    BeamIntegration *bi = theRule->getBeamIntegration();
    if (bi == 0) {
        opserr << "beam integration is null\n";
        return 0;
    }

    // check sections
    const ID &secTags = theRule->getSectionTags();
    SectionForceDeformation **sections =
        new SectionForceDeformation *[secTags.Size()];
    for (int i = 0; i < secTags.Size(); i++) {
        sections[i] = OPS_getSectionForceDeformation(secTags(i));
        if (sections[i] == 0) {
            opserr << "section " << secTags(i) << "not found\n";
            return 0;
        }
    }

    Element *theEle = new ForceBeamColumn3d(
        iData[0], iData[1], iData[2], secTags.Size(), sections, *bi,
        *theTransf, mass, maxIter, tol, numSub, subFac,theDamping);
    delete[] sections;
    return theEle;
}

// constructor:
// invoked by a FEM_ObjectBroker, recvSelf() needs to be invoked on this object.
ForceBeamColumn3d::ForceBeamColumn3d(): 
  Element(0,ELE_TAG_ForceBeamColumn3d), connectedExternalNodes(2), 
  beamIntegr(0), numSections(0), sections(0), crdTransf(0),
  rho(0.0), maxIters(0), tol(0.0),
  initialFlag(0),
  kv(NEBD,NEBD), Se(NEBD),
  kvcommit(NEBD,NEBD), Secommit(NEBD),
  fs(0), vs(0), Ssr(0), vscommit(0),
  numEleLoads(0), sizeEleLoads(0), eleLoads(0), eleLoadFactors(0), load(12),
  Ki(0), isTorsion(false), maxSubdivisions(1), subdivideFactor(1.0), parameterID(0),
  theDamping(0)
{
  load.Zero();

  theNodes[0] = 0;  
  theNodes[1] = 0;
}

// constructor which takes the unique element tag, sections,
// and the node ID's of its nodal end points. 
// allocates the necessary space needed by each object
ForceBeamColumn3d::ForceBeamColumn3d (int tag, int nodeI, int nodeJ,
				      int numSec, SectionForceDeformation **sec,
				      BeamIntegration &bi,
				      CrdTransf &coordTransf, double massDensPerUnitLength,
				      int maxNumIters, double tolerance,
				      int maxNumSub, double subFac,
				      Damping *damping):
  Element(tag,ELE_TAG_ForceBeamColumn3d), connectedExternalNodes(2),
  beamIntegr(0), numSections(0), sections(0), crdTransf(0),
  rho(massDensPerUnitLength),maxIters(maxNumIters), tol(tolerance), 
  initialFlag(0),
  kv(NEBD,NEBD), Se(NEBD), 
  kvcommit(NEBD,NEBD), Secommit(NEBD),
  fs(0), vs(0),Ssr(0), vscommit(0),
  numEleLoads(0), sizeEleLoads(0), eleLoads(0), eleLoadFactors(0), load(12),
  Ki(0), isTorsion(false), maxSubdivisions(maxNumSub), subdivideFactor(subFac), parameterID(0),
  theDamping(0)
{
  if (maxSubdivisions < 1)
    maxSubdivisions = 1;
  if (subdivideFactor < 1.0)
    subdivideFactor = 1.0;
  
  load.Zero();
  
  theNodes[0] = 0;
  theNodes[1] = 0;

  connectedExternalNodes(0) = nodeI;
  connectedExternalNodes(1) = nodeJ;    

  beamIntegr = bi.getCopy();
  if (beamIntegr == 0) {
    opserr << "Error: ForceBeamColumn3d::ForceBeamColumn3d: could not create copy of beam integration object" << endln;
    exit(-1);
  }

  // get copy of the transformation object   
  crdTransf = coordTransf.getCopy3d(); 
  if (crdTransf == 0) {
    opserr << "Error: ForceBeamColumn3d::ForceBeamColumn3d: could not create copy of coordinate transformation object" << endln;
    exit(-1);
  }

  if (damping)
  {
    theDamping =(*damping).getCopy();
    
    if (!theDamping) {
      opserr << "Error: ForceBeamColumn3d::ForceBeamColumn3d: could not create copy of damping object\n";
      exit(-1);
    }
  }

  this->setSectionPointers(numSec, sec);
}

// ~ForceBeamColumn3d():
// 	destructor
//      delete must be invoked on any objects created by the object
ForceBeamColumn3d::~ForceBeamColumn3d()
{
  if (sections != 0) {
    for (int i=0; i < numSections; i++)
      if (sections[i] != 0)
	delete sections[i];
    delete [] sections;
  }
  
  if (sizeEleLoads != 0) {
      if (eleLoads != 0)
          delete[] eleLoads;

      if (eleLoadFactors != 0)
          delete[] eleLoadFactors;
  }

  if (fs != 0) 
    delete [] fs;
  
  if (vs != 0) 
    delete [] vs;
  
  if (Ssr != 0) 
    delete [] Ssr;
  
  if (vscommit != 0) 
    delete [] vscommit;
  
  if (crdTransf != 0)
    delete crdTransf;

  if (beamIntegr != 0)
    delete beamIntegr;
  
  if (Ki != 0)
    delete Ki;

	if (theDamping) delete theDamping;
}

int
ForceBeamColumn3d::getNumExternalNodes(void) const
{
  return 2;
}

const ID &
ForceBeamColumn3d::getExternalNodes(void) 
{
  return connectedExternalNodes;
}

Node **
ForceBeamColumn3d::getNodePtrs()
{
  return theNodes;
}

int
ForceBeamColumn3d::getNumDOF(void) 
{
  return NEGD;
}

void
ForceBeamColumn3d::setDomain(Domain *theDomain)
{
  // check Domain is not null - invoked when object removed from a domain
  if (theDomain == 0) {
    theNodes[0] = 0;
    theNodes[1] = 0;
    
    opserr << "ForceBeamColumn3d::setDomain:  theDomain = 0 ";
    exit(0); 
  }

  // get pointers to the nodes
  
  int Nd1 = connectedExternalNodes(0);  
  int Nd2 = connectedExternalNodes(1);
  
  theNodes[0] = theDomain->getNode(Nd1);
  theNodes[1] = theDomain->getNode(Nd2);  
  
  if (theNodes[0] == 0) {
    opserr << "ForceBeamColumn3d::setDomain: Nd1: ";
    opserr << Nd1 << "does not exist in model\n";
    exit(0);
  }
  
  if (theNodes[1] == 0) {
    opserr << "ForceBeamColumn3d::setDomain: Nd2: ";
    opserr << Nd2 << "does not exist in model\n";
    exit(0);
  }
  
  // call the DomainComponent class method 
  this->DomainComponent::setDomain(theDomain);
  
  // ensure connected nodes have correct number of dof's
  int dofNode1 = theNodes[0]->getNumberDOF();
  int dofNode2 = theNodes[1]->getNumberDOF();
  
  if ((dofNode1 != NND) || (dofNode2 != NND)) {
    opserr << "ForceBeamColumn3d::setDomain(): Nd2 or Nd1 incorrect dof ";
    exit(0);
  }
   
  // initialize the transformation
  if (crdTransf->initialize(theNodes[0], theNodes[1])) {
    opserr << "ForceBeamColumn3d::setDomain(): Error initializing coordinate transformation";  
    exit(0);
  }
    
  // initialize the damping
  if (theDamping && theDamping->setDomain(theDomain, NEBD)) {
    opserr << "ForceBeamColumn3d::setDomain(): Error initializing damping";  
    exit(0);
  }
    
  // get element length
  double L = crdTransf->getInitialLength();
  if (L == 0.0) {
    opserr << "ForceBeamColumn3d::setDomain() -- zero length for element with tag: " << this->getTag();  
    exit(0);
  }

  if (initialFlag == 0) 
    this->initializeSectionHistoryVariables();
}

int
ForceBeamColumn3d::setDamping(Domain *theDomain, Damping *damping)
{
  if (theDomain && damping)
  {
    if (theDamping) delete theDamping;

    theDamping =(*damping).getCopy();
    
    if (!theDamping) {
      opserr << "ForceBeamColumn3d::setDamping -- failed to get copy of damping\n";
      return -1;
    }
    if (theDamping->setDomain(theDomain, NEBD)) {
      opserr << "ForceBeamColumn3d::setDamping -- Error initializing damping\n";
      return -2;
    }
  }
  
  return 0;
}

int
ForceBeamColumn3d::commitState()
{
  int err = 0;
  int i = 0;

  // call element commitState to do any base class stuff
  if ((err = this->Element::commitState()) != 0) {
    opserr << "ForceBeamColumn3d::commitState () - failed in base class";
  }    
  
  do {
    vscommit[i] = vs[i];
    err = sections[i++]->commitState();
    
  } while (err == 0 && i < numSections);
  
  if (err)
    return err;
  
  // commit the transformation between coord. systems
  if ((err = crdTransf->commitState()) != 0)
    return err;
  
  // commit the element variables state
  kvcommit = kv;
  Secommit = Se;
  
  //   initialFlag = 0;  fmk - commented out, see what happens to Example3.1.tcl if uncommented
  //                         - i have not a clue why, ask remo if he ever gets in contact with us again!
  
  if (theDamping && (err = theDamping->commitState()))
    return err;

  return err;
}

int ForceBeamColumn3d::revertToLastCommit()
{
  int err;
  int i = 0;
  
  do {
    vs[i] = vscommit[i];
    err = sections[i]->revertToLastCommit();
    
    sections[i]->setTrialSectionDeformation(vs[i]);      
    Ssr[i] = sections[i]->getStressResultant();
    fs[i]  = sections[i]->getSectionFlexibility();
    
    i++;
  } while (err == 0 && i < numSections);
  
  
  if (err)
    return err;
  
  // revert the transformation to last commit
  if ((err = crdTransf->revertToLastCommit()) != 0)
    return err;
  
  // revert the element state to last commit
  Se   = Secommit;
  kv   = kvcommit;
  
  initialFlag = 0;
  // this->update();
  
  if (theDamping && (err = theDamping->revertToLastCommit()))
    return err;

  return err;
}

int ForceBeamColumn3d::revertToStart()
{
  // revert the sections state to start
  int err;
  int i = 0;
  
  do {
    fs[i].Zero();
    vs[i].Zero();
    Ssr[i].Zero();
    err = sections[i++]->revertToStart();
    
  } while (err == 0 && i < numSections);
  
  if (err)
    return err;
  
  // revert the transformation to start
  if ((err = crdTransf->revertToStart()) != 0)
    return err;
  
  // revert the element state to start
  Secommit.Zero();
  kvcommit.Zero();
  
  Se.Zero();
  kv.Zero();
  
  initialFlag = 0;
  // this->update();

  if (theDamping && (err = theDamping->revertToStart()))
    return err;

  return err;
}


const Matrix &
ForceBeamColumn3d::getInitialStiff(void)
{
  // check for quick return
  if (Ki != 0)
    return *Ki;

  static Matrix f(NEBD,NEBD);   // element flexibility matrix  
  this->getInitialFlexibility(f);
  
  static Matrix I(NEBD,NEBD);   // an identity matrix for matrix inverse  
  I.Zero();
  for (int i=0; i<NEBD; i++)
    I(i,i) = 1.0;
  
  // calculate element stiffness matrix
  // invert3by3Matrix(f, kv);
  static Matrix kvInit(NEBD, NEBD);
  if (f.Solve(I, kvInit) < 0)
    opserr << "ForceBeamColumn3d::getInitialStiff() -- could not invert flexibility for element with tag: " << this->getTag() << endln;

  if(theDamping) kvInit *= theDamping->getStiffnessMultiplier();

    Ki = new Matrix(crdTransf->getInitialGlobalStiffMatrix(kvInit));

    return *Ki;
  }

  const Matrix &
  ForceBeamColumn3d::getTangentStiff(void)
  {
    crdTransf->update();	// Will remove once we clean up the corotational 3d transformation -- MHS
    return crdTransf->getGlobalStiffMatrix(kv, Se);
  }

void
ForceBeamColumn3d::computeReactions(double *p0)
{
  int type;
  double L = crdTransf->getInitialLength();
  
  for (int i = 0; i < numEleLoads; i++) {
    
    double loadFactor = eleLoadFactors[i];
    const Vector &data = eleLoads[i]->getData(type, loadFactor);

    if (type == LOAD_TAG_Beam3dUniformLoad) {
      double wy = data(0)*loadFactor;  // Transverse
      double wz = data(1)*loadFactor;  // Transverse
      double wa = data(2)*loadFactor;  // Axial

      p0[0] -= wa*L;
      double V = 0.5*wy*L;
      p0[1] -= V;
      p0[2] -= V;
      V = 0.5*wz*L;
      p0[3] -= V;
      p0[4] -= V;
    }
    else if (type == LOAD_TAG_Beam3dPartialUniformLoad) {
      double wa = data(2)*loadFactor;  // Axial
      double wy = data(0)*loadFactor;  // Transverse
      double wz = data(1)*loadFactor;  // Transverse
      double a = data(3)*L;
      double b = data(4)*L;

      p0[0] -= wa*(b-a);
      double Fy = wy*(b-a);
      double c = a + 0.5*(b-a);
      p0[1] -= Fy*(1-c/L);
      p0[2] -= Fy*c/L;
      double Fz = wz*(b-a);
      p0[3] -= Fz*(1-c/L);
      p0[4] -= Fz*c/L;      
    }
    else if (type == LOAD_TAG_Beam3dPointLoad) {
      double Py = data(0)*loadFactor;
      double Pz = data(1)*loadFactor;
      double N  = data(2)*loadFactor;
      double aOverL = data(3);
      
      if (aOverL < 0.0 || aOverL > 1.0)
	continue;
      
      double V1 = Py*(1.0-aOverL);
      double V2 = Py*aOverL;
      p0[0] -= N;
      p0[1] -= V1;
      p0[2] -= V2;
      V1 = Pz*(1.0-aOverL);
      V2 = Pz*aOverL;
      p0[3] -= V1;
      p0[4] -= V2;
    }
  }
}

void
ForceBeamColumn3d::computeReactionSensitivity(double *dp0dh, int gradNumber)
{
  int type;
  double L = crdTransf->getInitialLength();
  
  double dLdh = crdTransf->getdLdh();

  for (int i = 0; i < numEleLoads; i++) {
    
    const Vector &data = eleLoads[i]->getData(type, 1.0);

    if (type == LOAD_TAG_Beam3dUniformLoad) {
      double wy = data(0)*1.0;  // Transverse
      double wz = data(1)*1.0;  // Transverse
      double wa = data(2)*1.0;  // Axial

      const Vector &sens = eleLoads[i]->getSensitivityData(gradNumber);
      double dwydh = sens(0);
      double dwzdh = sens(1);
      double dwadh = sens(2);
      
      //p0[0] -= wa*L;
      dp0dh[0] -= wa*dLdh + dwadh*L;

      //double V = 0.5*wy*L;
      //p0[1] -= V;
      //p0[2] -= V;
      double dVdh = 0.5*(wy*dLdh + dwydh*L);
      dp0dh[1] -= dVdh;
      dp0dh[2] -= dVdh;
      dVdh = 0.5*(wz*L + dwzdh*L);
      dp0dh[3] -= dVdh;
      dp0dh[4] -= dVdh;
    }
    else if (type == LOAD_TAG_Beam3dPointLoad) {
      double Py = data(0)*1.0;
      double Pz = data(1)*1.0;
      double N  = data(2)*1.0;
      double aOverL = data(3);

      if (aOverL < 0.0 || aOverL > 1.0)
	continue;
      
      const Vector &sens = eleLoads[i]->getSensitivityData(gradNumber);
      double dPydh = sens(0);
      double dPzdh = sens(1);
      double dNdh  = sens(2);
      double daLdh = sens(3);

      //double a = aOverL*L;
      
      //double V1 = Py*(1.0-aOverL);
      //double V2 = Py*aOverL;
      double dV1dh = Py*(0.0-daLdh) + dPydh*(1.0-aOverL);
      double dV2dh = Py*daLdh + dPydh*aOverL;
      
      //p0[0] -= N;
      //p0[1] -= V1;
      //p0[2] -= V2;
      dp0dh[0] -= dNdh;
      dp0dh[1] -= dV1dh;
      dp0dh[2] -= dV2dh;

      dV1dh = Pz*(0.0-daLdh) + dPzdh*(1.0-aOverL);
      dV2dh = Pz*daLdh + dPzdh*aOverL;
      dp0dh[3] -= dV1dh;
      dp0dh[4] -= dV2dh;
    }
  }
}

const Vector &
ForceBeamColumn3d::getResistingForce(void)
{
  // Will remove once we clean up the corotational 3d transformation -- MHS
  crdTransf->update();
  
  double p0[5];
  Vector p0Vec(p0, 5);
  p0Vec.Zero();
  
  if (numEleLoads > 0)
    this->computeReactions(p0);
  
  theVector =  crdTransf->getGlobalResistingForce(Se, p0Vec);
  
  if (rho != 0)
    theVector.addVector(1.0, load, -1.0);
  
  return theVector;
}

const Vector &
ForceBeamColumn3d::getDampingForce(void)
{
  crdTransf->update();

  return crdTransf->getGlobalResistingForce(theDamping->getDampingForce(), Vector(5));
}

void
  ForceBeamColumn3d::initializeSectionHistoryVariables (void)
{
  for (int i = 0; i < numSections; i++) {
    int order = sections[i]->getOrder();
    
    fs[i] = Matrix(order,order);
    vs[i] = Vector(order);
    Ssr[i] = Vector(order);
    
      vscommit[i] = Vector(order);
  }
}

  /********* NEWTON , SUBDIVIDE AND INITIAL ITERATIONS ********************
   */
  int
  ForceBeamColumn3d::update()
  {
    // if have completed a recvSelf() - do a revertToLastCommit
    // to get Ssr, etc. set correctly
    if (initialFlag == 2)
      this->revertToLastCommit();

    // update the transformation
    crdTransf->update();

    // get basic displacements and increments
    const Vector &v = crdTransf->getBasicTrialDisp();    

    static Vector dv(NEBD);
    dv = crdTransf->getBasicIncrDeltaDisp();    

    if (initialFlag != 0 && dv.Norm() <= DBL_EPSILON && numEleLoads == 0)
      return 0;

    static Vector vin(NEBD);
    vin = v;
    vin -= dv;
    double L = crdTransf->getInitialLength();
    double oneOverL  = 1.0/L;  

    double xi[maxNumSections];
    beamIntegr->getSectionLocations(numSections, L, xi);

    double wt[maxNumSections];
    beamIntegr->getSectionWeights(numSections, L, wt);

    static Vector vr(NEBD);       // element residual displacements
    static Matrix f(NEBD,NEBD);   // element flexibility matrix

    static Matrix I(NEBD,NEBD);   // an identity matrix for matrix inverse
    double dW;                    // section strain energy (work) norm 
    int i, j;

    I.Zero();
    for (i=0; i<NEBD; i++)
      I(i,i) = 1.0;

    int numSubdivide = 1;
    bool converged = false;
    static Vector dSe(NEBD);
    static Vector dvToDo(NEBD);
    static Vector dvTrial(NEBD);
    static Vector SeTrial(NEBD);
    static Matrix kvTrial(NEBD, NEBD);

    dvToDo = dv;
    dvTrial = dvToDo;

    //static double factor = 10;
    double factor = subdivideFactor;
    double dW0 = 0.0;

    //maxSubdivisions = 10;

    // fmk - modification to get compatible ele forces and deformations 
    //   for a change in deformation dV we try first a newton iteration, if
    //   that fails we try an initial flexibility iteration on first iteration 
    //   and then regular newton, if that fails we use the initial flexiblity
    //   for all iterations.
    //
    //   if they both fail we subdivide dV & try to get compatible forces
    //   and deformations. if they work and we have subdivided we apply
    //   the remaining dV.

    while (converged == false && numSubdivide <= maxSubdivisions) {

      // try regular newton (if l==0), or
      // initial tangent iterations (if l==1), or
      // initial tangent on first iteration then regular newton (if l==2)

      for (int l=0; l<3; l++) {

	//      if (l == 1) l = 2;
	SeTrial = Se;
	kvTrial = kv;
	for (i=0; i<numSections; i++) {
	  vsSubdivide[i] = vs[i];
	  fsSubdivide[i] = fs[i];
	  SsrSubdivide[i] = Ssr[i];
	}

	// calculate nodal force increments and update nodal forces      
	// dSe = kv * dv;
	dSe.addMatrixVector(0.0, kvTrial, dvTrial, 1.0);
	SeTrial += dSe;

	if (initialFlag != 2) {

	  int numIters = maxIters;
	  if (l == 1) 
	    numIters = 10*maxIters; // allow 10 times more iterations for initial tangent

	  for (j=0; j <numIters; j++) {

	    // initialize f and vr for integration
	    f.Zero();
	    vr.Zero();

	for (i=0; i<numSections; i++) {
	  
	  int order      = sections[i]->getOrder();
	  const ID &code = sections[i]->getType();
	  
	  static Vector Ss;
	  static Vector dSs;
	  static Vector dvs;
	  static Matrix fb;
	  
	  Ss.setData(workArea, order);
	  dSs.setData(&workArea[order], order);
	  dvs.setData(&workArea[2*order], order);
	  fb.setData(&workArea[3*order], order, NEBD);
	  
	  double xL  = xi[i];
	  double xL1 = xL-1.0;
	  double wtL = wt[i]*L;
	  
	  // calculate total section forces
	  // Ss = b*Se + bp*currDistrLoad;
	  // Ss.addMatrixVector(0.0, b[i], Se, 1.0);
	  int ii;
	  for (ii = 0; ii < order; ii++) {
	    switch(code(ii)) {
	    case SECTION_RESPONSE_P:
	      Ss(ii) = SeTrial(0);
	      break;
	    case SECTION_RESPONSE_MZ:
	      Ss(ii) = xL1*SeTrial(1) + xL*SeTrial(2);
	      break;
	    case SECTION_RESPONSE_VY:
	      Ss(ii) = oneOverL*(SeTrial(1)+SeTrial(2));
	      break;
	    case SECTION_RESPONSE_MY:
	      Ss(ii) = xL1*SeTrial(3) + xL*SeTrial(4);
	      break;
	    case SECTION_RESPONSE_VZ:
	      Ss(ii) = oneOverL*(SeTrial(3)+SeTrial(4));
	      break;
	    case SECTION_RESPONSE_T:
	      Ss(ii) = SeTrial(5);
	      break;
	    default:
	      Ss(ii) = 0.0;
	      break;
	    }
	  }
	  
	  // Add the effects of element loads, if present
          // s = b*q + sp
          if (numEleLoads > 0)
            this->computeSectionForces(Ss, i);
	  
	  // dSs = Ss - Ssr[i];
	  dSs = Ss;
	  dSs.addVector(1.0, SsrSubdivide[i], -1.0);
	  
	  // compute section deformation increments
	  if (l == 0) {
	    
	    //  regular newton 
	    //    vs += fs * dSs;     
	    
	    dvs.addMatrixVector(0.0, fsSubdivide[i], dSs, 1.0);
	    
	  } else if (l == 2) {
	    
	    //  newton with initial tangent if first iteration
	    //    vs += fs0 * dSs;     
	    //  otherwise regular newton 
	    //    vs += fs * dSs;     
	    
	    if (j == 0) {
	      const Matrix &fs0 = sections[i]->getInitialFlexibility();
	      
	      dvs.addMatrixVector(0.0, fs0, dSs, 1.0);
	    } else
	      dvs.addMatrixVector(0.0, fsSubdivide[i], dSs, 1.0);
	    
	  } else {
	    
	    //  newton with initial tangent
	    //    vs += fs0 * dSs;     
	    
	    const Matrix &fs0 = sections[i]->getInitialFlexibility();
	    dvs.addMatrixVector(0.0, fs0, dSs, 1.0);
	  }
	  
	  // set section deformations
	  if (initialFlag != 0)
	    vsSubdivide[i] += dvs;
	  
	  if ( sections[i]->setTrialSectionDeformation(vsSubdivide[i]) < 0) {
	    opserr << "ForceBeamColumn3d::update() - section failed in setTrial\n";
	    return -1;
	  }
	  
	  // get section resisting forces
	  SsrSubdivide[i] = sections[i]->getStressResultant();
	  
	  // get section flexibility matrix
	  // FRANK 
	  fsSubdivide[i] = sections[i]->getSectionFlexibility();
	  
	  /*
	    const Matrix &sectionStiff = sections[i]->getSectionTangent();
	    int n = sectionStiff.noRows();
	    Matrix I(n,n); I.Zero(); for (int l=0; l<n; l++) I(l,l) = 1.0;
	    Matrix sectionFlex(n,n);
	    sectionStiff.SolveSVD(I, sectionFlex, 1.0e-6);
	    fsSubdivide[i] = sectionFlex;	    
	  */
	  
	  // calculate section residual deformations
	  // dvs = fs * (Ss - Ssr);
	  dSs = Ss;
	  dSs.addVector(1.0, SsrSubdivide[i], -1.0);  // dSs = Ss - Ssr[i];
	  
	  dvs.addMatrixVector(0.0, fsSubdivide[i], dSs, 1.0);
	  
	  // integrate element flexibility matrix
	  // f = f + (b^ fs * b) * wtL;
	  //f.addMatrixTripleProduct(1.0, b[i], fs[i], wtL);
	  int jj;
	  const Matrix &fSec = fsSubdivide[i];
	  fb.Zero();
	  double tmp;
	  for (ii = 0; ii < order; ii++) {
	    switch(code(ii)) {
	    case SECTION_RESPONSE_P:
	      for (jj = 0; jj < order; jj++)
		fb(jj,0) += fSec(jj,ii)*wtL;
	      break;
	    case SECTION_RESPONSE_MZ:
	      for (jj = 0; jj < order; jj++) {
		tmp = fSec(jj,ii)*wtL;
		    fb(jj,1) += xL1*tmp;
		    fb(jj,2) += xL*tmp;
		  }
		  break;
		case SECTION_RESPONSE_VY:
		  for (jj = 0; jj < order; jj++) {
		    tmp = oneOverL*fSec(jj,ii)*wtL;
		    fb(jj,1) += tmp;
		    fb(jj,2) += tmp;
		  }
		  break;
		case SECTION_RESPONSE_MY:
		  for (jj = 0; jj < order; jj++) {
		    tmp = fSec(jj,ii)*wtL;
		    fb(jj,3) += xL1*tmp;
		    fb(jj,4) += xL*tmp;
		  }
		  break;
		case SECTION_RESPONSE_VZ:
		  for (jj = 0; jj < order; jj++) {
		    tmp = oneOverL*fSec(jj,ii)*wtL;
		    fb(jj,3) += tmp;
		    fb(jj,4) += tmp;
		  }
		  break;
		case SECTION_RESPONSE_T:
		  for (jj = 0; jj < order; jj++)
		    fb(jj,5) += fSec(jj,ii)*wtL;
		  break;
		default:
		  break;
		}
	      }

	      for (ii = 0; ii < order; ii++) {
		switch (code(ii)) {
		case SECTION_RESPONSE_P:
		  for (jj = 0; jj < NEBD; jj++)
		    f(0,jj) += fb(ii,jj);
		  break;
		case SECTION_RESPONSE_MZ:
		  for (jj = 0; jj < NEBD; jj++) {
		    tmp = fb(ii,jj);
		    f(1,jj) += xL1*tmp;
		    f(2,jj) += xL*tmp;
		  }
		  break;
		case SECTION_RESPONSE_VY:
		  for (jj = 0; jj < NEBD; jj++) {
		    tmp = oneOverL*fb(ii,jj);
		    f(1,jj) += tmp;
		    f(2,jj) += tmp;
		  }
		  break;
		case SECTION_RESPONSE_MY:
		  for (jj = 0; jj < NEBD; jj++) {
		    tmp = fb(ii,jj);
		    f(3,jj) += xL1*tmp;
		    f(4,jj) += xL*tmp;
		  }
		  break;
		case SECTION_RESPONSE_VZ:
		  for (jj = 0; jj < NEBD; jj++) {
		    tmp = oneOverL*fb(ii,jj);
		    f(3,jj) += tmp;
		    f(4,jj) += tmp;
		  }
		  break;
		case SECTION_RESPONSE_T:
		  for (jj = 0; jj < NEBD; jj++)
		    f(5,jj) += fb(ii,jj);
		  break;
		default:
		  break;
		}
	      }

	      // integrate residual deformations
	      // vr += (b^ (vs + dvs)) * wtL;
	      //vr.addMatrixTransposeVector(1.0, b[i], vs[i] + dvs, wtL);
	      dvs.addVector(1.0, vsSubdivide[i], 1.0);
	      double dei;
	      for (ii = 0; ii < order; ii++) {
		dei = dvs(ii)*wtL;
		switch(code(ii)) {
		case SECTION_RESPONSE_P:
		  vr(0) += dei;
		  break;
		case SECTION_RESPONSE_MZ:
		  vr(1) += xL1*dei; vr(2) += xL*dei;
		  break;
		case SECTION_RESPONSE_VY:
		  tmp = oneOverL*dei;
		  vr(1) += tmp; vr(2) += tmp;
		  break;
		case SECTION_RESPONSE_MY:
		  vr(3) += xL1*dei; vr(4) += xL*dei;
		  break;
		case SECTION_RESPONSE_VZ:
		  tmp = oneOverL*dei;
		  vr(3) += tmp; vr(4) += tmp;
		  break;
		case SECTION_RESPONSE_T:
		  vr(5) += dei;
		  break;
		default:
		  break;
		}
	      }
	    }

	    if (!isTorsion) {
	      f(5,5) = DefaultLoverGJ;
	      vr(5) = SeTrial(5)*DefaultLoverGJ;
	    }

	    // calculate element stiffness matrix
	    // invert3by3Matrix(f, kv);	  
	    // FRANK
	    //	  if (f.SolveSVD(I, kvTrial, 1.0e-12) < 0)
	    if (f.Solve(I, kvTrial) < 0)
	      opserr << "ForceBeamColumn3d::update() -- could not invert flexibility for element with tag: " << this->getTag() << endln;;
	    
	    // dv = vin + dvTrial  - vr
	    dv = vin;
	    dv += dvTrial;
	    dv -= vr;

	    // dv.addVector(1.0, vr, -1.0);

	    // dSe = kv * dv;
	    dSe.addMatrixVector(0.0, kvTrial, dv, 1.0);

	    dW = dv ^ dSe; 
	    if (dW0 == 0.0) 
	      dW0 = dW;

	    SeTrial += dSe;

	    // check for convergence of this interval
	    if (fabs(dW) < tol) { 

	      // set the target displacement
	      dvToDo -= dvTrial;
	      vin += dvTrial;

	      // check if we have got to where we wanted
	      if (dvToDo.Norm() <= DBL_EPSILON) {
		converged = true;

	      } else {  // we convreged but we have more to do

		// reset variables for start of next subdivision
		dvTrial = dvToDo;
		numSubdivide = 1;  // NOTE setting subdivide to 1 again maybe too much
	      }

	      // set kv, vs and Se values
	      kv = kvTrial;
	      Se = SeTrial;

	      for (int k=0; k<numSections; k++) {
		vs[k] = vsSubdivide[k];
		fs[k] = fsSubdivide[k];
		Ssr[k] = SsrSubdivide[k];
	      }

	      // break out of j & l loops
	      j = numIters+1;
	      l = 4;

	    } else {   //  if (fabs(dW) < tol) { 

	      // if we have failed to convrege for all of our newton schemes
	      // - reduce step size by the factor specified
	      if (j == (numIters-1) && (l == 2)) {
		dvTrial /= factor;
		numSubdivide++;
	      }
	    }

	  } // for (j=0; j<numIters; j++)
	} // if (initialFlag != 2)
      } // for (int l=0; l<2; l++)
    } // while (converged == false)

  if (theDamping)
  {
    kv *= theDamping->getStiffnessMultiplier();
    theDamping->update(Se);
  }

    // if fail to converge we return an error flag & print an error message

    if (converged == false) {
      opserr << "WARNING - ForceBeamColumn3d::update - failed to get compatible ";
      opserr << "element forces & deformations for element: ";
      opserr << this->getTag() << "(dW: << " << dW << ", dW0: " << dW0 << ")\n";

      /*
      opserr << "Section Tangent Condition Numbers: ";
      for (int i=0; i<numSections; i++) {
	const Matrix &sectionStiff = sections[i]->getSectionTangent();
	double conditionNumber = sectionStiff.conditionNumber();
	opserr << conditionNumber << " ";
      }
      opserr << endln;
      */

      return -1;
    }

    initialFlag = 1;

    return 0;
  }

  void ForceBeamColumn3d::getForceInterpolatMatrix(double xi, Matrix &b, const ID &code)
  {
    b.Zero();

    double L = crdTransf->getInitialLength();
    for (int i = 0; i < code.Size(); i++) {
      switch (code(i)) {
      case SECTION_RESPONSE_MZ:		// Moment, Mz, interpolation
	b(i,1) = xi - 1.0;
	b(i,2) = xi;
	break;
      case SECTION_RESPONSE_P:		// Axial, P, interpolation
	b(i,0) = 1.0;
	break;
      case SECTION_RESPONSE_VY:		// Shear, Vy, interpolation
	b(i,1) = b(i,2) = 1.0/L;
	break;
      case SECTION_RESPONSE_MY:              // Moment, My, interpolation
	b(i,3) = xi - 1.0;
	b(i,4) = xi;
	break;
      case SECTION_RESPONSE_VZ:              // Shear, Vz, interpolation
	b(i,3) = b(i,4) = 1.0/L;
	break;
      case SECTION_RESPONSE_T:               // Torque, T, interpolation
	b(i,5) = 1.0;
	break;
      default:
	break;
      }
    }
  }

  void ForceBeamColumn3d::getDistrLoadInterpolatMatrix(double xi, Matrix &bp, const ID &code)
  {
    bp.Zero();

    double L = crdTransf->getInitialLength();
    for (int i = 0; i < code.Size(); i++) {
      switch (code(i)) {
      case SECTION_RESPONSE_MZ:		// Moment, Mz, interpolation
	bp(i,1) = xi*(xi-1)*L*L/2;
	break;
      case SECTION_RESPONSE_P:		// Axial, P, interpolation
	bp(i,0) = (1-xi)*L;
	break;
      case SECTION_RESPONSE_VY:		// Shear, Vy, interpolation
	bp(i,1) = (xi-0.5)*L;
	break;
      case SECTION_RESPONSE_MY:              // Moment, My, interpolation
	bp(i,2) = xi*(1-xi)*L*L/2;
	break;
      case SECTION_RESPONSE_VZ:              // Shear, Vz, interpolation
	bp(i,2) = (0.5-xi)*L;
	break;
      case SECTION_RESPONSE_T:               // Torsion, T, interpolation
	break;
      default:
	break;
      }
    }
  }

  const Matrix &
  ForceBeamColumn3d::getMass(void)
  { 
    theMatrix.Zero();

    double L = crdTransf->getInitialLength();
    if (rho != 0.0)
      theMatrix(0,0) = theMatrix(1,1) = theMatrix(2,2) =
	theMatrix(6,6) = theMatrix(7,7) = theMatrix(8,8) = 0.5*L*rho;

    return theMatrix;
  }

  void 
  ForceBeamColumn3d::zeroLoad(void)
  {
    load.Zero();
    
    // This is a semi-hack -- MHS
    numEleLoads = 0;
    
    return;
  }

int
ForceBeamColumn3d::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  if (numEleLoads == sizeEleLoads) {

    //
    // create larger arrays, copy old, delete old & set as new
    //

    ElementalLoad ** theNextEleLoads = new ElementalLoad *[sizeEleLoads+1];
    double *theNextEleLoadFactors = new double[sizeEleLoads+1];
    for (int i=0; i<numEleLoads; i++) {
      theNextEleLoads[i] = eleLoads[i];
      theNextEleLoadFactors[i] = eleLoadFactors[i];
    }
    delete [] eleLoads;
    delete [] eleLoadFactors;
    eleLoads = theNextEleLoads;
    eleLoadFactors = theNextEleLoadFactors;  

    // increment array size
    sizeEleLoads+=1;
  }

  eleLoadFactors[numEleLoads] = loadFactor;
  eleLoads[numEleLoads] = theLoad;
  numEleLoads++;

  return 0;
}

void
ForceBeamColumn3d::computeSectionForces(Vector &sp, int isec)
{
  int type;

  double L = crdTransf->getInitialLength();

  double xi[maxNumSections];
  beamIntegr->getSectionLocations(numSections, L, xi);
  double x = xi[isec]*L;

  int order = sections[isec]->getOrder();
  const ID &code = sections[isec]->getType();

  for (int i = 0; i < numEleLoads; i++) {

    double loadFactor = eleLoadFactors[i];
    const Vector &data = eleLoads[i]->getData(type, loadFactor);

    if (type == LOAD_TAG_Beam3dUniformLoad) {
      double wy = data(0)*loadFactor;  // Transverse
      double wz = data(1)*loadFactor;  // Transverse
      double wa = data(2)*loadFactor;  // Axial

      for (int ii = 0; ii < order; ii++) {
	
	switch(code(ii)) {
	case SECTION_RESPONSE_P:
	  sp(ii) += wa*(L-x);
	  break;
	case SECTION_RESPONSE_MZ:
	  sp(ii) += wy*0.5*x*(x-L);
	  break;
	case SECTION_RESPONSE_VY:
	  sp(ii) += wy*(x-0.5*L);
	  break;
	case SECTION_RESPONSE_MY:
	  sp(ii) += wz*0.5*x*(L-x);
	  break;
	case SECTION_RESPONSE_VZ:
	  sp(ii) += wz*(0.5*L-x);
	  break;
	default:
	  break;
	}
      }
    }
    else if (type == LOAD_TAG_Beam3dPartialUniformLoad) {
      double wa = data(2)*loadFactor;  // Axial
      double wy = data(0)*loadFactor;  // Transverse
      double wz = data(1)*loadFactor;  // Transverse
      double a = data(3)*L;
      double b = data(4)*L;

      double Fa = wa*(b-a); // resultant axial load
      double Fy = wy*(b-a); // resultant transverse load
      double Fz = wz*(b-a); // resultant transverse load
      double c = a + 0.5*(b-a);
      double VyI = Fy*(1-c/L);
      double VyJ = Fy*c/L;
      double VzI = Fz*(1-c/L);
      double VzJ = Fz*c/L;      

      for (int ii = 0; ii < order; ii++) {
	
	if (x <= a) {
	  switch(code(ii)) {
	  case SECTION_RESPONSE_P:
	    sp(ii) += Fa;
	    break;
	  case SECTION_RESPONSE_MZ:
	    sp(ii) -= VyI*x;
	    break;
	  case SECTION_RESPONSE_MY:
	    sp(ii) += VzI*x;
	    break;	    
	  case SECTION_RESPONSE_VY:
	    sp(ii) -= VyI;
	    break;
	  case SECTION_RESPONSE_VZ:
	    sp(ii) -= VzI;
	    break;	    
	  default:
	    break;
	  }
	}
	else if (x >= b) {
	  switch(code(ii)) {
	  case SECTION_RESPONSE_MZ:
	    sp(ii) += VyJ*(x-L);
	    break;
	  case SECTION_RESPONSE_MY:
	    sp(ii) -= VzJ*(x-L);
	    break;	    
	  case SECTION_RESPONSE_VY:
	    sp(ii) += VyJ;
	    break;
	  case SECTION_RESPONSE_VZ:
	    sp(ii) += VzJ;	    
	    break;
	  default:
	    break;
	  }
	}
	else {
	  switch(code(ii)) {
	  case SECTION_RESPONSE_P:
	    sp(ii) += Fa-wa*(x-a);
	    break;
	  case SECTION_RESPONSE_MZ:
	    sp(ii) += -VyI*x + 0.5*wy*x*x + wy*a*(0.5*a-x);
	    break;
	  case SECTION_RESPONSE_MY:
	    sp(ii) += VzI*x - 0.5*wz*x*x - wz*a*(0.5*a-x);
	    break;	    
	  case SECTION_RESPONSE_VY:
	    sp(ii) += -VyI + wy*(x-a);
	    break;
	  case SECTION_RESPONSE_VZ:
	    sp(ii) += -VzI + wz*(x-a);	    
	    break;
	  default:
	    break;
	  }
	}
      }
    }
    else if (type == LOAD_TAG_Beam3dPointLoad) {
      double Py = data(0)*loadFactor;
      double Pz = data(1)*loadFactor;
      double N  = data(2)*loadFactor;
      double aOverL = data(3);
      
      if (aOverL < 0.0 || aOverL > 1.0)
	continue;
      
      double a = aOverL*L;
      
      double Vy1 = Py*(1.0-aOverL);
      double Vy2 = Py*aOverL;
      
      double Vz1 = Pz*(1.0-aOverL);
      double Vz2 = Pz*aOverL;

      for (int ii = 0; ii < order; ii++) {
	
	if (x <= a) {
	  switch(code(ii)) {
	  case SECTION_RESPONSE_P:
	    sp(ii) += N;
	    break;
	  case SECTION_RESPONSE_MZ:
	    sp(ii) -= x*Vy1;
	    break;
	  case SECTION_RESPONSE_VY:
	    sp(ii) -= Vy1;
	    break;
	  case SECTION_RESPONSE_MY:
	    sp(ii) += x*Vz1;
	    break;
	  case SECTION_RESPONSE_VZ:
	    sp(ii) -= Vz1;
	    break;
	  default:
	    break;
	  }
	}
	else {
	  switch(code(ii)) {
	  case SECTION_RESPONSE_MZ:
	    sp(ii) -= (L-x)*Vy2;
	    break;
	  case SECTION_RESPONSE_VY:
	    sp(ii) += Vy2;
	    break;
	  case SECTION_RESPONSE_MY:
	    sp(ii) += (L-x)*Vz2;
	    break;
	  case SECTION_RESPONSE_VZ:
	    sp(ii) += Vz2;
	    break;
	  default:
	    break;
	  }
	}
      }
    }
    else {
      opserr << "ForceBeamColumn3d::addLoad -- load type unknown for element with tag: " <<
	this->getTag() << endln;
    }
  }
  
  // Don't think we need to do this anymore -- MHS
  //this->update(); // quick fix -- MHS
}

void
ForceBeamColumn3d::computeSectionForceSensitivity(Vector &dspdh, int isec,
						  int gradNumber)
{
  int type;

  double L = crdTransf->getInitialLength();
  double dLdh = crdTransf->getdLdh();

  double xi[maxNumSections];
  beamIntegr->getSectionLocations(numSections, L, xi);

  double dxidh[maxNumSections];
  beamIntegr->getLocationsDeriv(numSections, L, dLdh, dxidh);

  double x = xi[isec]*L;
  double dxdh = xi[isec]*dLdh + dxidh[isec]*L;
      
  int order = sections[isec]->getOrder();
  const ID &code = sections[isec]->getType();

  for (int i = 0; i < numEleLoads; i++) {

    const Vector &data = eleLoads[i]->getData(type, 1.0);
    
    if (type == LOAD_TAG_Beam3dUniformLoad) {
      double wy = data(0)*1.0;  // Transverse
      double wz = data(1)*1.0;  // Transverse
      double wa = data(2)*1.0;  // Axial

      const Vector &sens = eleLoads[i]->getSensitivityData(gradNumber);
      double dwydh = sens(0);
      double dwzdh = sens(1);
      double dwadh = sens(2);

      for (int ii = 0; ii < order; ii++) {
	
	switch(code(ii)) {
	case SECTION_RESPONSE_P:
	  //sp(ii) += wa*(L-x);
	  dspdh(ii) += dwadh*(L-x) + wa*(dLdh-dxdh);
	  break;
	case SECTION_RESPONSE_MZ:
	  //sp(ii) += wy*0.5*x*(x-L);
	  //dspdh(ii) += 0.5*(dwydh*x*(x-L) + wy*dxdh*(x-L) + wy*x*(dxdh-dLdh));
	  dspdh(ii) += 0.5*(dwydh*x*(x-L) + wy*(dxdh*(2*x-L)-x*dLdh));
	  break;
	case SECTION_RESPONSE_VY:
	  //sp(ii) += wy*(x-0.5*L);
	  dspdh(ii) += dwydh*(x-0.5*L) + wy*(dxdh-0.5*dLdh);
	  break;
    case SECTION_RESPONSE_MY:
        //sp(ii) += wz*0.5*x*(L-x);
        //dspdh(ii) += 0.5*(dwzdh*x*(L-x) + wz*dxdh*(L-x) + wz*x*(dLdh-dxdh));
        dspdh(ii) += 0.5*(dwzdh*x*(L-x) + wz*(dxdh*(L-2*x) + x*dLdh));
        break;
    case SECTION_RESPONSE_VZ:
        //sp(ii) += wz*(x-0.5*L);
        dspdh(ii) += dwzdh*(0.5*L-x) + wz*(0.5*dLdh-dxdh);
        break;
    default:
	  break;
	}
      }
    }
    else if (type == LOAD_TAG_Beam3dPointLoad) {
      double Py = data(0)*1.0;
      double Pz = data(1)*1.0;
      double N = data(2)*1.0;
      double aOverL = data(3);
      
      if (aOverL < 0.0 || aOverL > 1.0)
	continue;
      
      const Vector &sens = eleLoads[i]->getSensitivityData(gradNumber);
      double dPydh = sens(0);
      double dPzdh = sens(1);
      double dNdh = sens(2);
      double daLdh = sens(3);

      double a = aOverL*L;

      double Vy1 = Py*(1.0-aOverL);
      double Vy2 = Py*aOverL;
      double dVy1dh = Py*(0.0-daLdh) + dPydh*(1.0-aOverL);
      double dVy2dh = Py*daLdh + dPydh*aOverL;

      double Vz1 = Pz*(1.0 - aOverL);
      double Vz2 = Pz*aOverL;
      double dVz1dh = Pz*(0.0 - daLdh) + dPzdh*(1.0 - aOverL);
      double dVz2dh = Pz*daLdh + dPzdh*aOverL;

      for (int ii = 0; ii < order; ii++) {
	
	if (x <= a) {
	  switch(code(ii)) {
	  case SECTION_RESPONSE_P:
	    //sp(ii) += N;
	    dspdh(ii) += dNdh;
	    break;
	  case SECTION_RESPONSE_MZ:
	    //sp(ii) -= x*Vy1;
	    dspdh(ii) -= (dxdh*Vy1 + x*dVy1dh);
	    break;
	  case SECTION_RESPONSE_VY:
	    //sp(ii) -= Vy1;
	    dspdh(ii) -= dVy1dh;
	    break;
      case SECTION_RESPONSE_MY:
          //sp(ii) += x*Vz1;
          dspdh(ii) += (dxdh*Vz1 + x*dVz1dh);
          break;
      case SECTION_RESPONSE_VZ:
          //sp(ii) -= Vz1;
          dspdh(ii) -= dVz1dh;
          break;
      default:
	    break;
	  }
	}
	else {
	  switch(code(ii)) {
	  case SECTION_RESPONSE_MZ:
	    //sp(ii) -= (L-x)*Vy2;
	    dspdh(ii) -= (dLdh-dxdh)*Vy2 + (L-x)*dVy2dh;
	    break;
	  case SECTION_RESPONSE_VY:
	    //sp(ii) += Vy2;
	    dspdh(ii) += dVy2dh;
	    break;
      case SECTION_RESPONSE_MY:
          //sp(ii) += (L-x)*Vz2;
          dspdh(ii) += (dLdh-dxdh)*Vz2 + (L-x)*dVz2dh;
          break;
      case SECTION_RESPONSE_VZ:
          //sp(ii) += Vz2;
          dspdh(ii) += dVz2dh;
          break;
      default:
	    break;
	  }
	}
      }
    }
    else {
      opserr << "ForceBeamColumn3d::computeSectionForceSensitivity -- load type unknown for element with tag: " <<
	this->getTag() << endln;
    }
  }
}

  int 
  ForceBeamColumn3d::addInertiaLoadToUnbalance(const Vector &accel)
  {
    // Check for a quick return
    if (rho == 0.0)
      return 0;

    // get R * accel from the nodes
    const Vector &Raccel1 = theNodes[0]->getRV(accel);
    const Vector &Raccel2 = theNodes[1]->getRV(accel);    

    double L = crdTransf->getInitialLength();
    double m = 0.5*rho*L;

    load(0) -= m*Raccel1(0);
    load(1) -= m*Raccel1(1);
    load(2) -= m*Raccel1(2);
    load(6) -= m*Raccel2(0);
    load(7) -= m*Raccel2(1);
    load(8) -= m*Raccel2(2);

    return 0;
  }

  const Vector &
  ForceBeamColumn3d::getResistingForceIncInertia()
  {	
    // Compute the current resisting force
    theVector = this->getResistingForce();

    if (theDamping) theVector += this->getDampingForce();

    if (rho != 0.0) {
      const Vector &accel1 = theNodes[0]->getTrialAccel();
      const Vector &accel2 = theNodes[1]->getTrialAccel();

      double L = crdTransf->getInitialLength();
      double m = 0.5*rho*L;

      theVector(0) += m*accel1(0);
      theVector(1) += m*accel1(1);
      theVector(2) += m*accel1(2);
      theVector(6) += m*accel2(0);
      theVector(7) += m*accel2(1);
      theVector(8) += m*accel2(2);

      // add the damping forces if rayleigh damping
      if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
	theVector += this->getRayleighDampingForces();

    } else {

      // add the damping forces if rayleigh damping
      if (betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
	theVector += this->getRayleighDampingForces();
    }

    return theVector;
  }

  int
  ForceBeamColumn3d::sendSelf(int commitTag, Channel &theChannel)
  {  
    // place the integer data into an ID
    int dbTag = this->getDbTag();
    int i, j , k;
    int loc = 0;

    static ID idData(15);  
    idData(0) = this->getTag();
    idData(1) = connectedExternalNodes(0);
    idData(2) = connectedExternalNodes(1);
    idData(3) = numSections;
    idData(4) = maxIters;
    idData(5) = initialFlag;
    idData(6) = (isTorsion) ? 1 : 0;
    idData(13) = maxSubdivisions;
    
    idData(7) = crdTransf->getClassTag();
    int crdTransfDbTag  = crdTransf->getDbTag();
    if (crdTransfDbTag  == 0) {
      crdTransfDbTag = theChannel.getDbTag();
      if (crdTransfDbTag  != 0) 
	crdTransf->setDbTag(crdTransfDbTag);
    }
    idData(8) = crdTransfDbTag;


    idData(9) = beamIntegr->getClassTag();
    int beamIntegrDbTag  = beamIntegr->getDbTag();
    if (beamIntegrDbTag  == 0) {
      beamIntegrDbTag = theChannel.getDbTag();
      if (beamIntegrDbTag  != 0) 
	beamIntegr->setDbTag(beamIntegrDbTag);
    }
    idData(10) = beamIntegrDbTag;

    idData(11) = 0;
    idData(12) = 0;
    if (theDamping) {
      idData(11) = theDamping->getClassTag();
      int dbTag = theDamping->getDbTag();
      if (dbTag == 0) {
        dbTag = theChannel.getDbTag();
        if (dbTag != 0)
	        theDamping->setDbTag(dbTag);
	    }
      idData(12) = dbTag;
    }

    if (theChannel.sendID(dbTag, commitTag, idData) < 0) {
      opserr << "ForceBeamColumn3d::sendSelf() - failed to send ID data\n";
      return -1;
    }    

    // send the coordinate transformation

    if (crdTransf->sendSelf(commitTag, theChannel) < 0) {
      opserr << "ForceBeamColumn3d::sendSelf() - failed to send crdTranf\n";
      return -1;
    }      

    if (beamIntegr->sendSelf(commitTag, theChannel) < 0) {
      opserr << "ForceBeamColumn3d::sendSelf() - failed to send beamIntegr\n";
      return -1;
    }      

    //
    // send an ID for the sections containing each sections dbTag and classTag
    // if section ha no dbTag get one and assign it
    //

    ID idSections(2*numSections);
    loc = 0;
    for (i = 0; i<numSections; i++) {
      int sectClassTag = sections[i]->getClassTag();
      int sectDbTag = sections[i]->getDbTag();
      if (sectDbTag == 0) {
	sectDbTag = theChannel.getDbTag();
	sections[i]->setDbTag(sectDbTag);
      }

      idSections(loc) = sectClassTag;
      idSections(loc+1) = sectDbTag;
      loc += 2;
    }

    if (theChannel.sendID(dbTag, commitTag, idSections) < 0)  {
      opserr << "ForceBeamColumn3d::sendSelf() - failed to send ID data\n";
      return -1;
    }    
 
    //
    // send the sections
    //

    for (j = 0; j<numSections; j++) {
      if (sections[j]->sendSelf(commitTag, theChannel) < 0) {
	opserr << "ForceBeamColumn3d::sendSelf() - section " <<
	  j << "failed to send itself\n";
	return -1;
      }
    }

    // into a vector place distrLoadCommit, rho, UeCommit, Secommit and kvcommit
    int secDefSize = 0;
    for (i = 0; i < numSections; i++) {
       int size = sections[i]->getOrder();
       secDefSize   += size;
    }

    Vector dData(1+1+1+NEBD+NEBD*NEBD+secDefSize + 4); 
    loc = 0;

    // place double variables into Vector
    dData(loc++) = rho;
    dData(loc++) = tol;
    dData(loc++) = subdivideFactor;
  
    // put  distrLoadCommit into the Vector
    //  for (i=0; i<NL; i++) 
    //dData(loc++) = distrLoadcommit(i);

    // place kvcommit into vector
    for (i=0; i<NEBD; i++) 
      dData(loc++) = Secommit(i);

    // place kvcommit into vector
    for (i=0; i<NEBD; i++) 
       for (j=0; j<NEBD; j++)
	  dData(loc++) = kvcommit(i,j);

    // place vscommit into vector
    for (k=0; k<numSections; k++)
       for (i=0; i<sections[k]->getOrder(); i++)
	  dData(loc++) = (vscommit[k])(i);

    // send damping coefficients
    dData(loc++) = alphaM;
    dData(loc++) = betaK;
    dData(loc++) = betaK0;
    dData(loc++) = betaKc;
    
    if (theChannel.sendVector(dbTag, commitTag, dData) < 0) {
       opserr << "ForceBeamColumn3d::sendSelf() - failed to send Vector data\n";

       return -1;
    }    

    // Ask the Damping to send itself
    if (theDamping && theDamping->sendSelf(commitTag, theChannel) < 0) {
        opserr << "ForceBeamColumn3d::sendSelf -- could not send Damping\n";
        return -1;
    }

    return 0;
  }    

  int
  ForceBeamColumn3d::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
  {
    //
    // receive the integer data containing tag, numSections and coord transformation info
    //
    int dbTag = this->getDbTag();
    int i,j,k;

    static ID idData(15); // one bigger than needed 

    if (theChannel.recvID(dbTag, commitTag, idData) < 0)  {
      opserr << "ForceBeamColumn3d::recvSelf() - failed to recv ID data\n";

      return -1;
    }    
 
    this->setTag(idData(0));
    connectedExternalNodes(0) = idData(1);
    connectedExternalNodes(1) = idData(2);
    maxIters = idData(4);
    initialFlag = idData(5);
    isTorsion = (idData(6) == 1) ? true : false;
    maxSubdivisions = idData(13);
    
    int crdTransfClassTag = idData(7);
    int crdTransfDbTag = idData(8);

    int beamIntegrClassTag = idData(9);
    int beamIntegrDbTag = idData(10);

    // create a new crdTransf object if one needed
    if (crdTransf == 0 || crdTransf->getClassTag() != crdTransfClassTag) {
	if (crdTransf != 0)
	    delete crdTransf;

	crdTransf = theBroker.getNewCrdTransf(crdTransfClassTag);

	if (crdTransf == 0) {
	    opserr << "ForceBeamColumn3d::recvSelf() - failed to obtain a CrdTrans object with classTag" <<
	      crdTransfClassTag << endln;
	    return -2;	  
	}
    }

    crdTransf->setDbTag(crdTransfDbTag);
    // invoke recvSelf on the crdTransf obkject
    if (crdTransf->recvSelf(commitTag, theChannel, theBroker) < 0)  
    {
       opserr << "ForceBeamColumn3d::sendSelf() - failed to recv crdTranf\n";
       return -3;
    }      

    // create a new beamIntegr object if one needed
    if (beamIntegr == 0 || beamIntegr->getClassTag() != beamIntegrClassTag) {
	if (beamIntegr != 0)
	    delete beamIntegr;

	beamIntegr = theBroker.getNewBeamIntegration(beamIntegrClassTag);

	if (beamIntegr == 0) {opserr << "ForceBeamColumn3d::recvSelf() - failed to obtain the beam integration object with classTag" <<
	    beamIntegrClassTag << endln;
	  exit(-1);
	}
    }

    beamIntegr->setDbTag(beamIntegrDbTag);

    // invoke recvSelf on the beamIntegr object
    if (beamIntegr->recvSelf(commitTag, theChannel, theBroker) < 0)  
    {
       opserr << "ForceBeamColumn3d::sendSelf() - failed to recv beam integration\n";

       return -3;
    }      

    //
    // recv an ID for the sections containing each sections dbTag and classTag
    //

    ID idSections(2*idData(3));
    int loc = 0;

    if (theChannel.recvID(dbTag, commitTag, idSections) < 0)  {
      opserr << "ForceBeamColumn3d::recvSelf() - failed to recv ID data\n";
      return -1;
    }    

    //
    // now receive the sections
    //
    if (numSections != idData(3)) {

      //
      // we do not have correct number of sections, must delete the old and create
      // new ones before can recvSelf on the sections
      //

      // delete the old
      if (numSections != 0) {
	for (int i=0; i<numSections; i++)
	  delete sections[i];
	delete [] sections;
      }

      // create a section and recvSelf on it
      numSections = idData(3);

      // Delete the old
      if (vscommit != 0)
	delete [] vscommit;

      // Allocate the right number
      vscommit = new Vector[numSections];
      if (vscommit == 0) {
	opserr << "ForceBeamColumn3d::recvSelf -- failed to allocate vscommit array\n";
	return -1;
      }

      // Delete the old
      if (fs != 0)
	delete [] fs;

      // Allocate the right number
      fs = new Matrix[numSections];  
      if (fs == 0) {
	opserr << "ForceBeamColumn3d::recvSelf -- failed to allocate fs array\n";
	return -1;
      }

      // Delete the old
      if (vs != 0)
	delete [] vs;

      // Allocate the right number
      vs = new Vector[numSections];  
      if (vs == 0) {
	opserr << "ForceBeamColumn3d::recvSelf -- failed to allocate vs array\n";
	return -1;
      }

      // Delete the old
      if (Ssr != 0)
	delete [] Ssr;

      // Allocate the right number
      Ssr = new Vector[numSections];  
      if (Ssr == 0) {
	opserr << "ForceBeamColumn3d::recvSelf -- failed to allocate Ssr array\n";

	return -1;
      }

      // create a new array to hold pointers
      sections = new SectionForceDeformation *[idData(3)];
      if (sections == 0) {
	opserr << "ForceBeamColumn3d::recvSelf() - out of memory creating sections array of size" <<
	  idData(3) << endln;
	exit(-1);
      }    

      loc = 0;

      for (i=0; i<numSections; i++) {
	int sectClassTag = idSections(loc);
	int sectDbTag = idSections(loc+1);
	loc += 2;
	sections[i] = theBroker.getNewSection(sectClassTag);
	if (sections[i] == 0) {
	  opserr << "ForceBeamColumn3d::recvSelf() - Broker could not create Section of class type" <<
	    sectClassTag << endln;
	  exit(-1);
	}
	
	sections[i]->setDbTag(sectDbTag);
	if (sections[i]->recvSelf(commitTag, theChannel, theBroker) < 0) {
	  opserr << "ForceBeamColumn3d::recvSelf() - section " << 
	    i << " failed to recv itself\n";
	  return -1;
	} 
      }

      this->initializeSectionHistoryVariables();

    } else {

      // 
      // for each existing section, check it is of correct type
      // (if not delete old & create a new one) then recvSelf on it
      //

      loc = 0;
      for (i=0; i<numSections; i++) {
	int sectClassTag = idSections(loc);
	int sectDbTag = idSections(loc+1);
	loc += 2;

	// check of correct type
	if (sections[i]->getClassTag() !=  sectClassTag) {
	  // delete the old section[i] and create a new one
	  delete sections[i];
	  sections[i] = theBroker.getNewSection(sectClassTag);
	  if (sections[i] == 0) {
	    opserr << "ForceBeamColumn3d::recvSelf() - Broker could not create Section of class type " 
		   << sectClassTag << endln;;
	    exit(-1);
	  }
	}

	// recvvSelf on it
	sections[i]->setDbTag(sectDbTag);
	if (sections[i]->recvSelf(commitTag, theChannel, theBroker) < 0) {
	  opserr << "ForceBeamColumn3d::recvSelf() - section " <<
	    i << " failed to recv itself\n";
	  return -1;
	}     
      }
    }

    // into a vector place distrLoadCommit, rho, UeCommit, Secommit and kvcommit
    int secDefSize = 0;
    for (int ii = 0; ii < numSections; ii++) {
       int size = sections[ii]->getOrder();
       secDefSize   += size;
    }

    Vector dData(1+1+NEBD+NEBD*NEBD+secDefSize+4);   

    if (theChannel.recvVector(dbTag, commitTag, dData) < 0)  {
      opserr << "ForceBeamColumn3d::recvSelf() - failed to send Vector data\n";
      return -1;
    }    

    loc = 0;

    // place double variables into Vector
    rho = dData(loc++);
    tol = dData(loc++);
    subdivideFactor = dData(loc++);
  
    // put  distrLoadCommit into the Vector
    //for (i=0; i<NL; i++) 
    // distrLoad(i) = dData(loc++);

    // place kvcommit into vector
    for (i=0; i<NEBD; i++) 
      Secommit(i) = dData(loc++);

    // place kvcommit into vector
    for (i=0; i<NEBD; i++) 
       for (j=0; j<NEBD; j++)
	  kvcommit(i,j) = dData(loc++);

    kv   = kvcommit;
    Se   = Secommit;

    for (k = 0; k < numSections; k++) {
      int order = sections[k]->getOrder();

      // place vscommit into vector
      vscommit[k] = Vector(order);
      for (i = 0; i < order; i++)
	(vscommit[k])(i) = dData(loc++);
    }

    // set damping coefficients
    alphaM = dData(loc++);
    betaK = dData(loc++);
    betaK0 = dData(loc++);
    betaKc = dData(loc++);
    
    initialFlag = 2;  

    // Check if the Damping is null; if so, get a new one
    int dmpTag = (int)idData(11);
    if (dmpTag) {
      if (theDamping == 0) {
        theDamping = theBroker.getNewDamping(dmpTag);
        if (theDamping == 0) {
          opserr << "ForceBeamColumn3d::recvSelf -- could not get a Damping\n";
          exit(-1);
        }
      }
    
      // Check that the Damping is of the right type; if not, delete
      // the current one and get a new one of the right type
      if (theDamping->getClassTag() != dmpTag) {
        delete theDamping;
        theDamping = theBroker.getNewDamping(dmpTag);
        if (theDamping == 0) {
          opserr << "ForceBeamColumn3d::recvSelf -- could not get a Damping\n";
          exit(-1);
        }
      }
    
      // Now, receive the Damping
      theDamping->setDbTag((int)idData(12));
      if (theDamping->recvSelf(commitTag, theChannel, theBroker) < 0) {
        opserr << "ForceBeamColumn3d::recvSelf -- could not receive Damping\n";
        return -1;
      }
    }
    else {
      if (theDamping) {
        delete theDamping;
        theDamping = 0;
      }
    }
    
    return 0;
  }

  int
  ForceBeamColumn3d::getInitialFlexibility(Matrix &fe)
  {
    fe.Zero();

    double L = crdTransf->getInitialLength();
    double oneOverL  = 1.0/L;  

    double xi[maxNumSections];
    beamIntegr->getSectionLocations(numSections, L, xi);

    double wt[maxNumSections];
    beamIntegr->getSectionWeights(numSections, L, wt);

    for (int i = 0; i < numSections; i++) {

      int order      = sections[i]->getOrder();
      const ID &code = sections[i]->getType();

      Matrix fb(workArea, order, NEBD);

      double xL  = xi[i];
      double xL1 = xL-1.0;
      double wtL = wt[i]*L;

      const Matrix &fSec = sections[i]->getInitialFlexibility();
      fb.Zero();
      double tmp;
      int ii, jj;
      for (ii = 0; ii < order; ii++) {
	switch(code(ii)) {
	case SECTION_RESPONSE_P:
	  for (jj = 0; jj < order; jj++)
	    fb(jj,0) += fSec(jj,ii)*wtL;
	  break;
	case SECTION_RESPONSE_MZ:
	  for (jj = 0; jj < order; jj++) {
	    tmp = fSec(jj,ii)*wtL;
	    fb(jj,1) += xL1*tmp;
	    fb(jj,2) += xL*tmp;
	  }
	  break;
	case SECTION_RESPONSE_VY:
	  for (jj = 0; jj < order; jj++) {
	    tmp = oneOverL*fSec(jj,ii)*wtL;
	    fb(jj,1) += tmp;
	    fb(jj,2) += tmp;
	  }
	  break;
	case SECTION_RESPONSE_MY:
	  for (jj = 0; jj < order; jj++) {
	    tmp = fSec(jj,ii)*wtL;
	    fb(jj,3) += xL1*tmp;
	    fb(jj,4) += xL*tmp;
	  }
	  break;
	case SECTION_RESPONSE_VZ:
	  for (jj = 0; jj < order; jj++) {
	    tmp = oneOverL*fSec(jj,ii)*wtL;
	    fb(jj,3) += tmp;
	    fb(jj,4) += tmp;
	  }
	  break;
	case SECTION_RESPONSE_T:
	  for (jj = 0; jj < order; jj++)
	    fb(jj,5) += fSec(jj,ii)*wtL;
	  break;
	default:
	  break;
	}
      }
      for (ii = 0; ii < order; ii++) {
	switch (code(ii)) {
	case SECTION_RESPONSE_P:
	  for (jj = 0; jj < NEBD; jj++)
	    fe(0,jj) += fb(ii,jj);
	  break;
	case SECTION_RESPONSE_MZ:
	  for (jj = 0; jj < NEBD; jj++) {
	    tmp = fb(ii,jj);
	    fe(1,jj) += xL1*tmp;
	    fe(2,jj) += xL*tmp;
	  }
	  break;
	case SECTION_RESPONSE_VY:
	  for (jj = 0; jj < NEBD; jj++) {
	    tmp = oneOverL*fb(ii,jj);
	    fe(1,jj) += tmp;
	    fe(2,jj) += tmp;
	  }
	  break;
	case SECTION_RESPONSE_MY:
	  for (jj = 0; jj < NEBD; jj++) {
	    tmp = fb(ii,jj);
	    fe(3,jj) += xL1*tmp;
	    fe(4,jj) += xL*tmp;
	  }
	  break;
	case SECTION_RESPONSE_VZ:
	  for (jj = 0; jj < NEBD; jj++) {
	    tmp = oneOverL*fb(ii,jj);
	    fe(3,jj) += tmp;
	    fe(4,jj) += tmp;
	  }
	  break;
	case SECTION_RESPONSE_T:
	  for (jj = 0; jj < NEBD; jj++)
	    fe(5,jj) += fb(ii,jj);
	  break;
	default:
	  break;
	}
      }
    }

    if (!isTorsion)
      fe(5,5) = DefaultLoverGJ;

    return 0;
  }

int
ForceBeamColumn3d::getInitialDeformations(Vector &v0)
{
  v0.Zero();
  if (numEleLoads < 1)
      return 0;

  double L = crdTransf->getInitialLength();
  double oneOverL = 1.0 / L;

  double xi[maxNumSections];
  beamIntegr->getSectionLocations(numSections, L, xi);

  double wt[maxNumSections];
  beamIntegr->getSectionWeights(numSections, L, wt);

  for (int i = 0; i < numSections; i++) {

      int order = sections[i]->getOrder();
      const ID &code = sections[i]->getType();

      double xL = xi[i];
      double xL1 = xL - 1.0;
      double wtL = wt[i] * L;

      static Vector sp;
      sp.setData(workArea, order);
      sp.Zero();

      this->computeSectionForces(sp, i);

      const Matrix &fse = sections[i]->getInitialFlexibility();

      static Vector e;
      e.setData(&workArea[order], order);

      e.addMatrixVector(0.0, fse, sp, 1.0);

      double dei, tmp;
      for (int ii = 0; ii < order; ii++) {
          dei = e(ii)*wtL;
          switch (code(ii)) {
          case SECTION_RESPONSE_P:
              v0(0) += dei;
              break;
          case SECTION_RESPONSE_MZ:
              v0(1) += xL1*dei; v0(2) += xL*dei;
              break;
          case SECTION_RESPONSE_VY:
              tmp = oneOverL*dei;
              v0(1) += tmp; v0(2) += tmp;
              break;
          case SECTION_RESPONSE_MY:
              v0(3) += xL1*dei; v0(4) += xL*dei;
              break;
          case SECTION_RESPONSE_VZ:
              tmp = oneOverL*dei;
              v0(3) += tmp; v0(4) += tmp;
              break;
          default:
              break;
          }
      }
  }

  return 0;
}

  void
  ForceBeamColumn3d::compSectionDisplacements(Vector sectionCoords[],
					      Vector sectionDispls[]) const
  {
     // get basic displacements and increments
     static Vector ub(NEBD);
     ub = crdTransf->getBasicTrialDisp();    

     double L = crdTransf->getInitialLength();

     // get integration point positions and weights
     static double pts[maxNumSections];
     beamIntegr->getSectionLocations(numSections, L, pts);

     // setup Vandermode and CBDI influence matrices
     int i;
     double xi;

     // get CBDI influence matrix
     Matrix ls(numSections, numSections);
     getCBDIinfluenceMatrix(numSections, pts, L, ls);

     // get section curvatures
     Vector kappa_y(numSections);  // curvature
     Vector kappa_z(numSections);  // curvature
     static Vector vs;                // section deformations 

     for (i=0; i<numSections; i++) {
	 // THIS IS VERY INEFFICIENT ... CAN CHANGE IF RUNS TOO SLOW
	 int sectionKey1 = 0;
	 int sectionKey2 = 0;
	 const ID &code = sections[i]->getType();
	 int j;
	 for (j = 0; j < code.Size(); j++)
	 {
	     if (code(j) == SECTION_RESPONSE_MZ)
		 sectionKey1 = j;
	     if (code(j) == SECTION_RESPONSE_MY)
		 sectionKey2 = j;
	 }
	 if (sectionKey1 == 0) {
	   opserr << "FATAL ForceBeamColumn3d::compSectionResponse - section does not provide Mz response\n";
	   exit(-1);
	 }
	 if (sectionKey2 == 0) {
	   opserr << "FATAL ForceBeamColumn3d::compSectionResponse - section does not provide My response\n";
	   exit(-1);
	 }

	 // get section deformations
	 vs = sections[i]->getSectionDeformation();

	 kappa_z(i) = vs(sectionKey1);
	 kappa_y(i) = vs(sectionKey2); 
     }

     //cout << "kappa_y: " << kappa_y;   
     //cout << "kappa_z: " << kappa_z;   

     Vector v(numSections), w(numSections);
     static Vector xl(NDM), uxb(NDM);
     static Vector xg(NDM), uxg(NDM); 
     // double theta;                             // angle of twist of the sections

     // v = ls * kappa_z;  
     v.addMatrixVector (0.0, ls, kappa_z, 1.0);  
     // w = ls * kappa_y *  (-1);  
     w.addMatrixVector (0.0, ls, kappa_y, -1.0);

     for (i=0; i<numSections; i++)
     {
	xi = pts[i];

	xl(0) = xi * L;
	xl(1) = 0;
	xl(2) = 0;

	// get section global coordinates
	sectionCoords[i] = crdTransf->getPointGlobalCoordFromLocal(xl);

	// compute section displacements
	//theta  = xi * ub(5); // consider linear variation for angle of twist. CHANGE LATER!!!!!!!!!!
	uxb(0) = xi * ub(0); // consider linear variation for axial displacement. CHANGE LATER!!!!!!!!!!
	uxb(1) = v(i);
	uxb(2) = w(i);

	// get section displacements in global system 
	sectionDispls[i] = crdTransf->getPointGlobalDisplFromBasic(xi, uxb);
     }	       
    return;	       
  }

  void
  ForceBeamColumn3d::Print(OPS_Stream &s, int flag)
  {
    // flags with negative values are used by GSA
    if (flag == -1) { 
      int eleTag = this->getTag();
      s << "EL_BEAM\t" << eleTag << "\t";
      s << sections[0]->getTag() << "\t" << sections[numSections-1]->getTag(); 
      s  << "\t" << connectedExternalNodes(0) << "\t" << connectedExternalNodes(1);
      s << "\t0\t0.0000000\n";
    }  

    // flags with negative values are used by GSA  
    else if (flag < -1) {
      int eleTag = this->getTag();
      int counter = (flag +1) * -1;
      double P  = Secommit(0);
      double MZ1 = Secommit(1);
      double MZ2 = Secommit(2);
      double MY1 = Secommit(3);
      double MY2 = Secommit(4);
      double L = crdTransf->getInitialLength();
      double VY = (MZ1+MZ2)/L;
      theVector(1) =  VY;
      theVector(4) = -VY;
      double VZ = (MY1+MY2)/L;
      double T  = Secommit(5);

      double p0[5]; p0[0] = p0[1] = p0[2] = p0[3] = p0[4] = 0.0;
      if (numEleLoads > 0)
          this->computeReactions(p0);

      s << "FORCE\t" << eleTag << "\t" << counter << "\t0";
      s << "\t" << -P+p0[0] << "\t"  <<  VY+p0[1] << "\t"  << -VZ+p0[3]  << endln;
      s << "FORCE\t" << eleTag << "\t" << counter << "\t1";
      s << "\t"  << P  << ' '  << -VY+p0[2] << ' ' << VZ+p0[4] << endln;
      s << "MOMENT\t" << eleTag << "\t" << counter << "\t0";
      s << "\t" << -T << "\t"  << MY1 << "\t" << MZ1 << endln;
      s << "MOMENT\t" << eleTag << "\t" << counter << "\t1";
      s << "\t" << T << ' ' << MY2 << ' '  <<  MZ2 << endln;
    }

    // flag set to 2 used to print everything .. used for viewing data for UCSD renderer  
    else if (flag == 2) {
       static Vector xAxis(3);
       static Vector yAxis(3);
       static Vector zAxis(3);


       crdTransf->getLocalAxes(xAxis, yAxis, zAxis);

       s << "#ForceBeamColumn3D\n";
       s << "#LocalAxis " << xAxis(0) << " " << xAxis(1) << " " << xAxis(2);
       s << " " << yAxis(0) << " " << yAxis(1) << " " << yAxis(2) << " ";
       s << zAxis(0) << " " << zAxis(1) << " " << zAxis(2) << endln;

       const Vector &node1Crd = theNodes[0]->getCrds();
       const Vector &node2Crd = theNodes[1]->getCrds();	
       const Vector &node1Disp = theNodes[0]->getDisp();
       const Vector &node2Disp = theNodes[1]->getDisp();    

       s << "#NODE " << node1Crd(0) << " " << node1Crd(1) << " " << node1Crd(2)
	 << " " << node1Disp(0) << " " << node1Disp(1) << " " << node1Disp(2)
	 << " " << node1Disp(3) << " " << node1Disp(4) << " " << node1Disp(5) << endln;

       s << "#NODE " << node2Crd(0) << " " << node2Crd(1) << " " << node2Crd(2)
	 << " " << node2Disp(0) << " " << node2Disp(1) << " " << node2Disp(2)
	 << " " << node2Disp(3) << " " << node2Disp(4) << " " << node2Disp(5) << endln;

       double P  = Secommit(0);
       double MZ1 = Secommit(1);
       double MZ2 = Secommit(2);
       double MY1 = Secommit(3);
       double MY2 = Secommit(4);
       double L = crdTransf->getInitialLength();
       double VY = (MZ1+MZ2)/L;
       theVector(1) =  VY;
       theVector(4) = -VY;
       double VZ = (MY1+MY2)/L;
       double T  = Secommit(5);

       double p0[5]; p0[0] = p0[1] = p0[2] = p0[3] = p0[4] = 0.0;
       if (numEleLoads > 0)
           this->computeReactions(p0);

       s << "#END_FORCES " << -P+p0[0] << ' '  <<  VY+p0[1] << ' '  << -VZ+p0[3] << ' ' 
	 << -T << ' '  << MY1 << ' ' << MZ1 << endln;
       s << "#END_FORCES "  << P  << ' '  << -VY+p0[2] << ' ' << VZ+p0[4] << ' '  
	 << T << ' ' << MY2 << ' '  <<  MZ2 << endln;

       // plastic hinge rotation
       static Vector vp(6);
       static Matrix fe(6,6);
       this->getInitialFlexibility(fe);
       vp = crdTransf->getBasicTrialDisp();
       vp.addMatrixVector(1.0, fe, Se, -1.0);
       s << "#PLASTIC_HINGE_ROTATION " << vp[1] << " " << vp[2] << " " << vp[3] << " " << vp[4] 
	 << " " << 0.1*L << " " << 0.1*L << endln;

       // allocate array of vectors to store section coordinates and displacements
       static int maxNumSections = 0;
       static Vector *coords = 0;
       static Vector *displs = 0;
       if (maxNumSections < numSections) {
	 if (coords != 0) 
	   delete [] coords;
	 if (displs != 0)
	   delete [] displs;

	 coords = new Vector [numSections];
	 displs = new Vector [numSections];

	 if (!coords) {
	   opserr << "ForceBeamColumn3d::Print() -- failed to allocate coords array";   
	   exit(-1);
	 }

	 int i;
	 for (i = 0; i < numSections; i++)
	   coords[i] = Vector(NDM);

	 if (!displs) {
	   opserr << "ForceBeamColumn3d::Print() -- failed to allocate coords array";   
	   exit(-1);
	 }

	 for (i = 0; i < numSections; i++)
	   displs[i] = Vector(NDM);

	 maxNumSections = numSections;
       }

       // compute section location & displacements
       this->compSectionDisplacements(coords, displs);

       // spit out the section location & invoke print on the scetion
       for (int i=0; i<numSections; i++) {
	 s << "#SECTION " << (coords[i])(0) << " " << (coords[i])(1) << " " << (coords[i])(2);       
	 s << " " << (displs[i])(0) << " " << (displs[i])(1) << " " << (displs[i])(2) << endln;
	 sections[i]->Print(s, flag); 
       }
     }

    if (flag == OPS_PRINT_CURRENTSTATE) {
       s << "\nElement: " << this->getTag() << " Type: ForceBeamColumn3d ";
       s << "\tConnected Nodes: " << connectedExternalNodes ;
       s << "\tNumber of Sections: " << numSections;
       s << "\tMass density: " << rho << endln;
       beamIntegr->Print(s, flag);
       double P  = Secommit(0);
       double MZ1 = Secommit(1);
       double MZ2 = Secommit(2);
       double MY1 = Secommit(3);
       double MY2 = Secommit(4);
       double L = crdTransf->getInitialLength();
       double VY = (MZ1+MZ2)/L;
       theVector(1) =  VY;
       theVector(4) = -VY;
       double VZ = (MY1+MY2)/L;
       double T  = Secommit(5);

       double p0[5]; p0[0] = p0[1] = p0[2] = p0[3] = p0[4] = 0.0;
       if (numEleLoads > 0)
           this->computeReactions(p0);

       s << "\tEnd 1 Forces (P MZ VY MY VZ T): "
	 << -P+p0[0] << " " << MZ1 << " " <<  VY+p0[1] << " " 
	 << MY1 << " " << -VZ+p0[3] << " " << T << endln;
       s << "\tEnd 2 Forces (P MZ VY MY VZ T): "
	 << P        << " " << MZ2 << " " << -VY+p0[2] << " " 
	 << MY2 << " " <<  VZ+p0[4] << " " << -T << endln;

       for (int i = 0; i < numSections; i++) {
	 //opserr << "Section Type: " << theSections[i]->getClassTag() << endln;
	 sections[i]->Print(s,flag);
       }
    }

	 if (flag == OPS_PRINT_PRINTMODEL_JSON) {
		 s << "\t\t\t{";
		 s << "\"name\": " << this->getTag() << ", ";
		 s << "\"type\": \"ForceBeamColumn3d\", ";
		 s << "\"nodes\": [" << connectedExternalNodes(0) << ", " << connectedExternalNodes(1) << "], ";
		 s << "\"sections\": [";
		 for (int i = 0; i < numSections - 1; i++)
			 s << "\"" << sections[i]->getTag() << "\", ";
		 s << "\"" << sections[numSections - 1]->getTag() << "\"], ";
		 s << "\"integration\": ";
		 beamIntegr->Print(s, flag);
		 s << ", \"massperlength\": " << rho << ", ";
		 s << "\"crdTransformation\": \"" << crdTransf->getTag() << "\"}";
	 }
  }

  OPS_Stream &operator<<(OPS_Stream &s, ForceBeamColumn3d &E)
  {
    E.Print(s);
    return s;
  }

  int
  ForceBeamColumn3d::displaySelf(Renderer &theViewer, int displayMode, float fact, const char** displayModes, int numModes)
  {
    static Vector v1(3);
    static Vector v2(3);

    theNodes[0]->getDisplayCrds(v1, fact, displayMode);
    theNodes[1]->getDisplayCrds(v2, fact, displayMode);

    return theViewer.drawLine(v1, v2, 1.0, 1.0, this->getTag());
  }

  Response*
  ForceBeamColumn3d::setResponse(const char **argv, int argc, OPS_Stream &output)
  {
    Response *theResponse = 0;
    
    output.tag("ElementOutput");
    output.attr("eleType","ForceBeamColumn3d");
    output.attr("eleTag",this->getTag());
    output.attr("node1",connectedExternalNodes[0]);
    output.attr("node2",connectedExternalNodes[1]);

    //
    // we compare argv[0] for known response types 
    //

    // global force - 
    if (strcmp(argv[0],"forces") == 0 || strcmp(argv[0],"force") == 0
	|| strcmp(argv[0],"globalForce") == 0 || strcmp(argv[0],"globalForces") == 0) {

      output.tag("ResponseType","Px_1");
      output.tag("ResponseType","Py_1");
      output.tag("ResponseType","Pz_1");
      output.tag("ResponseType","Mx_1");
      output.tag("ResponseType","My_1");
      output.tag("ResponseType","Mz_1");
      output.tag("ResponseType","Px_2");
      output.tag("ResponseType","Py_2");
      output.tag("ResponseType","Pz_2");
      output.tag("ResponseType","Mx_2");
      output.tag("ResponseType","My_2");
      output.tag("ResponseType","Mz_2");


      theResponse = new ElementResponse(this, 1, theVector);

    // local force -
    }  else if (strcmp(argv[0],"localForce") == 0 || strcmp(argv[0],"localForces") == 0) {

      output.tag("ResponseType","N_1");
      output.tag("ResponseType","Vy_1");
      output.tag("ResponseType","Vz_1");
      output.tag("ResponseType","T_1");
      output.tag("ResponseType","My_1");
      output.tag("ResponseType","Mz_1");
      output.tag("ResponseType","N_2");
      output.tag("ResponseType","Vy_2");
      output.tag("ResponseType","Vz_2");
      output.tag("ResponseType","T_2");
      output.tag("ResponseType","My_2");
      output.tag("ResponseType","Mz_2");
      
      theResponse = new ElementResponse(this, 2, theVector);

    // basic force -
    } else if (strcmp(argv[0],"basicForce") == 0 || strcmp(argv[0],"basicForces") == 0) {
      
      output.tag("ResponseType","N");
      output.tag("ResponseType","Mz_1");
      output.tag("ResponseType","Mz_2");
      output.tag("ResponseType","My_1");
      output.tag("ResponseType","My_2");
      output.tag("ResponseType","T");            
      
      theResponse =  new ElementResponse(this, 7, Vector(6));

    // basic stiffness -
    } else if (strcmp(argv[0],"basicStiffness") == 0) {

      output.tag("ResponseType","N");
      output.tag("ResponseType","Mz_1");
      output.tag("ResponseType","Mz_2");
      output.tag("ResponseType","My_1");
      output.tag("ResponseType","My_2");
      output.tag("ResponseType","T");                  
      
      theResponse =  new ElementResponse(this, 19, Matrix(6,6));
      
    //global damping force -
    } else if (theDamping && (strcmp(argv[0],"globalDampingForce") == 0 || strcmp(argv[0],"globalDampingForces") == 0)) {

      output.tag("ResponseType","Px_1");
      output.tag("ResponseType","Py_1");
      output.tag("ResponseType","Pz_1");
      output.tag("ResponseType","Mx_1");
      output.tag("ResponseType","My_1");
      output.tag("ResponseType","Mz_1");
      output.tag("ResponseType","Px_2");
      output.tag("ResponseType","Py_2");
      output.tag("ResponseType","Pz_2");
      output.tag("ResponseType","Mx_2");
      output.tag("ResponseType","My_2");
      output.tag("ResponseType","Mz_2");


      theResponse = new ElementResponse(this, 21, theVector);

    // local damping force -
    } else if (theDamping && (strcmp(argv[0],"localDampingForce") == 0 || strcmp(argv[0],"localDampingForces") == 0)) {

      output.tag("ResponseType","N_1");
      output.tag("ResponseType","Vy_1");
      output.tag("ResponseType","Vz_1");
      output.tag("ResponseType","T_1");
      output.tag("ResponseType","My_1");
      output.tag("ResponseType","Mz_1");
      output.tag("ResponseType","N_2");
      output.tag("ResponseType","Vy_2");
      output.tag("ResponseType","Vz_2");
      output.tag("ResponseType","T_2");
      output.tag("ResponseType","My_2");
      output.tag("ResponseType","Mz_2");
      
      theResponse = new ElementResponse(this, 22, theVector);

    } else if (theDamping && (strcmp(argv[0],"basicDampingForce") == 0 || strcmp(argv[0],"basicDampingForces") == 0)) {

      output.tag("ResponseType","N");
      output.tag("ResponseType","Mz_1");
      output.tag("ResponseType","Mz_2");
      output.tag("ResponseType","My_1");
      output.tag("ResponseType","My_2");
      output.tag("ResponseType","T");
    
      theResponse = new ElementResponse(this, 23, Vector(6));
      
    // chord rotation -
    }  else if (strcmp(argv[0],"chordRotation") == 0 || strcmp(argv[0],"chordDeformation") == 0 
		|| strcmp(argv[0],"basicDeformation") == 0) {
      
      output.tag("ResponseType","eps");
      output.tag("ResponseType","thetaZ_1");
      output.tag("ResponseType","thetaZ_2");
      output.tag("ResponseType","thetaY_1");
      output.tag("ResponseType","thetaY_2");
      output.tag("ResponseType","thetaX");
      
      theResponse = new ElementResponse(this, 3, Vector(6));
      
      // plastic rotation -
    } else if (strcmp(argv[0],"plasticRotation") == 0 || strcmp(argv[0],"plasticDeformation") == 0) {
      
      output.tag("ResponseType","epsP");
      output.tag("ResponseType","thetaZP_1");
      output.tag("ResponseType","thetaZP_2");
      output.tag("ResponseType","thetaYP_1");
      output.tag("ResponseType","thetaYP_2");
      output.tag("ResponseType","thetaXP");
      
      theResponse = new ElementResponse(this, 4, Vector(6));
      
      // point of inflection
    } else if (strcmp(argv[0],"inflectionPoint") == 0) {
      theResponse = new ElementResponse(this, 5, Vector(2));
      
      // tangent drift
    } else if (strcmp(argv[0],"tangentDrift") == 0) {
      theResponse = new ElementResponse(this, 6, Vector(4));
      
    } else if (strcmp(argv[0],"getRemCriteria1") == 0) {
      theResponse = new ElementResponse(this, 77, Vector(2));

    } else if (strcmp(argv[0],"getRemCriteria2") == 0) {
      theResponse = new ElementResponse(this, 8, Vector(2), ID(6));

    } else if (strcmp(argv[0],"RayleighForces") == 0 || 
	       strcmp(argv[0],"rayleighForces") == 0) {

      theResponse = new ElementResponse(this, 12, theVector);

    } else if (strcmp(argv[0],"sections") ==0) { 
      CompositeResponse *theCResponse = new CompositeResponse();
      int numResponse = 0;
      double xi[maxNumSections];
      double L = crdTransf->getInitialLength();
      beamIntegr->getSectionLocations(numSections, L, xi);
      
      for (int i=0; i<numSections; i++) {
	
	output.tag("GaussPointOutput");
	output.attr("number",i+1);
	output.attr("eta",xi[i]*L);
	
	Response *theSectionResponse = sections[i]->setResponse(&argv[1], argc-1, output);
	
	if (theSectionResponse != 0) {
	  numResponse = theCResponse->addResponse(theSectionResponse);
	}
      }
      
      if (numResponse == 0) // no valid responses found
	delete theCResponse;
	  else
	    theResponse = theCResponse;
    }

    else if (strcmp(argv[0],"integrationPoints") == 0)
      theResponse = new ElementResponse(this, 10, Vector(numSections));

    else if (strcmp(argv[0],"integrationWeights") == 0)
      theResponse = new ElementResponse(this, 11, Vector(numSections));

    else if (strcmp(argv[0],"sectionTags") == 0)
      theResponse = new ElementResponse(this, 110, ID(numSections));  
    
    else if (strcmp(argv[0],"sectionDisplacements") == 0) {
      if (argc > 1 && strcmp(argv[1],"local") == 0)
	theResponse = new ElementResponse(this, 1111, Matrix(numSections,3));
      else
	theResponse = new ElementResponse(this, 111, Matrix(numSections,3));
    }
    
    else if (strcmp(argv[0],"cbdiDisplacements") == 0)
      theResponse = new ElementResponse(this, 112, Matrix(1,3));

    // section response -
    else if (strcmp(argv[0],"sectionX") == 0) {
      if (argc > 2) {
	float sectionLoc = atof(argv[1]);
	
	double xi[maxNumSections];
	double L = crdTransf->getInitialLength();
	beamIntegr->getSectionLocations(numSections, L, xi);
	
	sectionLoc /= L;
	
	float minDistance = fabs(xi[0]-sectionLoc);
	int sectionNum = 0;
	for (int i = 1; i < numSections; i++) {
	  if (fabs(xi[i]-sectionLoc) < minDistance) {
	    minDistance = fabs(xi[i]-sectionLoc);
	    sectionNum = i;
	  }
	}
	
	output.tag("GaussPointOutput");
	output.attr("number",sectionNum+1);
	output.attr("eta",xi[sectionNum]*L);
	
	theResponse = sections[sectionNum]->setResponse(&argv[2], argc-2, output);
      }
    }

    else if (strcmp(argv[0],"section") == 0) { 

      if (argc > 1) {

	int sectionNum = atoi(argv[1]);
	
	if (sectionNum > 0 && sectionNum <= numSections && argc > 2) {
	  double xi[maxNumSections];
	  double L = crdTransf->getInitialLength();
	  beamIntegr->getSectionLocations(numSections, L, xi);
	  
	  output.tag("GaussPointOutput");
	  output.attr("number",sectionNum);
	  output.attr("eta",2.0*xi[sectionNum-1]-1.0);

	  if (strcmp(argv[2],"dsdh") != 0) {
	    theResponse = sections[sectionNum-1]->setResponse(&argv[2], argc-2, output);
	  } else {
	    int order = sections[sectionNum-1]->getOrder();
	    theResponse = new ElementResponse(this, 76, Vector(order));
	    Information &info = theResponse->getInformation();
	    info.theInt = sectionNum;
	  }
	  
	  output.endTag();
	  
	} else if (sectionNum == 0) { // argv[1] was not an int, we want all sections, 

	  CompositeResponse *theCResponse = new CompositeResponse();
	  int numResponse = 0;
	  double xi[maxNumSections];
	  double L = crdTransf->getInitialLength();
	  beamIntegr->getSectionLocations(numSections, L, xi);
	  
	  for (int i=0; i<numSections; i++) {
	    
	    output.tag("GaussPointOutput");
	    output.attr("number",i+1);
	    output.attr("eta",xi[i]*L);
	    
	    Response *theSectionResponse = sections[i]->setResponse(&argv[1], argc-1, output);
	    
	    if (theSectionResponse != 0) {
	      numResponse = theCResponse->addResponse(theSectionResponse);
	    }
	  }
	  
	  if (numResponse == 0) // no valid responses found
	    delete theCResponse;
	  else
	    theResponse = theCResponse;
	}
      }
    }
	//by SAJalali
	else if (strcmp(argv[0], "energy") == 0)
	{
		theResponse = new ElementResponse(this, 10, 0.0);
	}

    if (theResponse == 0) {
      theResponse = crdTransf->setResponse(argv, argc, output);
    }
    
    output.endTag();

    return theResponse;
}

int 
ForceBeamColumn3d::getResponse(int responseID, Information &eleInfo)
{
  static Vector vp(6);
  static Matrix fe(6,6);

  if (responseID == 1)
    return eleInfo.setVector(this->getResistingForce());
  
  else if (responseID == 2) {
    double p0[5]; p0[0] = p0[1] = p0[2] = p0[3] = p0[4] = 0.0;
    if (numEleLoads > 0)
      this->computeReactions(p0);
    // Axial
    double N = Se(0);
    theVector(6) =  N;
    theVector(0) = -N+p0[0];
    
    // Torsion
    double T = Se(5);
    theVector(9) =  T;
    theVector(3) = -T;
    
    // Moments about z and shears along y
    double M1 = Se(1);
    double M2 = Se(2);
    theVector(5)  = M1;
    theVector(11) = M2;
    double L = crdTransf->getInitialLength();
    double V = (M1+M2)/L;
    theVector(1) =  V+p0[1];
    theVector(7) = -V+p0[2];
    
    // Moments about y and shears along z
    M1 = Se(3);
    M2 = Se(4);
    theVector(4)  = M1;
    theVector(10) = M2;
    V = (M1+M2)/L;
    theVector(2) = -V+p0[3];
    theVector(8) =  V+p0[4];
      
    return eleInfo.setVector(theVector);

  }
      
  else if (responseID == 21)
    return eleInfo.setVector(this->getDampingForce());

  else if (responseID == 22) {
    Vector Sd(NEBD);
    Sd = theDamping->getDampingForce();
    // Axial
    double N = Sd(0);
    theVector(6) =  N;
    theVector(0) = -N;
    
    // Torsion
    double T = Sd(5);
    theVector(9) =  T;
    theVector(3) = -T;
    
    // Moments about z and shears along y
    double M1 = Sd(1);
    double M2 = Sd(2);
    theVector(5)  = M1;
    theVector(11) = M2;
    double L = crdTransf->getInitialLength();
    double V = (M1+M2)/L;
    theVector(1) =  V;
    theVector(7) = -V;
    
    // Moments about y and shears along z
    M1 = Sd(3);
    M2 = Sd(4);
    theVector(4)  = M1;
    theVector(10) = M2;
    V = (M1+M2)/L;
    theVector(2) = -V;
    theVector(8) =  V;
      
    return eleInfo.setVector(theVector);
  }

  else if (responseID == 23)
    return eleInfo.setVector(theDamping->getDampingForce());
      
  // Chord rotation
  else if (responseID == 3) {
    vp = crdTransf->getBasicTrialDisp();
    return eleInfo.setVector(vp);
  }

  else if (responseID == 7)
    return eleInfo.setVector(Se);

  else if (responseID == 19)
    return eleInfo.setMatrix(kv);
  
  // Plastic rotation
  else if (responseID == 4) {
    this->getInitialFlexibility(fe);
    vp = crdTransf->getBasicTrialDisp();
    vp.addMatrixVector(1.0, fe, Se, -1.0);
    Vector v0(6);
    this->getInitialDeformations(v0);
    vp.addVector(1.0, v0, -1.0);
    return eleInfo.setVector(vp);
  }

  else if (responseID == 10) {
    double L = crdTransf->getInitialLength();
    double pts[maxNumSections];
    beamIntegr->getSectionLocations(numSections, L, pts);
    Vector locs(numSections);
    for (int i = 0; i < numSections; i++)
      locs(i) = pts[i]*L;
    return eleInfo.setVector(locs);
  }

  else if (responseID == 11) {
    double L = crdTransf->getInitialLength();
    double wts[maxNumSections];
    beamIntegr->getSectionWeights(numSections, L, wts);
    Vector weights(numSections);
    for (int i = 0; i < numSections; i++)
      weights(i) = wts[i]*L;
    return eleInfo.setVector(weights);
  }

  else if (responseID == 110) {
    ID tags(numSections);
    for (int i = 0; i < numSections; i++)
      tags(i) = sections[i]->getTag();
    return eleInfo.setID(tags);
  }
  
  else if (responseID == 111 || responseID == 1111) {
    double L = crdTransf->getInitialLength();
    double pts[maxNumSections];
    beamIntegr->getSectionLocations(numSections, L, pts);
    // CBDI influence matrix
    Matrix ls(numSections, numSections);
    getCBDIinfluenceMatrix(numSections, pts, L, ls);
    // Curvature vector
    Vector kappaz(numSections); // about section z
    Vector kappay(numSections); // about section y
    for (int i = 0; i < numSections; i++) {
      const ID &code = sections[i]->getType();
      const Vector &e = sections[i]->getSectionDeformation();
      int order = sections[i]->getOrder();
      for (int j = 0; j < order; j++) {
	if (code(j) == SECTION_RESPONSE_MZ)
	  kappaz(i) += e(j);
	if (code(j) == SECTION_RESPONSE_MY)
	  kappay(i) += e(j);
      }
    }
    // Displacement vector
    Vector dispsy(numSections); // along local y
    Vector dispsz(numSections); // along local z    
    dispsy.addMatrixVector(0.0, ls, kappaz,  1.0);
    dispsz.addMatrixVector(0.0, ls, kappay, -1.0);    
    beamIntegr->getSectionLocations(numSections, L, pts);
    static Vector uxb(3);
    static Vector uxg(3);
    Matrix disps(numSections,3);
    vp = crdTransf->getBasicTrialDisp();
    for (int i = 0; i < numSections; i++) {
      uxb(0) = pts[i]*vp(0); // linear shape function
      uxb(1) = dispsy(i);
      uxb(2) = dispsz(i);
      if (responseID == 111)
	uxg = crdTransf->getPointGlobalDisplFromBasic(pts[i],uxb);
      else
	uxg = crdTransf->getPointLocalDisplFromBasic(pts[i],uxb);
      disps(i,0) = uxg(0);
      disps(i,1) = uxg(1);
      disps(i,2) = uxg(2);            
    }
    return eleInfo.setMatrix(disps);
  }

  else if (responseID == 112) {
    double L = crdTransf->getInitialLength();
    double ipts[maxNumSections];
    beamIntegr->getSectionLocations(numSections, L, ipts);
    // CBDI influence matrix
    double pts[1];
    pts[0] = eleInfo.theDouble;
    Matrix ls(1, numSections);
    getCBDIinfluenceMatrix(1, pts, numSections, ipts, L, ls);
    // Curvature vector
    Vector kappaz(numSections); // about section z
    Vector kappay(numSections); // about section y    
    for (int i = 0; i < numSections; i++) {
      const ID &code = sections[i]->getType();
      const Vector &e = sections[i]->getSectionDeformation();
      int order = sections[i]->getOrder();
      for (int j = 0; j < order; j++) {
	if (code(j) == SECTION_RESPONSE_MZ)
	  kappaz(i) += e(j);
	if (code(j) == SECTION_RESPONSE_MY)
	  kappay(i) += e(j);
      }
    }
    // Displacement vector
    Vector dispsy(1); // along local y
    Vector dispsz(1); // along local z    
    dispsy.addMatrixVector(0.0, ls, kappaz,  1.0);
    dispsz.addMatrixVector(0.0, ls, kappay, -1.0);
    static Vector uxb(3);
    static Vector uxg(3);
    Matrix disps(1,3);
    vp = crdTransf->getBasicTrialDisp();
    uxb(0) = pts[0]*vp(0); // linear shape function
    uxb(1) = dispsy(0);
    uxb(2) = dispsz(0);      
    uxg = crdTransf->getPointGlobalDisplFromBasic(pts[0],uxb);
    disps(0,0) = uxg(0);
    disps(0,1) = uxg(1);
    disps(0,2) = uxg(2);            

    return eleInfo.setMatrix(disps);
  }

  else if (responseID == 12)
    return eleInfo.setVector(this->getRayleighDampingForces());

  // Point of inflection
  else if (responseID == 5) {
    static Vector LI(2);
    LI(0) = 0.0;
    LI(1) = 0.0;

    double L = crdTransf->getInitialLength();

    if (fabs(Se(1)+Se(2)) > DBL_EPSILON)
      LI(0) = Se(1)/(Se(1)+Se(2))*L;

    if (fabs(Se(3)+Se(4)) > DBL_EPSILON)
      LI(1) = Se(3)/(Se(3)+Se(4))*L;

    return eleInfo.setVector(LI);
  }

  // Tangent drift
  else if (responseID == 6) {
    double d2z = 0.0;
    double d2y = 0.0;
    double d3z = 0.0;
    double d3y = 0.0;

    double L = crdTransf->getInitialLength();

    double wts[maxNumSections];
    beamIntegr->getSectionWeights(numSections, L, wts);

    double pts[maxNumSections];
    beamIntegr->getSectionLocations(numSections, L, pts);

    // Location of inflection point from node I
    double LIz = 0.0;
    if (fabs(Se(1)+Se(2)) > DBL_EPSILON)
      LIz = Se(1)/(Se(1)+Se(2))*L;

    double LIy = 0.0;
    if (fabs(Se(3)+Se(4)) > DBL_EPSILON)
      LIy = Se(3)/(Se(3)+Se(4))*L;

    int i;
    for (i = 0; i < numSections; i++) {
      double x = pts[i]*L;
      const ID &type = sections[i]->getType();
      int order = sections[i]->getOrder();
      double kappa = 0.0;
      if (x < LIz) {
	for (int j = 0; j < order; j++)
	  if (type(j) == SECTION_RESPONSE_MZ)
	    kappa += vs[i](j);
	double b = -LIz+x;
	d2z += (wts[i]*L)*kappa*b;
      }
      kappa = 0.0;
      if (x < LIy) {
	for (int j = 0; j < order; j++)
	  if (type(j) == SECTION_RESPONSE_MY)
	    kappa += vs[i](j);
	double b = -LIy+x;
	d2y += (wts[i]*L)*kappa*b;
      }
    }

    for (i = numSections-1; i >= 0; i--) {
      double x = pts[i]*L;
      const ID &type = sections[i]->getType();
      int order = sections[i]->getOrder();
      double kappa = 0.0;
      if (x > LIz) {
	for (int j = 0; j < order; j++)
	  if (type(j) == SECTION_RESPONSE_MZ)
	    kappa += vs[i](j);
	double b = x-LIz;
	d3z += (wts[i]*L)*kappa*b;
      }
      kappa = 0.0;
      if (x > LIy) {
	for (int j = 0; j < order; j++)
	  if (type(j) == SECTION_RESPONSE_MY)
	    kappa += vs[i](j);
	double b = x-LIy;
	d3y += (wts[i]*L)*kappa*b;
      }
    }

    static Vector d(4);
    d(0) = d2z;
    d(1) = d3z;
    d(2) = d2y;
    d(3) = d3y;

    return eleInfo.setVector(d);

  } else if (responseID == 77) { // Why is this here?
    return -1;
  } else if (responseID == 8) {

    ID *eleInfoID = eleInfo.theID;

    int compID = (*eleInfoID)(0);
    int critID = (*eleInfoID)(1);
    int nTagbotn11 = (*eleInfoID)(2);
    int nTagmidn11 = (*eleInfoID)(3);
    int nTagtopn11 = (*eleInfoID)(4);
    int globgrav11 = (*eleInfoID)(5);

    const char* filenamewall = eleInfo.theString;

    // int returns
    double value = 0.0;
    double checkvalue1 = 0.0;

    if (critID == 7) {
      Domain *theDomain = this->getDomain();

      double oofwallresp;
      // determine the in plane horizontal deformation axis
      // and the out of plane horizontal deformation axis
      Node *theNode1a = theDomain->getNode(nTagbotn11);
      Node *theNode3a = theDomain->getNode(nTagtopn11); 
      const Vector &crdIa1 = theNode1a->getCrds();
      const Vector &crdJa1 = theNode3a->getCrds();
      int indwdir1;
      int indwdir2;
      if (globgrav11==1) {
	indwdir1=1;
	indwdir2=2;
      }
      else if (globgrav11==2) {
	indwdir1=0;
	indwdir2=2;		  
      }
      else if (globgrav11==3) {
	indwdir1=0;
	indwdir2=1;		  
      }

      double dir1a1=crdJa1(indwdir1)-crdIa1(indwdir1);
      double dir2a1=crdJa1(indwdir2)-crdIa1(indwdir2);
      double dirsumdum=sqrt(dir1a1*dir1a1+dir2a1*dir2a1);
      double dir1inp=dir1a1/dirsumdum;		
      double dir2inp=dir2a1/dirsumdum;
      
      double dir1oop=-dir2inp;
      double dir2oop=dir1inp;
      
      Node *theNode1 = theDomain->getNode(nTagbotn11);
      const Vector &theResponsewall = theNode1->getTrialDisp();
      double valbotinfn=theResponsewall(indwdir1)*dir1inp+theResponsewall(indwdir2)*dir2inp;
      double valbotoutfn=theResponsewall(indwdir1)*dir1oop+theResponsewall(indwdir2)*dir2oop;
      
      Node *theNode2 = theDomain->getNode(nTagmidn11);
      const Vector &theResponsewall2 = theNode2->getTrialDisp();
      double valmidinfn=theResponsewall2(indwdir1)*dir1inp+theResponsewall2(indwdir2)*dir2inp;
      double valmidoutfn=theResponsewall2(indwdir1)*dir1oop+theResponsewall2(indwdir2)*dir2oop;
      
      Node *theNode3 = theDomain->getNode(nTagtopn11);
      const Vector &theResponsewall3 = theNode3->getTrialDisp();
      double valtopinfn=theResponsewall3(indwdir1)*dir1inp+theResponsewall3(indwdir2)*dir2inp;
      double valtopoutfn=theResponsewall3(indwdir1)*dir1oop+theResponsewall3(indwdir2)*dir2oop;
      
      value = sqrt(pow((valtopinfn-valbotinfn),2.0));
      double valoutchck=valmidoutfn-(valtopoutfn+valbotoutfn)/2.0;
      oofwallresp = sqrt(pow(valoutchck,2.0));
      //        
      double outplanevaldat; // variable for input value
      double inplanevaldat; // variable for input value
      double outplanevaldat1; // variable for input value
      double inplanevaldat1; // variable for input value
      ifstream indata;

      if (filenamewall!=NULL) {
	//		
	indata.open(filenamewall); // opens the file
	if(!indata) { // file couldn't be opened
		opserr << "ForceBeamColumn3d::getResponse()"
			<< " file for infill wall (" << filenamewall << " could not be opened" << endln;
	  return -1;
	}
	checkvalue1=0.0;
	int counterdum=0;
	while ( !indata.eof() ) { // keep reading until end-of-file
	  counterdum=counterdum+1;
	  indata >> outplanevaldat >> inplanevaldat ; // sets EOF flag if no value found
	  if (counterdum!=1)  {
	    if (oofwallresp >= outplanevaldat1 && oofwallresp <= outplanevaldat)  {
	      checkvalue1= inplanevaldat1+(oofwallresp-outplanevaldat1)/(outplanevaldat-outplanevaldat1)*(inplanevaldat-inplanevaldat1);
	      break; 
	    }
	  }
	  indata >> outplanevaldat1 >> inplanevaldat1;
	  if (oofwallresp >= outplanevaldat && oofwallresp <= outplanevaldat1)  {
	    checkvalue1= inplanevaldat+(oofwallresp-outplanevaldat)/(outplanevaldat1-outplanevaldat)*(inplanevaldat1-inplanevaldat);
	    break;
	  }
	}
	indata.close();
      }

      static Vector result8(2);
      result8(0) = value;
      result8(1) = checkvalue1;      
      
      return eleInfo.setVector(result8);
    }

    return -1;
  }
  //by SAJalali
  else if (responseID == 10) {
	  double xi[maxNumSections];
	  double L = crdTransf->getInitialLength();
	  beamIntegr->getSectionWeights(numSections, L, xi);
	  double energy = 0;
	  for (int i = 0; i < numSections; i++) {
		  energy += sections[i]->getEnergy()*xi[i] * L;
	  }
	  return eleInfo.setDouble(energy);
  }

  else
    return -1;
}

int 
ForceBeamColumn3d::getResponseSensitivity(int responseID, int gradNumber,
					  Information &eleInfo)
{
  // Basic deformation sensitivity
  if (responseID == 3) {  
    const Vector &dvdh = crdTransf->getBasicDisplSensitivity(gradNumber);
    return eleInfo.setVector(dvdh);
  }

  // Basic force sensitivity
  else if (responseID == 7) {
    static Vector dqdh(6);

    const Vector &dvdh = crdTransf->getBasicDisplSensitivity(gradNumber);

    dqdh.addMatrixVector(0.0, kv, dvdh, 1.0);

    dqdh.addVector(1.0, this->computedqdh(gradNumber), 1.0);
    //opserr << "FBC2d: " << gradNumber;

    return eleInfo.setVector(dqdh);
  }

  // dsdh
  else if (responseID == 76) {

    int sectionNum = eleInfo.theInt;
    int order = sections[sectionNum-1]->getOrder();

    Vector dsdh(order);
    dsdh.Zero();

    if (numEleLoads > 0) {
      this->computeSectionForceSensitivity(dsdh, sectionNum-1, gradNumber);
    }
    //opserr << "FBC3d::getRespSens dspdh: " << dsdh;
    static Vector dqdh(6);

    const Vector &dvdh = crdTransf->getBasicDisplSensitivity(gradNumber);

    dqdh.addMatrixVector(0.0, kv, dvdh, 1.0);

    dqdh.addVector(1.0, this->computedqdh(gradNumber), 1.0);

    //opserr << "FBC3d::getRespSens dqdh: " << dqdh;
 
    double L = crdTransf->getInitialLength();
    double oneOverL  = 1.0/L;  
    double pts[maxNumSections];
    beamIntegr->getSectionLocations(numSections, L, pts);
    
    const ID &code = sections[sectionNum-1]->getType();
      
    double xL  = pts[sectionNum-1];
    double xL1 = xL-1.0;
    
    for (int ii = 0; ii < order; ii++) {
      switch(code(ii)) {
      case SECTION_RESPONSE_P:
	dsdh(ii) += dqdh(0);
	break;
      case SECTION_RESPONSE_MZ:
	dsdh(ii) +=  xL1*dqdh(1) + xL*dqdh(2);
	break;
      case SECTION_RESPONSE_VY:
	dsdh(ii) += oneOverL*(dqdh(1)+dqdh(2));
	break;
      case SECTION_RESPONSE_MY:
	dsdh(ii) +=  xL1*dqdh(3) + xL*dqdh(4);
	break;
      case SECTION_RESPONSE_VZ:
	dsdh(ii) += oneOverL*(dqdh(3)+dqdh(4));
	break;
      case SECTION_RESPONSE_T:
	dsdh(ii) += dqdh(5);
	break;
      default:
	dsdh(ii) += 0.0;
	break;
      }
    }
    
    double dLdh = crdTransf->getdLdh();
    double d1oLdh = crdTransf->getd1overLdh();
    
    double dptsdh[maxNumSections];
    beamIntegr->getLocationsDeriv(numSections, L, dLdh, dptsdh);
    double dxLdh = dptsdh[sectionNum-1];// - xL/L*dLdh;

    for (int j = 0; j < order; j++) {
      switch (code(j)) {
      case SECTION_RESPONSE_MZ:
	dsdh(j) += dxLdh*(Se(1)+Se(2));
	//dsdh(j) -= dLdh*xL/L*(Se(1)+Se(2));
	break;
      case SECTION_RESPONSE_VY:
	dsdh(j) += d1oLdh*(Se(1)+Se(2));
	break;
      case SECTION_RESPONSE_MY:
	dsdh(j) += dxLdh*(Se(3)+Se(4));
	break;
      case SECTION_RESPONSE_VZ:
	dsdh(j) += d1oLdh*(Se(3)+Se(4));
	break;
      default:
	break;
      }
    }

    /*
    opserr << "FBC3d::getRespSens dsdh=b*dqdh+dspdh: " << dsdh;

    dsdh.Zero();
    if (numEleLoads > 0) {
      this->computeSectionForceSensitivity(dsdh, sectionNum-1, gradNumber);
    }
    const Matrix &ks = sections[sectionNum-1]->getSectionTangent();
    const Vector &dedh  = sections[sectionNum-1]->getSectionDeformationSensitivity(gradNumber);
    dsdh.addMatrixVector(1.0, ks, dedh, 1.0);
    dsdh.addVector(1.0, sections[sectionNum-1]->getStressResultantSensitivity(gradNumber, true), 1.0);

    opserr << "FBC3d::getRespSens dsdh=b*dqdh+dspdh: " << dsdh;
    */

    return eleInfo.setVector(dsdh);
  }

  // Plastic deformation sensitivity
  else if (responseID == 4) {
    static Vector dvpdh(6);

    const Vector &dvdh = crdTransf->getBasicDisplSensitivity(gradNumber);

    dvpdh = dvdh;
    //opserr << dvpdh;

    static Matrix fe(6,6);
    this->getInitialFlexibility(fe);

    const Vector &dqdh = this->computedqdh(gradNumber);

    dvpdh.addMatrixVector(1.0, fe, dqdh, -1.0);
    //opserr << dvpdh;

    static Matrix fek(6,6);
    fek.addMatrixProduct(0.0, fe, kv, 1.0);

    dvpdh.addMatrixVector(1.0, fek, dvdh, -1.0);
    //opserr << dvpdh;

    const Matrix &dfedh = this->computedfedh(gradNumber);

    dvpdh.addMatrixVector(1.0, dfedh, Se, -1.0);
    //opserr << dvpdh << endln;
    //opserr << dfedh << endln;

    //opserr << dqdh + kv*dvdh << endln;

    return eleInfo.setVector(dvpdh);
  }

  else
    return -1;
}

int
ForceBeamColumn3d::setParameter(const char **argv, int argc, Parameter &param)
{
  if (argc < 1)
    return -1;

  int result = -1;

  // If the parameter belongs to the element itself
  if (strcmp(argv[0],"rho") == 0) {
    param.setValue(rho);
    return param.addObject(1, this);
  }

  // section response -
  if (strstr(argv[0],"sectionX") != 0) {
    if (argc > 2) {
      float sectionLoc = atof(argv[1]);

      double xi[maxNumSections];
      double L = crdTransf->getInitialLength();
      beamIntegr->getSectionLocations(numSections, L, xi);
      
      sectionLoc /= L;

      float minDistance = fabs(xi[0]-sectionLoc);
      int sectionNum = 0;
      for (int i = 1; i < numSections; i++) {
	if (fabs(xi[i]-sectionLoc) < minDistance) {
	  minDistance = fabs(xi[i]-sectionLoc);
	  sectionNum = i;
	}
      }

      return sections[sectionNum]->setParameter(&argv[2], argc-2, param);
    }
  }

  // If the parameter belongs to a particular section or lower
  if (strstr(argv[0],"section") != 0) {
    
    if (argc < 3)
      return -1;
    
    // Get section number
    int sectionNum = atoi(argv[1]);
   
    if (sectionNum > 0 && sectionNum <= numSections)
      return sections[sectionNum-1]->setParameter(&argv[2], argc-2, param);

    else
      return -1;
  }
  
  // If the parameter belongs to all sections or lower
  if (strstr(argv[0],"allSections") != 0) {
    
    if (argc < 2)
      return -1;
    
    int ok;
    for (int i = 0; i < numSections; i++) {
      ok = sections[i]->setParameter(&argv[1], argc-1, param);
      if (ok != -1)
	result = ok;
    }

    return result;
  }

  if (strstr(argv[0],"integration") != 0) {
    
    if (argc < 2)
      return -1;

    return beamIntegr->setParameter(&argv[1], argc-1, param);
  }

  // Default, send to everything
  int ok;

  for (int i = 0; i < numSections; i++) {
    ok = sections[i]->setParameter(argv, argc, param);
    if (ok != -1)
      result = ok;
  }

  ok = beamIntegr->setParameter(argv, argc, param);
  if (ok != -1)
    result = ok;

  return result;
}

int
ForceBeamColumn3d::updateParameter (int parameterID, Information &info)
{
  if (parameterID == 1) {    
    this->rho = info.theDouble;
    return 0;
  }
  else
    return -1;
}

int
ForceBeamColumn3d::activateParameter(int passedParameterID)
{
  parameterID = passedParameterID;

  return 0;  
}

const Matrix&
ForceBeamColumn3d::getKiSensitivity(int gradNumber)
{
  theMatrix.Zero();
  return theMatrix;
}

const Matrix&
ForceBeamColumn3d::getMassSensitivity(int gradNumber)
{
    theMatrix.Zero();

    double L = crdTransf->getInitialLength();
    if (rho != 0.0 && parameterID == 1)
      theMatrix(0,0) = theMatrix(1,1) = theMatrix(2,2) =
	theMatrix(6,6) = theMatrix(7,7) = theMatrix(8,8) = 0.5*L;

    return theMatrix;
}

const Vector&
ForceBeamColumn3d::getResistingForceSensitivity(int gradNumber)
{
  static Vector dqdh(6);
  dqdh = this->computedqdh(gradNumber);

  // Transform forces
  double dp0dh[6]; 
  dp0dh[0] = 0.0; dp0dh[1] = 0.0; dp0dh[2] = 0.0;
  dp0dh[3] = 0.0; dp0dh[4] = 0.0; dp0dh[5] = 0.0;
  this->computeReactionSensitivity(dp0dh, gradNumber);
  Vector dp0dhVec(dp0dh, 6);

  static Vector P(12);
  P.Zero();

  if (crdTransf->isShapeSensitivity()) {
  // dAdh^T q
    P = crdTransf->getGlobalResistingForceShapeSensitivity(Se, dp0dhVec, gradNumber);
    // k dAdh u
    const Vector &dAdh_u = crdTransf->getBasicTrialDispShapeSensitivity();
    dqdh.addMatrixVector(1.0, kv, dAdh_u, 1.0);
  }

  // A^T (dqdh + k dAdh u)
  P += crdTransf->getGlobalResistingForce(dqdh, dp0dhVec);

  return P;
}

int
ForceBeamColumn3d::commitSensitivity(int gradNumber, int numGrads)
{
  int err = 0;

  double L = crdTransf->getInitialLength();
  double oneOverL = 1.0/L;
  
  double pts[maxNumSections];
  beamIntegr->getSectionLocations(numSections, L, pts);
  
  double wts[maxNumSections];
  beamIntegr->getSectionWeights(numSections, L, wts);

  double dLdh = crdTransf->getdLdh();

  double dptsdh[maxNumSections];
  beamIntegr->getLocationsDeriv(numSections, L, dLdh, dptsdh);

  double d1oLdh = crdTransf->getd1overLdh();

  static Vector dqdh(6);
  dqdh = this->computedqdh(gradNumber);

  // dvdh = A dudh + dAdh u
  const Vector &dvdh = crdTransf->getBasicDisplSensitivity(gradNumber);
  dqdh.addMatrixVector(1.0, kv, dvdh, 1.0);  // A dudh

  if (crdTransf->isShapeSensitivity()) {
    //const Vector &dAdh_u = crdTransf->getBasicTrialDispShapeSensitivity(gradNumber);
    //dqdh.addMatrixVector(1.0, kv, dAdh_u, 1.0);  // dAdh u
  }

  // Loop over integration points
  for (int i = 0; i < numSections; i++) {

    int order = sections[i]->getOrder();
    const ID &code = sections[i]->getType();
    
    double xL  = pts[i];
    double xL1 = xL-1.0;

    double dxLdh  = dptsdh[i];    

    Vector ds(workArea, order);
    ds.Zero();

    // Add sensitivity wrt element loads
    if (numEleLoads > 0) {
      this->computeSectionForceSensitivity(ds, i, gradNumber);
    }

    int j;
    for (j = 0; j < order; j++) {
      switch(code(j)) {
      case SECTION_RESPONSE_P:
	ds(j) += dqdh(0);
	break;
      case SECTION_RESPONSE_MZ:
	ds(j) += xL1*dqdh(1) + xL*dqdh(2);
	break;
      case SECTION_RESPONSE_VY:
	ds(j) += oneOverL*(dqdh(1)+dqdh(2));
	break;
      case SECTION_RESPONSE_MY:
	ds(j) += xL1*dqdh(3) + xL*dqdh(4);
	break;
      case SECTION_RESPONSE_VZ:
	ds(j) += oneOverL*(dqdh(3)+dqdh(4));
	break;
      case SECTION_RESPONSE_T:
	ds(j) += dqdh(5);
	break;
      default:
	ds(j) += 0.0;
	break;
      }
    }

    const Vector &dsdh = sections[i]->getStressResultantSensitivity(gradNumber,true);
    ds -= dsdh;

    for (j = 0; j < order; j++) {
      switch (code(j)) {
      case SECTION_RESPONSE_MZ:
	ds(j) += dxLdh*(Se(1)+Se(2));
	break;
      case SECTION_RESPONSE_VY:
	ds(j) += d1oLdh*(Se(1)+Se(2));
	break;
      case SECTION_RESPONSE_MY:
	ds(j) += dxLdh*(Se(3)+Se(4));
	break;
      case SECTION_RESPONSE_VZ:
	ds(j) += d1oLdh*(Se(3)+Se(4));
	break;
      default:
	break;
      }
    }

    Vector de(&workArea[order], order);
    const Matrix &fs = sections[i]->getSectionFlexibility();
    de.addMatrixVector(0.0, fs, ds, 1.0);

    err += sections[i]->commitSensitivity(de, gradNumber, numGrads);
  }

  return err;
}

const Vector &
ForceBeamColumn3d::computedqdh(int gradNumber)
{
  //opserr << "FBC3d::computedqdh " << gradNumber << endln;

  double L = crdTransf->getInitialLength();
  double oneOverL = 1.0/L;
  
  double pts[maxNumSections];
  beamIntegr->getSectionLocations(numSections, L, pts);
  
  double wts[maxNumSections];
  beamIntegr->getSectionWeights(numSections, L, wts);

  double dLdh = crdTransf->getdLdh();

  double dptsdh[maxNumSections];
  beamIntegr->getLocationsDeriv(numSections, L, dLdh, dptsdh);

  double dwtsdh[maxNumSections];
  beamIntegr->getWeightsDeriv(numSections, L, dLdh, dwtsdh);

  double d1oLdh = crdTransf->getd1overLdh();

  static Vector dvdh(6);
  dvdh.Zero();

  // Loop over the integration points
  for (int i = 0; i < numSections; i++) {

    int order = sections[i]->getOrder();
    const ID &code = sections[i]->getType();
    
    double xL  = pts[i];
    double xL1 = xL-1.0;
    double wtL = wts[i]*L;
    
    double dxLdh  = dptsdh[i];// - xL/L*dLdh;
    double dwtLdh = wts[i]*dLdh + dwtsdh[i]*L;

    //opserr << dptsdh[i] << ' ' << dwtsdh[i] << endln;

    // Get section stress resultant gradient
    Vector dsdh(&workArea[order], order);
    dsdh = sections[i]->getStressResultantSensitivity(gradNumber,true);
    //opserr << "FBC2d::dqdh -- " << gradNumber << ' ' << dsdh;
    
    Vector dspdh(&workArea[2*order], order);
    dspdh.Zero();
    // Add sensitivity wrt element loads
    if (numEleLoads > 0) {
      this->computeSectionForceSensitivity(dspdh, i, gradNumber);
      //opserr << "FBC2d::dspdh -- " << i << ' ' << dsdh;
    }
    dsdh.addVector(1.0, dspdh, -1.0);

    int j;
    for (j = 0; j < order; j++) {
      switch (code(j)) {
      case SECTION_RESPONSE_MZ:
	dsdh(j) -= dxLdh*(Se(1)+Se(2));
	break;
      case SECTION_RESPONSE_VY:
	dsdh(j) -= d1oLdh*(Se(1)+Se(2));
	break;
      case SECTION_RESPONSE_MY:
	dsdh(j) -= dxLdh*(Se(3)+Se(4));
	break;
      case SECTION_RESPONSE_VZ:
	dsdh(j) -= d1oLdh*(Se(3)+Se(4));
	break;
      default:
	break;
      }
    }

    Vector dedh(workArea, order);
    const Matrix &fs = sections[i]->getSectionFlexibility();
    dedh.addMatrixVector(0.0, fs, dsdh, 1.0);

    for (j = 0; j < order; j++) {
      double dei = dedh(j)*wtL;
      switch(code(j)) {
      case SECTION_RESPONSE_P:
	dvdh(0) += dei; 
	break;
      case SECTION_RESPONSE_MZ:
	dvdh(1) += xL1*dei; 
	dvdh(2) += xL*dei;
	break;
      case SECTION_RESPONSE_VY:
	dei = oneOverL*dei;
	dvdh(1) += dei;
	dvdh(2) += dei;
	break;
      case SECTION_RESPONSE_MY:
	dvdh(3) += xL1*dei; 
	dvdh(4) += xL*dei;
	break;
      case SECTION_RESPONSE_VZ:
	dei = oneOverL*dei;
	dvdh(3) += dei;
	dvdh(4) += dei;
	break;
      case SECTION_RESPONSE_T:
	dvdh(5) += dei; 
	break;
      default:
	break;
      }
    }

    const Vector &e = vs[i];
    for (j = 0; j < order; j++) {
      switch(code(j)) {
      case SECTION_RESPONSE_P:
	dvdh(0) -= e(j)*dwtLdh;
	break;
      case SECTION_RESPONSE_MZ:
	dvdh(1) -= xL1*e(j)*dwtLdh;
	dvdh(2) -= xL*e(j)*dwtLdh;
	
	dvdh(1) -= dxLdh*e(j)*wtL;
	dvdh(2) -= dxLdh*e(j)*wtL;
	break;
      case SECTION_RESPONSE_VY:
	dvdh(1) -= oneOverL*e(j)*dwtLdh;
	dvdh(2) -= oneOverL*e(j)*dwtLdh;

	dvdh(1) -= d1oLdh*e(j)*wtL;
	dvdh(2) -= d1oLdh*e(j)*wtL;
	break;
      case SECTION_RESPONSE_MY:
	dvdh(3) -= xL1*e(j)*dwtLdh;
	dvdh(4) -= xL*e(j)*dwtLdh;
	
	dvdh(3) -= dxLdh*e(j)*wtL;
	dvdh(4) -= dxLdh*e(j)*wtL;
	break;
      case SECTION_RESPONSE_VZ:
	dvdh(3) -= oneOverL*e(j)*dwtLdh;
	dvdh(4) -= oneOverL*e(j)*dwtLdh;

	dvdh(3) -= d1oLdh*e(j)*wtL;
	dvdh(4) -= d1oLdh*e(j)*wtL;
	break;
      case SECTION_RESPONSE_T:
	dvdh(5) -= e(j)*dwtLdh;
	break;
      default:
	break;
      }
    }
  }

  static Matrix dfedh(6,6);
  dfedh.Zero();

  //opserr << "dfedh: " << dfedh << endln;

  static Vector dqdh(6);
  dqdh.addMatrixVector(0.0, kv, dvdh, 1.0);
  
  //opserr << "dqdh: " << dqdh << endln;

  return dqdh;
}

const Matrix&
ForceBeamColumn3d::computedfedh(int gradNumber)
{
  static Matrix dfedh(6,6);

  dfedh.Zero();

  double L = crdTransf->getInitialLength();
  double oneOverL  = 1.0/L;  

  double dLdh = crdTransf->getdLdh();
  double d1oLdh = crdTransf->getd1overLdh();

  double xi[maxNumSections];
  beamIntegr->getSectionLocations(numSections, L, xi);
  
  double wt[maxNumSections];
  beamIntegr->getSectionWeights(numSections, L, wt);

  double dptsdh[maxNumSections];
  beamIntegr->getLocationsDeriv(numSections, L, dLdh, dptsdh);

  double dwtsdh[maxNumSections];
  beamIntegr->getWeightsDeriv(numSections, L, dLdh, dwtsdh);

  for (int i = 0; i < numSections; i++) {

    int order      = sections[i]->getOrder();
    const ID &code = sections[i]->getType();
    
    Matrix fb(workArea, order, NEBD);
    Matrix fb2(&workArea[order*NEBD], order, NEBD);

    double xL  = xi[i];
    double xL1 = xL-1.0;
    double wtL = wt[i]*L;

    double dxLdh  = dptsdh[i];
    double dwtLdh = wt[i]*dLdh + dwtsdh[i]*L;

    const Matrix &fs = sections[i]->getInitialFlexibility();
    const Matrix &dfsdh = sections[i]->getInitialFlexibilitySensitivity(gradNumber);
    fb.Zero();
    fb2.Zero();

    double tmp;
    int ii, jj;
    for (ii = 0; ii < order; ii++) {
      switch(code(ii)) {
      case SECTION_RESPONSE_P:
	for (jj = 0; jj < order; jj++) {
	  fb(jj,0) += dfsdh(jj,ii)*wtL; // 1

	  //fb(jj,0) += fs(jj,ii)*dwtLdh; // 3

	  //fb2(jj,0) += fs(jj,ii)*wtL; // 4
	}
	break;
      case SECTION_RESPONSE_MZ:
	for (jj = 0; jj < order; jj++) {
	  tmp = dfsdh(jj,ii)*wtL; // 1
	  fb(jj,1) += xL1*tmp;
	  fb(jj,2) += xL*tmp;

	  tmp = fs(jj,ii)*wtL; // 2
	  //fb(jj,1) += dxLdh*tmp;
	  //fb(jj,2) += dxLdh*tmp;

	  tmp = fs(jj,ii)*dwtLdh; // 3
	  //fb(jj,1) += xL1*tmp;
	  //fb(jj,2) += xL*tmp;

	  tmp = fs(jj,ii)*wtL; // 4
	  //fb2(jj,1) += xL1*tmp;
	  //fb2(jj,2) += xL*tmp;
	}
	break;
      case SECTION_RESPONSE_VY:
	for (jj = 0; jj < order; jj++) {
	  tmp = oneOverL*dfsdh(jj,ii)*wtL;
	  fb(jj,1) += tmp;
	  fb(jj,2) += tmp;

	  // Need to complete for dLdh != 0
	}
	break;
      default:
	break;
      }
    }
    for (ii = 0; ii < order; ii++) {
      switch (code(ii)) {
      case SECTION_RESPONSE_P:
	for (jj = 0; jj < NEBD; jj++)
	  dfedh(0,jj) += fb(ii,jj);
	break;
      case SECTION_RESPONSE_MZ:
	for (jj = 0; jj < NEBD; jj++) {
	  tmp = fb(ii,jj); // 1,2,3
	  dfedh(1,jj) += xL1*tmp;
	  dfedh(2,jj) += xL*tmp;

	  tmp = fb2(ii,jj); // 4
	  //dfedh(1,jj) += dxLdh*tmp;
	  //dfedh(2,jj) += dxLdh*tmp;
	}
	break;
      case SECTION_RESPONSE_VY:
	for (jj = 0; jj < NEBD; jj++) {
	  tmp = oneOverL*fb(ii,jj);
	  dfedh(1,jj) += tmp;
	  dfedh(2,jj) += tmp;

	  // Need to complete for dLdh != 0
	}
	break;
      default:
	break;
      }
    }
  }
  
  return dfedh;
}

void
ForceBeamColumn3d::setSectionPointers(int numSec, SectionForceDeformation **secPtrs)
{
  if (numSec > maxNumSections) {
    opserr << "Error: ForceBeamColumn3d::setSectionPointers -- max number of sections exceeded";
  }
  
  numSections = numSec;
  
  if (secPtrs == 0) {
    opserr << "Error: ForceBeamColumn3d::setSectionPointers -- invalid section pointer";
  }	  
  
  sections = new SectionForceDeformation *[numSections];
  if (sections == 0) {
    opserr << "Error: ForceBeamColumn3d::setSectionPointers -- could not allocate section pointers";
  }  
  
  for (int i = 0; i < numSections; i++) {
    
    if (secPtrs[i] == 0) {
      opserr << "Error: ForceBeamColumn3d::setSectionPointers -- null section pointer " << i << endln;
    }
    
    sections[i] = secPtrs[i]->getCopy();
    
    if (sections[i] == 0) {
      opserr << "Error: ForceBeamColumn3d::setSectionPointers -- could not create copy of section " << i << endln;
    }

    int order = sections[i]->getOrder();
    const ID &code = sections[i]->getType();
    for (int j = 0; j < order; j++) {
      if (code(j) == SECTION_RESPONSE_T)
	isTorsion = true;
    }
  }
  
  if (!isTorsion)
    opserr << "ForceBeamColumn3d::ForceBeamColumn3d -- no torsion detected in sections, " <<
      "continuing with element torsional stiffness GJ/L = " << 1.0/DefaultLoverGJ;

  // allocate section flexibility matrices and section deformation vectors
  fs  = new Matrix [numSections];
  if (fs == 0) {
    opserr << "ForceBeamColumn3d::setSectionPointers -- failed to allocate fs array";
  }
  
  vs = new Vector [numSections];
  if (vs == 0) {
    opserr << "ForceBeamColumn3d::setSectionPointers -- failed to allocate vs array";
  }
  
  Ssr  = new Vector [numSections];
  if (Ssr == 0) {
    opserr << "ForceBeamColumn3d::setSectionPointers -- failed to allocate Ssr array";
  }
  
  vscommit = new Vector [numSections];
  if (vscommit == 0) {
    opserr << "ForceBeamColumn3d::setSectionPointers -- failed to allocate vscommit array";   
  }
  
}
