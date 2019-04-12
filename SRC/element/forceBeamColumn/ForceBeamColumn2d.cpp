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
                                                                        
// $Revision: 1.42 $
// $Date: 2010-05-13 00:16:33 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/forceBeamColumn/ForceBeamColumn2d.cpp,v $

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
#include <ForceBeamColumn2d.h>
#include <MatrixUtil.h>
#include <Domain.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Renderer.h>
#include <math.h>
#include <elementAPI.h>
#include <string>
#include <ElementResponse.h>
#include <CompositeResponse.h>
#include <ElementalLoad.h>
#include <ElementIter.h>
#include <map>

Matrix ForceBeamColumn2d::theMatrix(6,6);
Vector ForceBeamColumn2d::theVector(6);
double ForceBeamColumn2d::workArea[200];

Vector *ForceBeamColumn2d::vsSubdivide = 0;
Matrix *ForceBeamColumn2d::fsSubdivide = 0;
Vector *ForceBeamColumn2d::SsrSubdivide = 0;

void* OPS_ForceBeamColumn2d()
{
    if(OPS_GetNumRemainingInputArgs() < 5) {
	opserr<<"insufficient arguments:eleTag,iNode,jNode,transfTag,integrationTag\n";
	return 0;
    }

    int ndm = OPS_GetNDM();
    int ndf = OPS_GetNDF();
    if(ndm != 2 || ndf != 3) {
	opserr<<"ndm must be 2 and ndf must be 3\n";
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
    double mass = 0.0, tol=1e-12;
    int maxIter = 10;
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
	} else if(strcmp(type,"-mass") == 0) {
	    if(OPS_GetNumRemainingInputArgs() > 0) {
		if(OPS_GetDoubleInput(&numData,&mass) < 0) {
		    opserr << "WARNING invalid mass\n";
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
	    delete [] sections;
	    return 0;
	}
    }

    Element *theEle =  new ForceBeamColumn2d(iData[0],iData[1],iData[2],secTags.Size(),sections,
					     *bi,*theTransf,mass,maxIter,tol);
    delete [] sections;
    return theEle;
}

void* OPS_ForceBeamColumn2d(const ID &info)
{
    // data
    int iData[5];
    int numData;
    double mass = 0.0, tol=1e-12;
    int maxIter = 10;

    // regular element, not in a mesh, get tags
    if (info.Size() == 0) {
	if(OPS_GetNumRemainingInputArgs() < 5) {
	    opserr<<"insufficient arguments:eleTag,iNode,jNode,transfTag,integrationTag\n";
	    return 0;
	}

	int ndm = OPS_GetNDM();
	int ndf = OPS_GetNDF();
	if(ndm != 2 || ndf != 3) {
	    opserr<<"ndm must be 2 and ndf must be 3\n";
	    return 0;
	}

	// inputs:
	numData = 3;
	if(OPS_GetIntInput(&numData,&iData[0]) < 0) {
	    opserr << "WARNING invalid int inputs\n";
	    return 0;
	}
    }

    // regular element, or in a mesh
    if (info.Size()==0 || info(0)==1) {
	if(OPS_GetNumRemainingInputArgs() < 2) {
	    opserr<<"insufficient arguments: transfTag,integrationTag\n";
	    return 0;
	}

	numData = 2;
	if(OPS_GetIntInput(&numData,&iData[3]) < 0) {
	    opserr << "WARNING invalid int inputs\n";
	    return 0;
	}

	// options
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
	    } else if(strcmp(type,"-mass") == 0) {
		if(OPS_GetNumRemainingInputArgs() > 0) {
		    if(OPS_GetDoubleInput(&numData,&mass) < 0) {
			opserr << "WARNING invalid mass\n";
			return 0;
		    }
		}
	    }
	}
    }

    // store data for different mesh
    static std::map<int, Vector> meshdata;
    if (info.Size()>0 && info(0)==1) {
	if (info.Size() < 2) {
	    opserr << "WARNING: need info -- inmesh, meshtag\n";
	    return 0;
	}

	// save the data for a mesh
	Vector& mdata = meshdata[info(1)];
	mdata.resize(5);
	mdata(0) = iData[3];
	mdata(1) = iData[4];
	mdata(2) = mass;
	mdata(3) = tol;
	mdata(4) = maxIter;
	return &meshdata;

    } else if (info.Size()>0 && info(0)==2) {
	if (info.Size() < 5) {
	    opserr << "WARNING: need info -- inmesh, meshtag, eleTag, nd1, nd2\n";
	    return 0;
	}

	// get the data for a mesh
	Vector& mdata = meshdata[info(1)];
	if (mdata.Size() < 5) return 0;

	iData[0] = info(2);
	iData[1] = info(3);
	iData[2] = info(4);
	iData[3] = mdata(0);
	iData[4] = mdata(1);
	mass = mdata(2);
	tol = mdata(3);
	maxIter = mdata(4);
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
	    delete [] sections;
	    return 0;
	}
    }

    Element *theEle =  new ForceBeamColumn2d(iData[0],iData[1],iData[2],secTags.Size(),sections,
					     *bi,*theTransf,mass,maxIter,tol);
    delete [] sections;
    return theEle;
}

int OPS_ForceBeamColumn2d(Domain& theDomain, const ID& elenodes, ID& eletags)
{
    if(OPS_GetNumRemainingInputArgs() < 2) {
	opserr<<"insufficient arguments:transfTag,integrationTag\n";
	return -1;
    }

    // inputs: 
    int iData[2];
    int numData = 2;
    if(OPS_GetIntInput(&numData,&iData[0]) < 0) return -1;

    // options
    double mass = 0.0, tol=1e-12;
    int maxIter = 10;
    numData = 1;
    while(OPS_GetNumRemainingInputArgs() > 0) {
	const char* type = OPS_GetString();
	if(strcmp(type,"-iter") == 0) {
	    if(OPS_GetNumRemainingInputArgs() > 1) {
		if(OPS_GetIntInput(&numData,&maxIter) < 0) return -1;
		if(OPS_GetDoubleInput(&numData,&tol) < 0) return -1;
	    }
	} else if(strcmp(type,"-mass") == 0) {
	    if(OPS_GetNumRemainingInputArgs() > 0) {
		if(OPS_GetDoubleInput(&numData,&mass) < 0) return -1;
	    }
	}
    }

    // check transf
    CrdTransf* theTransf = OPS_getCrdTransf(iData[0]);
    if(theTransf == 0) {
	opserr<<"coord transfomration not found\n";
	return -1;
    }

    // check beam integrataion
    BeamIntegrationRule* theRule = OPS_getBeamIntegrationRule(iData[1]);
    if(theRule == 0) {
	opserr<<"beam integration not found\n";
	return -1;
    }
    BeamIntegration* bi = theRule->getBeamIntegration();
    if(bi == 0) {
	opserr<<"beam integration is null\n";
	return -1;
    }

    // check sections
    const ID& secTags = theRule->getSectionTags();
    SectionForceDeformation** sections = new SectionForceDeformation *[secTags.Size()];
    for(int i=0; i<secTags.Size(); i++) {
	sections[i] = OPS_getSectionForceDeformation(secTags(i));
	if(sections[i] == 0) {
	    opserr<<"section "<<secTags(i)<<"not found\n";
		delete [] sections;
	    return -1;
	}
    }

    // create elements
    ElementIter& theEles = theDomain.getElements();
    Element* theEle = theEles();
    int currTag = 0;
    if (theEle != 0) {
	currTag = theEle->getTag();
    }
    eletags.resize(elenodes.Size()/2);
    for (int i=0; i<elenodes.Size()/2; i++) {
	theEle = new ForceBeamColumn2d(--currTag,elenodes(2*i),elenodes(2*i+1),secTags.Size(),
				       sections,*bi,*theTransf,mass,maxIter,tol);
	if (theEle == 0) {
	    opserr<<"WARNING: run out of memory for creating element\n";
	    return -1;
	}
	if (theDomain.addElement(theEle) == false) {
	    opserr<<"WARNING: failed to add element to domain\n";
	    delete theEle;
	    return -1;
	}
	eletags(i) = currTag;
    }
    
    delete [] sections;
    return 0;
}

// constructor:
// invoked by a FEM_ObjectBroker, recvSelf() needs to be invoked on this object.
ForceBeamColumn2d::ForceBeamColumn2d(): 
  Element(0,ELE_TAG_ForceBeamColumn2d), connectedExternalNodes(2), 
  beamIntegr(0), numSections(0), sections(0), crdTransf(0),
  rho(0.0), maxIters(0), tol(0.0),
  initialFlag(0),
  kv(NEBD,NEBD), Se(NEBD),
  kvcommit(NEBD,NEBD), Secommit(NEBD),
  fs(0), vs(0), Ssr(0), vscommit(0), 
  numEleLoads(0), sizeEleLoads(0), eleLoads(0), eleLoadFactors(0), load(6),
  Ki(0), parameterID(0)
{
  load.Zero();

  theNodes[0] = 0;  
  theNodes[1] = 0;

  if (vsSubdivide == 0)
    vsSubdivide  = new Vector [maxNumSections];
  if (fsSubdivide == 0)
    fsSubdivide  = new Matrix [maxNumSections];
  if (SsrSubdivide == 0)
    SsrSubdivide  = new Vector [maxNumSections];
  if (!vsSubdivide || !fsSubdivide || !SsrSubdivide) {
    opserr << "ForceBeamColumn2d::ForceBeamColumn2d() -- failed to allocate Subdivide arrays";   
    exit(-1);
  }
}

// constructor which takes the unique element tag, sections,
// and the node ID's of its nodal end points. 
// allocates the necessary space needed by each object
ForceBeamColumn2d::ForceBeamColumn2d (int tag, int nodeI, int nodeJ,
				      int numSec, SectionForceDeformation **sec,
				      BeamIntegration &bi,
				      CrdTransf &coordTransf, double massDensPerUnitLength,
				      int maxNumIters, double tolerance):
  Element(tag,ELE_TAG_ForceBeamColumn2d), connectedExternalNodes(2),
  beamIntegr(0), numSections(0), sections(0), crdTransf(0),
  rho(massDensPerUnitLength),maxIters(maxNumIters), tol(tolerance), 
  initialFlag(0),
  kv(NEBD,NEBD), Se(NEBD), 
  kvcommit(NEBD,NEBD), Secommit(NEBD),
  fs(0), vs(0),Ssr(0), vscommit(0), 
  numEleLoads(0), sizeEleLoads(0), eleLoads(0), eleLoadFactors(0), load(6),
  Ki(0), parameterID(0)
{
  load.Zero();
  
  theNodes[0] = 0;
  theNodes[1] = 0;

  connectedExternalNodes(0) = nodeI;
  connectedExternalNodes(1) = nodeJ;    
  
  beamIntegr = bi.getCopy();
  if (beamIntegr == 0) {
    opserr << "Error: ForceBeamColumn2d::ForceBeamColumn2d: could not create copy of beam integration object" << endln;
    exit(-1);
  }
  
  // get copy of the transformation object   
  crdTransf = coordTransf.getCopy2d(); 
  if (crdTransf == 0) {
    opserr << "Error: ForceBeamColumn2d::ForceBeamColumn2d: could not create copy of coordinate transformation object" << endln;
    exit(-1);
  }

  this->setSectionPointers(numSec, sec);
  
  if (vsSubdivide == 0)
    vsSubdivide  = new Vector [maxNumSections];
  if (fsSubdivide == 0)
    fsSubdivide  = new Matrix [maxNumSections];
  if (SsrSubdivide == 0)
    SsrSubdivide  = new Vector [maxNumSections];
  if (!vsSubdivide || !fsSubdivide || !SsrSubdivide) {
    opserr << "ForceBeamColumn2d::ForceBeamColumn2d() -- failed to allocate Subdivide arrays";   
    exit(-1);
  }
}

// ~ForceBeamColumn2d():
// 	destructor
//      delete must be invoked on any objects created by the object
ForceBeamColumn2d::~ForceBeamColumn2d()
{
  if (sections != 0) {
    for (int i=0; i < numSections; i++)
      if (sections[i] != 0)
	delete sections[i];
    delete [] sections;
  }

  if (sizeEleLoads != 0) {
    if (eleLoads != 0)
      delete [] eleLoads;

    if (eleLoadFactors != 0)
      delete [] eleLoadFactors;
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
}

int
ForceBeamColumn2d::getNumExternalNodes(void) const
{
  return 2;
}

const ID &
ForceBeamColumn2d::getExternalNodes(void) 
{
  return connectedExternalNodes;
}

Node **
ForceBeamColumn2d::getNodePtrs()
{
  return theNodes;
}

int
ForceBeamColumn2d::getNumDOF(void) 
{
  return NEGD;
}

void
ForceBeamColumn2d::setDomain(Domain *theDomain)
{
  // check Domain is not null - invoked when object removed from a domain
  if (theDomain == 0) {
    theNodes[0] = 0;
    theNodes[1] = 0;
    
    opserr << "ForceBeamColumn2d::setDomain:  theDomain = 0 ";
    exit(0); 
  }

  // get pointers to the nodes
  
  int Nd1 = connectedExternalNodes(0);  
  int Nd2 = connectedExternalNodes(1);
  
  theNodes[0] = theDomain->getNode(Nd1);
  theNodes[1] = theDomain->getNode(Nd2);  
  
  if (theNodes[0] == 0) {
    opserr << "ForceBeamColumn2d::setDomain: Nd1: ";
    opserr << Nd1 << "does not exist in model\n";
    exit(0);
  }
  
  if (theNodes[1] == 0) {
    opserr << "ForceBeamColumn2d::setDomain: Nd2: ";
    opserr << Nd2 << "does not exist in model\n";
    exit(0);
  }
  
  // call the DomainComponent class method 
  this->DomainComponent::setDomain(theDomain);
  
  // ensure connected nodes have correct number of dof's
  int dofNode1 = theNodes[0]->getNumberDOF();
  int dofNode2 = theNodes[1]->getNumberDOF();
  
  if ((dofNode1 != NND) || (dofNode2 != NND)) {
    opserr << "ForceBeamColumn2d::setDomain(): Nd2 or Nd1 incorrect dof ";
    exit(0);
  }
   
  // initialize the transformation
  if (crdTransf->initialize(theNodes[0], theNodes[1])) {
    opserr << "ForceBeamColumn2d::setDomain(): Error initializing coordinate transformation";  
    exit(0);
  }
    
  // get element length
  double L = crdTransf->getInitialLength();
  if (L == 0.0) {
    opserr << "ForceBeamColumn2d::setDomain(): Zero element length:" << this->getTag();  
    exit(0);
  }

  if (initialFlag == 0) 
    this->initializeSectionHistoryVariables();
}

int
ForceBeamColumn2d::commitState()
{
  int err = 0;
  int i = 0;

  // call element commitState to do any base class stuff
  if ((err = this->Element::commitState()) != 0) {
    opserr << "ForceBeamColumn2d::commitState () - failed in base class";
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
  
  return err;
}

int ForceBeamColumn2d::revertToLastCommit()
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
  
  return err;
}

int ForceBeamColumn2d::revertToStart()
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
  return err;
}


const Matrix &
ForceBeamColumn2d::getInitialStiff(void)
{
  // check for quick return
  if (Ki != 0)
    return *Ki;

  /*
  else
    Ki = new Matrix(this->getTangentStiff());
  */

  static Matrix f(NEBD, NEBD);   // element flexibility matrix  
  this->getInitialFlexibility(f);

  /*
  static Matrix I(NEBD,NEBD);   // an identity matrix for matrix inverse  
  I.Zero();
  for (int i=0; i<NEBD; i++)
    I(i,i) = 1.0;
  
  // calculate element stiffness matrix
  // invert3by3Matrix(f, kv);

  static Matrix kvInit(NEBD, NEBD);
  if (f.Solve(I, kvInit) < 0)
    opserr << "ForceBeamColumn2d::getInitialStiff() -- could not invert flexibility\n";
  */

  static Matrix kvInit(NEBD, NEBD);
  f.Invert(kvInit);
  Ki = new Matrix(crdTransf->getInitialGlobalStiffMatrix(kvInit));
  return *Ki;
}

const Matrix &
ForceBeamColumn2d::getTangentStiff(void)
{
  crdTransf->update();	// Will remove once we clean up the corotational 2d transformation -- MHS
  return crdTransf->getGlobalStiffMatrix(kv, Se);
}
    
void
ForceBeamColumn2d::computeReactions(double *p0)
{
  int type;
  double L = crdTransf->getInitialLength();
  
  for (int i = 0; i < numEleLoads; i++) {
    
    double loadFactor = eleLoadFactors[i];
    const Vector &data = eleLoads[i]->getData(type, loadFactor);

    if (type == LOAD_TAG_Beam2dUniformLoad) {
      double wa = data(1)*loadFactor;  // Axial
      double wy = data(0)*loadFactor;  // Transverse
      
      p0[0] -= wa*L;
      double V = 0.5*wy*L;
      p0[1] -= V;
      p0[2] -= V;
    }
    else if (type == LOAD_TAG_Beam2dPartialUniformLoad) {
      double waa = data(2)*loadFactor;  // Axial
      double wab = data(3)*loadFactor;  // Axial
      double wya = data(0)*loadFactor;  // Transverse
      double wyb = data(1)*loadFactor;  // Transverse      
      double a = data(4)*L;
      double b = data(5)*L;

      p0[0] -= waa*(b-a) + 0.5*(wab-waa)*(b-a);
      double Fy = wya*(b-a);
      double c = a + 0.5*(b-a);
      p0[1] -= Fy*(1-c/L);
      p0[2] -= Fy*c/L;
      Fy = 0.5*(wyb-wya)*(b-a);
      c = a + 2.0/3.0*(b-a);
      p0[1] -= Fy*(1-c/L);
      p0[2] -= Fy*c/L;      
    }
    else if (type == LOAD_TAG_Beam2dPointLoad) {
      double P = data(0)*loadFactor;
      double N = data(1)*loadFactor;
      double aOverL = data(2);
      
      if (aOverL < 0.0 || aOverL > 1.0)
	continue;
      
      double V1 = P*(1.0-aOverL);
      double V2 = P*aOverL;
      
      p0[0] -= N;
      p0[1] -= V1;
      p0[2] -= V2;
    }
  }
}

void
ForceBeamColumn2d::computeReactionSensitivity(double *dp0dh, int gradNumber)
{
  int type;
  double L = crdTransf->getInitialLength();
  
  double dLdh = crdTransf->getdLdh();

  for (int i = 0; i < numEleLoads; i++) {
    
    const Vector &data = eleLoads[i]->getData(type, 1.0);

    if (type == LOAD_TAG_Beam2dUniformLoad) {
      double wy = data(0)*1.0;  // Transverse
      double wa = data(1)*1.0;  // Axial

      const Vector &sens = eleLoads[i]->getSensitivityData(gradNumber);
      double dwydh = sens(0);
      double dwadh = sens(1);

      //p0[0] -= wa*L;
      dp0dh[0] -= wa*dLdh + dwadh*L;

      //double V = 0.5*wy*L;
      //p0[1] -= V;
      //p0[2] -= V;
      double dVdh = 0.5*(wy*dLdh + dwydh*L);
      dp0dh[1] -= dVdh;
      dp0dh[2] -= dVdh;
    }
    else if (type == LOAD_TAG_Beam2dPointLoad) {
      double P = data(0)*1.0;
      double N = data(1)*1.0;
      double aOverL = data(2);

      if (aOverL < 0.0 || aOverL > 1.0)
	continue;
      
      const Vector &sens = eleLoads[i]->getSensitivityData(gradNumber);
      double dPdh = sens(0);
      double dNdh = sens(1);
      double daLdh = sens(2);

      //double a = aOverL*L;
      
      //double V1 = P*(1.0-aOverL);
      //double V2 = P*aOverL;
      double dV1dh = P*(0.0-daLdh) + dPdh*(1.0-aOverL);
      double dV2dh = P*daLdh + dPdh*aOverL;
      
      //p0[0] -= N;
      //p0[1] -= V1;
      //p0[2] -= V2;
      dp0dh[0] -= dNdh;
      dp0dh[1] -= dV1dh;
      dp0dh[2] -= dV2dh;
    }
  }
}

const Vector &
ForceBeamColumn2d::getResistingForce(void)
{
  // Will remove once we clean up the corotational 2d transformation -- MHS
  crdTransf->update();
  
  double p0[3];
  Vector p0Vec(p0, 3);
  p0Vec.Zero();

  if (numEleLoads > 0)
    this->computeReactions(p0);
  // Compute the current resisting force
  theVector = crdTransf->getGlobalResistingForce(Se, p0Vec);
  
  if (rho != 0)
    theVector.addVector(1.0, load, -1.0);

  return theVector;
}

void
ForceBeamColumn2d::initializeSectionHistoryVariables (void)
{
  for (int i = 0; i < numSections; i++) {
    int order = sections[i]->getOrder();
    
	fs[i] = Matrix(order, order);
	vs[i] = Vector(order);
	Ssr[i] = Vector(order);
    
	vscommit[i] = Vector(order);
  }
}

/********* NEWTON , SUBDIVIDE AND INITIAL ITERATIONS ********************
 */
int
ForceBeamColumn2d::update()
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

  static double factor = 10;

  maxSubdivisions = 4;

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
    // initial tangent on first iteration then regular newton (if l==1), or 
    // initial tangent iterations (if l==2)

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

	  if (beamIntegr->addElasticFlexibility(L, f) < 0) {
	    vr(0) += f(0,0)*SeTrial(0);
	    vr(1) += f(1,1)*SeTrial(1) + f(1,2)*SeTrial(2);
	    vr(2) += f(2,1)*SeTrial(1) + f(2,2)*SeTrial(2);
	  }

	  double v0[3];
	  v0[0] = 0.0; v0[1] = 0.0; v0[2] = 0.0;

	  for (int ie = 0; ie < numEleLoads; ie++) 
	    beamIntegr->addElasticDeformations(eleLoads[ie], eleLoadFactors[ie], L, v0);

	  // Add effects of element loads
	  vr(0) += v0[0];
	  vr(1) += v0[1];
	  vr(2) += v0[2];
	  
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
		Ss(ii) =  xL1*SeTrial(1) + xL*SeTrial(2);
		break;
	      case SECTION_RESPONSE_VY:
		Ss(ii) = oneOverL*(SeTrial(1)+SeTrial(2));
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

	    if (sections[i]->setTrialSectionDeformation(vsSubdivide[i]) < 0) {
	      opserr << "ForceBeamColumn2d::update() - section failed in setTrial\n";
	      return -1;
	    }

	    // get section resisting forces
	    SsrSubdivide[i] = sections[i]->getStressResultant();

	    // get section flexibility matrix
	    fsSubdivide[i] = sections[i]->getSectionFlexibility();

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
	      default:
		break;
	      }
	    }
	  }
	  
	  // calculate element stiffness matrix
	  // invert3by3Matrix(f, kv);	  
	  if (f.Solve(I, kvTrial) < 0)
	    opserr << "ForceBeamColumn2d::update() -- could not invert flexibility\n";
				    
	  // dv = vin + dvTrial  - vr
	  dv = vin;
	  dv += dvTrial;
	  dv -= vr;
	  
	  // dv.addVector(1.0, vr, -1.0);
	  
	  // dSe = kv * dv;
	  dSe.addMatrixVector(0.0, kvTrial, dv, 1.0);

	  dW = dv ^ dSe; 
	  
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
	      //opserr << dvToDo << dvTrial << endln;
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


  // if fail to converge we return an error flag & print an error message

  if (converged == false) {
    opserr << "WARNING - ForceBeamColumn2d::update - failed to get compatible ";
    opserr << "element forces & deformations for element: ";
    opserr << this->getTag() << "(dW: << " << dW << ")\n";
    return -1;
  }

  initialFlag = 1;

  return 0;
}

void ForceBeamColumn2d::getForceInterpolatMatrix(double xi, Matrix &b, const ID &code)
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
    default:
      break;
    }
  }
}

void ForceBeamColumn2d::getDistrLoadInterpolatMatrix(double xi, Matrix &bp, const ID &code)
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
    default:
      break;
    }
  }
}

const Matrix &
ForceBeamColumn2d::getMass(void)
{ 
  theMatrix.Zero();
  
  double L = crdTransf->getInitialLength();
  if (rho != 0.0)
    theMatrix(0,0) = theMatrix(1,1) = theMatrix(3,3) = theMatrix(4,4) = 0.5*L*rho;
  
  return theMatrix;
}

void 
ForceBeamColumn2d::zeroLoad(void)
{
  load.Zero();
  
  // This is a semi-hack -- MHS
  numEleLoads = 0;

  return;
}

int
ForceBeamColumn2d::addLoad(ElementalLoad *theLoad, double loadFactor)
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
ForceBeamColumn2d::computeSectionForces(Vector &sp, int isec)
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

    if (type == LOAD_TAG_Beam2dUniformLoad) {
      double wa = data(1)*loadFactor;  // Axial
      double wy = data(0)*loadFactor;  // Transverse
      
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
	default:
	  break;
	}
      }
    }
    else if (type == LOAD_TAG_Beam2dPartialUniformLoad) {
      double waa = data(2)*loadFactor;  // Axial
      double wab = data(3)*loadFactor;  // Axial
      double wya = data(0)*loadFactor;  // Transverse
      double wyb = data(1)*loadFactor;  // Transverse
      double a = data(4)*L;
      double b = data(5)*L;

      double Fa = waa*(b-a) + 0.5*(wab-waa)*(b-a); // resultant axial load
      double Fy = wya*(b-a); // resultant transverse load
      double c = a + 0.5*(b-a);
      double VI = Fy*(1-c/L);
      double VJ = Fy*c/L;
      Fy = 0.5*(wyb-wya)*(b-a); // resultant transverse load
      c = a + 2.0/3.0*(b-a);
      VI += Fy*(1-c/L);
      VJ += Fy*c/L;      

      for (int ii = 0; ii < order; ii++) {
	
	if (x <= a) {
	  switch(code(ii)) {
	  case SECTION_RESPONSE_P:
	    sp(ii) += Fa;
	    break;
	  case SECTION_RESPONSE_MZ:
	    sp(ii) -= VI*x;
	    break;
	  case SECTION_RESPONSE_VY:
	    sp(ii) -= VI;
	    break;
	  default:
	    break;
	  }
	}
	else if (x >= b) {
	  switch(code(ii)) {
	  case SECTION_RESPONSE_MZ:
	    sp(ii) += VJ*(x-L);
	    break;
	  case SECTION_RESPONSE_VY:
	    sp(ii) += VJ;
	    break;
	  default:
	    break;
	  }
	}
	else {
	  double wx = wya + (wyb-wya)/(b-a)*(x-a);
	  switch(code(ii)) {
	  case SECTION_RESPONSE_P:
	    sp(ii) += Fa - waa*(x-a) - 0.5*(wab-waa)/(b-a)*(x-a)*(x-a);
	    break;
	  case SECTION_RESPONSE_MZ:
	    sp(ii) += -VI*x + wya*(x-a)*0.5*(x-a) + 0.5*(wx-wya)*(x-a)*(x-a)/3.0;
	    break;
	  case SECTION_RESPONSE_VY:
	    sp(ii) += -VI + wya*(x-a) + 0.5*(wx-wya)*(x-a);
	    break;
	  default:
	    break;
	  }
	}
      }
    }
    else if (type == LOAD_TAG_Beam2dPointLoad) {
      double P = data(0)*loadFactor;
      double N = data(1)*loadFactor;
      double aOverL = data(2);
      
      if (aOverL < 0.0 || aOverL > 1.0)
	continue;
      
      double a = aOverL*L;
      
      double V1 = P*(1.0-aOverL);
      double V2 = P*aOverL;
      
      for (int ii = 0; ii < order; ii++) {
	
	if (x <= a) {
	  switch(code(ii)) {
	  case SECTION_RESPONSE_P:
	    sp(ii) += N;
	    break;
	  case SECTION_RESPONSE_MZ:
	    sp(ii) -= x*V1;
	    break;
	  case SECTION_RESPONSE_VY:
	    sp(ii) -= V1;
	    break;
	  default:
	    break;
	  }
	}
	else {
	  switch(code(ii)) {
	  case SECTION_RESPONSE_MZ:
	    sp(ii) -= (L-x)*V2;
	    break;
	  case SECTION_RESPONSE_VY:
	    sp(ii) += V2;
	    break;
	  default:
	    break;
	  }
	}
      }
    }
    else {
      opserr << "ForceBeamColumn2d::addLoad -- load type unknown for element with tag: " <<
	this->getTag() << endln;
    }
  }
  
  // Don't think we need to do this anymore -- MHS
  //this->update(); // quick fix -- MHS
}

void
ForceBeamColumn2d::computeSectionForceSensitivity(Vector &dspdh, int isec,
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
    
    if (type == LOAD_TAG_Beam2dUniformLoad) {
      double wy = data(0)*1.0;  // Transverse
      double wa = data(1)*1.0;  // Axial

      const Vector &sens = eleLoads[i]->getSensitivityData(gradNumber);
      double dwydh = sens(0);
      double dwadh = sens(1);
      //opserr << wy << ' ' << dwydh << endln;
      for (int ii = 0; ii < order; ii++) {
	
	switch(code(ii)) {
	case SECTION_RESPONSE_P:
	  //sp(ii) += wa*(L-x);
	  dspdh(ii) += dwadh*(L-x) + wa*(dLdh-dxdh);
	  break;
	case SECTION_RESPONSE_MZ:
	  //sp(ii) += wy*0.5*x*(x-L);
	  //dspdh(ii) += 0.5 * (dwydh*x*(x-L) + wy*dxdh*(x-L) + wy*x*(dxdh-dLdh));
	  dspdh(ii) += 0.5 * (dwydh*x*(x-L) + wy*(dxdh*(2*x-L)-x*dLdh));
	  break;
	case SECTION_RESPONSE_VY:
	  //sp(ii) += wy*(x-0.5*L);
	  dspdh(ii) += dwydh*(x-0.5*L) + wy*(dxdh-0.5*dLdh);
	  break;
	default:
	  break;
	}
      }
    }
    else if (type == LOAD_TAG_Beam2dPartialUniformLoad) {
      // Needs to be done
      // Do nothing for now
    }
    else if (type == LOAD_TAG_Beam2dPointLoad) {
      double P = data(0)*1.0;
      double N = data(1)*1.0;
      double aOverL = data(2);

      if (aOverL < 0.0 || aOverL > 1.0)
	continue;
      
      const Vector &sens = eleLoads[i]->getSensitivityData(gradNumber);
      double dPdh = sens(0);
      double dNdh = sens(1);
      double daLdh = sens(2);

      double a = aOverL*L;

      double V1 = P*(1.0-aOverL);
      double V2 = P*aOverL;
      double dV1dh = P*(0.0-daLdh) + dPdh*(1.0-aOverL);
      double dV2dh = P*daLdh + dPdh*aOverL;

      for (int ii = 0; ii < order; ii++) {
	
	if (x <= a) {
	  switch(code(ii)) {
	  case SECTION_RESPONSE_P:
	    //sp(ii) += N;
	    dspdh(ii) += dNdh;
	    break;
	  case SECTION_RESPONSE_MZ:
	    //sp(ii) -= x*V1;
	    dspdh(ii) -= (dxdh*V1 + x*dV1dh);
	    break;
	  case SECTION_RESPONSE_VY:
	    //sp(ii) -= V1;
	    dspdh(ii) -= dV1dh;
	    break;
	  default:
	    break;
	  }
	}
	else {
	  switch(code(ii)) {
	  case SECTION_RESPONSE_MZ:
	    //sp(ii) -= (L-x)*V2;
	    dspdh(ii) -= (dLdh-dxdh)*V2 + (L-x)*dV2dh;
	    break;
	  case SECTION_RESPONSE_VY:
	    //sp(ii) += V2;
	    dspdh(ii) += dV2dh;
	    break;
	  default:
	    break;
	  }
	}
      }
    }
    else {
      opserr << "ForceBeamColumn2d::computeSectionForceSensitivity -- load type unknown for element with tag: " <<
	this->getTag() << endln;
    }
  }
}

int 
ForceBeamColumn2d::addInertiaLoadToUnbalance(const Vector &accel)
{
  // Check for a quick return
  if (rho == 0.0)
    return 0;

  // get R * accel from the nodes
  const Vector &Raccel1 = theNodes[0]->getRV(accel);
  const Vector &Raccel2 = theNodes[1]->getRV(accel);    

  double L = crdTransf->getInitialLength();
  double m = 0.5*rho*L;

  // Should be done through p0[0]
  load(0) -= m*Raccel1(0);
  load(1) -= m*Raccel1(1);
  load(3) -= m*Raccel2(0);
  load(4) -= m*Raccel2(1);

  return 0;
}

const Vector &
ForceBeamColumn2d::getResistingForceIncInertia()
{	
  // Compute the current resisting force
  theVector = this->getResistingForce();

  // Check for a quick return
  if (rho != 0.0) {
    const Vector &accel1 = theNodes[0]->getTrialAccel();
    const Vector &accel2 = theNodes[1]->getTrialAccel();
    
    double L = crdTransf->getInitialLength();
    double m = 0.5*rho*L;

    theVector(0) += m*accel1(0);
    theVector(1) += m*accel1(1);
    theVector(3) += m*accel2(0);
    theVector(4) += m*accel2(1);
    
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
ForceBeamColumn2d::sendSelf(int commitTag, Channel &theChannel)
{  
  // place the integer data into an ID
  int dbTag = this->getDbTag();
  int i, j , k;
  int loc = 0;

  static ID idData(11);  // one bigger than needed so no clash later
  idData(0) = this->getTag();
  idData(1) = connectedExternalNodes(0);
  idData(2) = connectedExternalNodes(1);
  idData(3) = numSections;
  idData(4) = maxIters;
  idData(5) = initialFlag;

  idData(6) = crdTransf->getClassTag();
  int crdTransfDbTag  = crdTransf->getDbTag();
  if (crdTransfDbTag  == 0) {
    crdTransfDbTag = theChannel.getDbTag();
    if (crdTransfDbTag  != 0) 
      crdTransf->setDbTag(crdTransfDbTag);
  }
  idData(7) = crdTransfDbTag;

  idData(8) = beamIntegr->getClassTag();
  int beamIntegrDbTag  = beamIntegr->getDbTag();
  if (beamIntegrDbTag  == 0) {
    beamIntegrDbTag = theChannel.getDbTag();
    if (beamIntegrDbTag  != 0) 
      beamIntegr->setDbTag(beamIntegrDbTag);
  }
  idData(9) = beamIntegrDbTag;

  if (theChannel.sendID(dbTag, commitTag, idData) < 0) {
    opserr << "ForceBeamColumn2d::sendSelf() - failed to send ID data\n";
    return -1;
  }    

  // send the coordinate transformation
  if (crdTransf->sendSelf(commitTag, theChannel) < 0) {
    opserr << "ForceBeamColumn2d::sendSelf() - failed to send crdTrans\n";
    return -1;
  }      

  // send the beam integration
  if (beamIntegr->sendSelf(commitTag, theChannel) < 0) {
    opserr << "ForceBeamColumn2d::sendSelf() - failed to send beamIntegr\n";
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
    opserr << "ForceBeamColumn2d::sendSelf() - failed to send ID data\n";
    return -1;
  }    

  //
  // send the sections
  //
  
  for (j = 0; j<numSections; j++) {
    if (sections[j]->sendSelf(commitTag, theChannel) < 0) {
      opserr << "ForceBeamColumn2d::sendSelf() - section " << 
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

  Vector dData(1+1+NEBD+NEBD*NEBD+secDefSize+4); 
  loc = 0;

  // place double variables into Vector
  dData(loc++) = rho;
  dData(loc++) = tol;
  
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
    opserr << "ForceBeamColumn2d::sendSelf() - failed to send Vector data\n";
     return -1;
  }    

  return 0;
}    

int
ForceBeamColumn2d::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  //
  // receive the integer data containing tag, numSections and coord transformation info
  //
  int dbTag = this->getDbTag();
  int i,j,k;
  
  static ID idData(11); // one bigger than needed 

  if (theChannel.recvID(dbTag, commitTag, idData) < 0)  {
    opserr << "ForceBeamColumn2d::recvSelf() - failed to recv ID data\n";
    return -1;
  }    

  this->setTag(idData(0));
  connectedExternalNodes(0) = idData(1);
  connectedExternalNodes(1) = idData(2);
  maxIters = idData(4);
  initialFlag = idData(5);
  
  int crdTransfClassTag = idData(6);
  int crdTransfDbTag = idData(7);

  int beamIntegrClassTag = idData(8);
  int beamIntegrDbTag = idData(9);

  // create a new crdTransf object if one needed
  if (crdTransf == 0 || crdTransf->getClassTag() != crdTransfClassTag) {
      if (crdTransf != 0)
	  delete crdTransf;

      crdTransf = theBroker.getNewCrdTransf(crdTransfClassTag);

      if (crdTransf == 0) {
	opserr << "ForceBeamColumn2d::recvSelf() - failed to obtain a CrdTrans object with classTag" <<
	  crdTransfClassTag << endln;
	exit(-1);
      }
  }

  crdTransf->setDbTag(crdTransfDbTag);

  // invoke recvSelf on the coordTransf object
  if (crdTransf->recvSelf(commitTag, theChannel, theBroker) < 0)  
  {
     opserr << "ForceBeamColumn2d::sendSelf() - failed to recv crdTranf\n";
	     		     
     return -3;
  }      

  // create a new beamIntegr object if one needed
  if (beamIntegr == 0 || beamIntegr->getClassTag() != beamIntegrClassTag) {
      if (beamIntegr != 0)
	  delete beamIntegr;

      beamIntegr = theBroker.getNewBeamIntegration(beamIntegrClassTag);

      if (beamIntegr == 0) {
	opserr << "ForceBeamColumn2d::recvSelf() - failed to obtain the beam integration object with classTag" <<
	  beamIntegrClassTag << endln;
	exit(-1);
      }
  }

  beamIntegr->setDbTag(beamIntegrDbTag);

  // invoke recvSelf on the beamIntegr object
  if (beamIntegr->recvSelf(commitTag, theChannel, theBroker) < 0)  
  {
     opserr << "ForceBeamColumn2d::sendSelf() - failed to recv beam integration\n";
	     		     
     return -3;
  }      
  
  //
  // recv an ID for the sections containing each sections dbTag and classTag
  //

  ID idSections(2*idData(3));
  int loc = 0;

  if (theChannel.recvID(dbTag, commitTag, idSections) < 0)  {
    opserr << "ForceBeamColumn2d::recvSelf() - failed to recv ID data\n";
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
      opserr << "ForceBeamColumn2d::recvSelf -- failed to allocate vscommit array\n";
			      
      return -1;
    }

    // Delete the old
    if (fs != 0)
      delete [] fs;

    // Allocate the right number
    fs = new Matrix[numSections];  
    if (fs == 0) {
      opserr << "ForceBeamColumn2d::recvSelf -- failed to allocate fs array\n";
      return -1;
    }

    // Delete the old
    if (vs != 0)
      delete [] vs;

    // Allocate the right number
    vs = new Vector[numSections];  
    if (vs == 0) {
      opserr << "ForceBeamColumn2d::recvSelf -- failed to allocate vs array\n";
      return -1;
    }

    // Delete the old
    if (Ssr != 0)
      delete [] Ssr;
    
    // Allocate the right number
    Ssr = new Vector[numSections];  
    if (Ssr == 0) {
      opserr << "ForceBeamColumn2d::recvSelf -- failed to allocate Ssr array\n";
      return -1;
    }

    // create a new array to hold pointers
    sections = new SectionForceDeformation *[idData(3)];
    if (sections == 0) {
      opserr << "ForceBeamColumn2d::recvSelf() - " << 
	"out of memory creating sections array of size" << idData(3) << endln;
      exit(-1);
    }    

    loc = 0;
    
    for (i=0; i<numSections; i++) {
      int sectClassTag = idSections(loc);
      int sectDbTag = idSections(loc+1);
      loc += 2;
      sections[i] = theBroker.getNewSection(sectClassTag);
      if (sections[i] == 0) {
	opserr << "ForceBeamColumn2d::recvSelf() - " << 
	  "Broker could not create Section of class type " << sectClassTag << endln;
	exit(-1);
      }
      sections[i]->setDbTag(sectDbTag);
      if (sections[i]->recvSelf(commitTag, theChannel, theBroker) < 0) {
	opserr << "ForceBeamColumn2d::recvSelf() - section " << 
	  i << "failed to recv itself\n";
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
	  opserr << "ForceBeamColumn2d::recvSelf() - Broker could not create Section of class type" <<
	    sectClassTag << endln;
	  return -1;
	}
      }

      // recvvSelf on it
      sections[i]->setDbTag(sectDbTag);
      if (sections[i]->recvSelf(commitTag, theChannel, theBroker) < 0) {
	opserr << "ForceBeamColumn2d::recvSelf() - section " << 
	  i << "failed to recv itself\n";
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
    opserr << "ForceBeamColumn2d::sendSelf() - failed to send Vector data\n";
    return -1;
  }    
  
  loc = 0;

  // place double variables into Vector
  rho = dData(loc++);
  tol = dData(loc++);

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

  return 0;
}

int
ForceBeamColumn2d::getInitialFlexibility(Matrix &fe)
{
  fe.Zero();
  
  double L = crdTransf->getInitialLength();
  double oneOverL  = 1.0/L;  
  
  // Flexibility from elastic interior
  beamIntegr->addElasticFlexibility(L, fe);
  
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
      default:
	break;
      }
    }
  }
  
  return 0;
}

int
ForceBeamColumn2d::getInitialDeformations(Vector &v0)
{
  v0.Zero();
  if (numEleLoads < 1)
    return 0;

  double L = crdTransf->getInitialLength();
  double oneOverL  = 1.0/L;

  double xi[maxNumSections];
  beamIntegr->getSectionLocations(numSections, L, xi);
  
  double wt[maxNumSections];
  beamIntegr->getSectionWeights(numSections, L, wt);
  
  for (int i = 0; i < numSections; i++) {
    
    int order      = sections[i]->getOrder();
    const ID &code = sections[i]->getType();
    
    double xL  = xi[i];
    double xL1 = xL-1.0;
    double wtL = wt[i]*L;

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
      switch(code(ii)) {
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
      default:
	break;
      }
    }    
  }

  return 0;
}

void ForceBeamColumn2d::compSectionDisplacements(Vector sectionCoords[], Vector sectionDispls[]) const
{
   // get basic displacements and increments
   static Vector ub(NEBD);
   ub = crdTransf->getBasicTrialDisp();    

   double L = crdTransf->getInitialLength();
  
   // get integration point positions and weights
   //   const Matrix &xi_pt  = quadRule.getIntegrPointCoords(numSections);
   // get integration point positions and weights
   static double xi_pts[maxNumSections];
   beamIntegr->getSectionLocations(numSections, L, xi_pts);

   // setup Vandermode and CBDI influence matrices
   int i;
   double xi;
 
   // get CBDI influence matrix
   Matrix ls(numSections, numSections);
   getCBDIinfluenceMatrix(numSections, xi_pts, L, ls);

   // get section curvatures
   Vector kappa(numSections);  // curvature
   static Vector vs;              // section deformations 

   for (i=0; i<numSections; i++)
   {
       // THIS IS VERY INEFFICIENT ... CAN CHANGE LATER
       int sectionKey = 0;
       const ID &code = sections[i]->getType();
       int ii;
       for (ii = 0; ii < code.Size(); ii++)
	   if (code(ii) == SECTION_RESPONSE_MZ)
	   {
	       sectionKey = ii;
	       break;
	   }

       if (ii == code.Size()) {
	 opserr << "FATAL NLBeamColumn2d::compSectionDispls - section does not provide Mz response\n";
	 exit(-1);
       }
			
       // get section deformations
       vs = sections[i]->getSectionDeformation();
       kappa(i) = vs(sectionKey);
   }

   Vector w(numSections);
   static Vector xl(NDM), uxb(NDM);
   static Vector xg(NDM), uxg(NDM); 

   // w = ls * kappa;  
   w.addMatrixVector (0.0, ls, kappa, 1.0);
   
   for (i=0; i<numSections; i++)
   {
      xi = xi_pts[i];

      xl(0) = xi * L;
      xl(1) = 0;

      // get section global coordinates
      sectionCoords[i] = crdTransf->getPointGlobalCoordFromLocal(xl);

      // compute section displacements
      uxb(0) = xi * ub(0); // consider linear variation for axial displacement. CHANGE LATER!!!!!!!!!!
      uxb(1) = w(i);
             
      // get section displacements in global system 
      sectionDispls[i] = crdTransf->getPointGlobalDisplFromBasic(xi, uxb);
   }	       
   return;	       
}

void
ForceBeamColumn2d::Print(OPS_Stream &s, int flag)
{
  if (flag == 2) {

    s << "#ForceBeamColumn2D\n";

    const Vector &node1Crd = theNodes[0]->getCrds();
    const Vector &node2Crd = theNodes[1]->getCrds();	
    const Vector &node1Disp = theNodes[0]->getDisp();
    const Vector &node2Disp = theNodes[1]->getDisp();    
    
    s << "#NODE " << node1Crd(0) << " " << node1Crd(1) 
      << " " << node1Disp(0) << " " << node1Disp(1) << " " << node1Disp(2) << endln;
    
    s << "#NODE " << node2Crd(0) << " " << node2Crd(1) 
      << " " << node2Disp(0) << " " << node2Disp(1) << " " << node2Disp(2) << endln;

    double P  = Secommit(0);
    double M1 = Secommit(1);
    double M2 = Secommit(2);
    double L = crdTransf->getInitialLength();
    double V = (M1+M2)/L;
    
    double p0[3]; p0[0] = 0.0; p0[1] = 0.0; p0[2] = 0.0;
    if (numEleLoads > 0)
      this->computeReactions(p0);

    s << "#END_FORCES " << -P+p0[0] << " " << V+p0[1] << " " << M1 << endln;
    s << "#END_FORCES " << P << " " << -V+p0[2] << " " << M2 << endln;

    // plastic hinge rotation
    static Vector vp(3);
    static Matrix fe(3,3);
    this->getInitialFlexibility(fe);
    vp = crdTransf->getBasicTrialDisp();
    vp.addMatrixVector(1.0, fe, Se, -1.0);
    s << "#PLASTIC_HINGE_ROTATION " << vp[1] << " " << vp[2] << " " << 0.1*L << " " << 0.1*L << endln;
  }

  if (flag == OPS_PRINT_CURRENTSTATE) {

    s << "\nElement: " << this->getTag() << " Type: ForceBeamColumn2d ";
    s << "\tConnected Nodes: " << connectedExternalNodes ;
    s << "\tNumber of Sections: " << numSections;
    s << "\tMass density: " << rho << endln;
    beamIntegr->Print(s, flag);
    double P  = Secommit(0);
    double M1 = Secommit(1);
    double M2 = Secommit(2);
    double L = crdTransf->getInitialLength();
    double V = (M1+M2)/L;
    theVector(1) = V;
    theVector(4) = -V;
    double p0[3]; p0[0] = 0.0; p0[1] = 0.0; p0[2] = 0.0;
    if (numEleLoads > 0)
      this->computeReactions(p0);

    s << "\tEnd 1 Forces (P V M): " << -P+p0[0] << " " << V+p0[1] << " " << M1 << endln;
    s << "\tEnd 2 Forces (P V M): " << P << " " << -V+p0[2] << " " << M2 << endln;

    if (flag == 1) { 
      for (int i = 0; i < numSections; i++)
	s << "\nSection "<<i<<" :" << *sections[i];
    }
  }

  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << "\t\t\t{";
	s << "\"name\": " << this->getTag() << ", ";
	s << "\"type\": \"ForceBeamColumn2d\", ";
    s << "\"nodes\": [" << connectedExternalNodes(0) << ", " << connectedExternalNodes(1) << "], ";
    s << "\"sections\": [";
    for (int i = 0; i < numSections-1; i++)
      s << "\"" << sections[i]->getTag() << "\", ";
    s << "\"" << sections[numSections-1]->getTag() << "\"], ";
    s << "\"integration\": ";
    beamIntegr->Print(s, flag);
	s << ", \"massperlength\": " << rho << ", ";
	s << "\"crdTransformation\": \"" << crdTransf->getTag() << "\"}";
  }  
}

OPS_Stream &operator<<(OPS_Stream &s, ForceBeamColumn2d &E)
{
  E.Print(s);
  return s;
}


void
ForceBeamColumn2d::setSectionPointers(int numSec, SectionForceDeformation **secPtrs)
{
  if (numSec > maxNumSections) {
    opserr << "Error: ForceBeamColumn2d::setSectionPointers -- max number of sections exceeded";
  }
  
  numSections = numSec;
  
  if (secPtrs == 0) {
    opserr << "Error: ForceBeamColumn2d::setSectionPointers -- invalid section pointer";
  }	  
  
  sections = new SectionForceDeformation *[numSections];
  if (sections == 0) {
    opserr << "Error: ForceBeamColumn2d::setSectionPointers -- could not allocate section pointers";
  }  
  
  for (int i = 0; i < numSections; i++) {
    
    if (secPtrs[i] == 0) {
      opserr << "Error: ForceBeamColumn2d::setSectionPointers -- null section pointer " << i << endln;
    }
    
    sections[i] = secPtrs[i]->getCopy();

    if (sections[i] == 0) {
      opserr << "Error: ForceBeamColumn2d::setSectionPointers -- could not create copy of section " << i << endln;
    }
  }
  
  // allocate section flexibility matrices and section deformation vectors
  fs  = new Matrix [numSections];
  if (fs == 0) {
    opserr << "ForceBeamColumn2d::setSectionPointers -- failed to allocate fs array";
  }
  
  vs = new Vector [numSections];
  if (vs == 0) {
    opserr << "ForceBeamColumn2d::setSectionPointers -- failed to allocate vs array";
  }
  
  Ssr  = new Vector [numSections];
  if (Ssr == 0) {
    opserr << "ForceBeamColumn2d::setSectionPointers -- failed to allocate Ssr array";
  }
  
  vscommit = new Vector [numSections];
  if (vscommit == 0) {
    opserr << "ForceBeamColumn2d::setSectionPointers -- failed to allocate vscommit array";   
  }
  
}

int
ForceBeamColumn2d::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **displayModes, int numModes)
{
  static Vector v1(3);
  static Vector v2(3);

  if (displayMode >= 0) {

    theNodes[0]->getDisplayCrds(v1, fact);
    theNodes[1]->getDisplayCrds(v2, fact);

  } else {

    theNodes[0]->getDisplayCrds(v1, 0.);
    theNodes[1]->getDisplayCrds(v2, 0.);

    // add eigenvector values
    int mode = displayMode  *  -1;
    const Matrix &eigen1 = theNodes[0]->getEigenvectors();
    const Matrix &eigen2 = theNodes[1]->getEigenvectors();
    if (eigen1.noCols() >= mode) {
      for (int i = 0; i < 2; i++) {
	v1(i) += eigen1(i,mode-1)*fact;
	v2(i) += eigen2(i,mode-1)*fact;    
      }    
    } 
  }


  
  return theViewer.drawLine (v1, v2, 1.0, 1.0, this->getTag());
}

Response*
ForceBeamColumn2d::setResponse(const char **argv, int argc, OPS_Stream &output)
{
  Response *theResponse = 0;

  output.tag("ElementOutput");
  output.attr("eleType","ForceBeamColumn2d");
  output.attr("eleTag",this->getTag());
  output.attr("node1",connectedExternalNodes[0]);
  output.attr("node2",connectedExternalNodes[1]);

  // global force - 
  if (strcmp(argv[0],"forces") == 0 || strcmp(argv[0],"force") == 0
      || strcmp(argv[0],"globalForce") == 0 || strcmp(argv[0],"globalForces") == 0) {

    output.tag("ResponseType","Px_1");
    output.tag("ResponseType","Py_1");
    output.tag("ResponseType","Mz_1");
    output.tag("ResponseType","Px_2");
    output.tag("ResponseType","Py_2");
    output.tag("ResponseType","Mz_2");

    theResponse =  new ElementResponse(this, 1, theVector);
  
  // local force -
  } else if (strcmp(argv[0],"localForce") == 0 || strcmp(argv[0],"localForces") == 0) {

    output.tag("ResponseType","N_1");
    output.tag("ResponseType","V_1");
    output.tag("ResponseType","M_1");
    output.tag("ResponseType","N_2");
    output.tag("ResponseType","V_2");
    output.tag("ResponseType","M_2");

    theResponse =  new ElementResponse(this, 2, theVector);

  // basic force -
  } else if (strcmp(argv[0],"basicForce") == 0 || strcmp(argv[0],"basicForces") == 0) {

    output.tag("ResponseType","N");
    output.tag("ResponseType","M_1");
    output.tag("ResponseType","M_2");

    theResponse =  new ElementResponse(this, 7, Vector(3));

  // basic stiffness -
  } else if (strcmp(argv[0],"basicStiffness") == 0) {

    output.tag("ResponseType","N");
    output.tag("ResponseType","M1");
    output.tag("ResponseType","M2");

    theResponse =  new ElementResponse(this, 19, Matrix(3,3));

  // chord rotation -
  } else if (strcmp(argv[0],"chordRotation") == 0 || strcmp(argv[0],"chordDeformation") == 0 
	     || strcmp(argv[0],"basicDeformation") == 0) {

    output.tag("ResponseType","eps");
    output.tag("ResponseType","theta_1");
    output.tag("ResponseType","theta_2");

    theResponse =  new ElementResponse(this, 3, Vector(3));
  
  // plastic rotation -
  } else if (strcmp(argv[0],"plasticRotation") == 0 || strcmp(argv[0],"plasticDeformation") == 0) {

    output.tag("ResponseType","epsP");
    output.tag("ResponseType","thetaP_1");
    output.tag("ResponseType","thetaP_2");

    theResponse =  new ElementResponse(this, 4, Vector(3));

  // point of inflection
  } else if (strcmp(argv[0],"inflectionPoint") == 0) {
    
    output.tag("ResponseType","inflectionPoint");

    theResponse =  new ElementResponse(this, 5, 0.0);
  
  // tangent drift
  } else if (strcmp(argv[0],"tangentDrift") == 0) {
    theResponse =  new ElementResponse(this, 6, Vector(2));

  // basic forces
  } else if (strcmp(argv[0],"basicForce") == 0)
    theResponse = new ElementResponse(this, 7, Se);

  /*
  // Curvature sensitivity
  else if (strcmp(argv[0],"dcurvdh") == 0)
    return new ElementResponse(this, 7, Vector(numSections));

  // basic deformation sensitivity
  else if (strcmp(argv[0],"dvdh") == 0)
    return new ElementResponse(this, 8, Vector(3));
  */

  // plastic deformation sensitivity
  else if (strcmp(argv[0],"dvpdh") == 0)
    return new ElementResponse(this, 9, Vector(3));

  // basic force sensitivity
  else if (strcmp(argv[0],"dqdh") == 0)
    return new ElementResponse(this, 12, Vector(3));

  else if (strcmp(argv[0],"integrationPoints") == 0)
    theResponse = new ElementResponse(this, 10, Vector(numSections));

  else if (strcmp(argv[0],"integrationWeights") == 0)
    theResponse = new ElementResponse(this, 11, Vector(numSections));

  else if (strcmp(argv[0],"RayleighForces") == 0 || strcmp(argv[0],"rayleighForces") ==0 || strcmp(argv[0],"dampingForces") == 0) { 

    output.tag("ResponseType","Px_1");
    output.tag("ResponseType","Py_1");
    output.tag("ResponseType","Mz_1");
    output.tag("ResponseType","Px_2");
    output.tag("ResponseType","Py_2");
    output.tag("ResponseType","Mz_2");

    theResponse =  new ElementResponse(this, 13, theVector);
  
    // section response -
  } else if (strstr(argv[0],"sectionX") != 0) {
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
      
      if (strcmp(argv[2],"dsdh") != 0) {
	theResponse = sections[sectionNum]->setResponse(&argv[2], argc-2, output);
      } else {
	int order = sections[sectionNum]->getOrder();
	theResponse = new ElementResponse(this, 76, Vector(order));
	Information &info = theResponse->getInformation();
	info.theInt = sectionNum;
      }
    }
  }

  // section response -
  else if (strstr(argv[0],"section") != 0) {

    if (argc > 1) {

      int sectionNum = atoi(argv[1]);

      if (sectionNum > 0 && sectionNum <= numSections && argc > 2) {

	double xi[maxNumSections];
	double L = crdTransf->getInitialLength();
	beamIntegr->getSectionLocations(numSections, L, xi);

	output.tag("GaussPointOutput");
	output.attr("number",sectionNum);
	output.attr("eta",xi[sectionNum-1]*L);
	
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

	  output.endTag();

	}

	if (numResponse == 0) // no valid responses found
	  delete theCResponse;
	else
	  theResponse = theCResponse;

      }
    }
  }
  
  output.endTag(); // ElementOutput

  return theResponse;
}

int 
ForceBeamColumn2d::getResponse(int responseID, Information &eleInfo)
{
  static Vector vp(3);
  static Matrix fe(3,3);

  if (responseID == 1)
    return eleInfo.setVector(this->getResistingForce());

  else if (responseID == 13)
    return eleInfo.setVector(this->getRayleighDampingForces());
  
  else if (responseID == 2) {
    double p0[3]; p0[0] = 0.0; p0[1] = 0.0; p0[2] = 0.0;
    if (numEleLoads > 0)
      this->computeReactions(p0);
    theVector(3) =  Se(0);
    theVector(0) = -Se(0)+p0[0];
    theVector(2) = Se(1);
    theVector(5) = Se(2);
    double V = (Se(1)+Se(2))/crdTransf->getInitialLength();
    theVector(1) =  V+p0[1];
    theVector(4) = -V+p0[2];
    return eleInfo.setVector(theVector);
  }
      
  // Chord rotation
  else if (responseID == 3) {
    vp = crdTransf->getBasicTrialDisp();
    return eleInfo.setVector(vp);
  }

  // Plastic rotation
  else if (responseID == 4) {
    this->getInitialFlexibility(fe);
    vp = crdTransf->getBasicTrialDisp();
    vp.addMatrixVector(1.0, fe, Se, -1.0);
    static Vector v0(3);
    this->getInitialDeformations(v0);
    vp.addVector(1.0, v0, -1.0);
    return eleInfo.setVector(vp);
  }

  // Point of inflection
  else if (responseID == 5) {
    double LI = 0.0;

    if (fabs(Se(1)+Se(2)) > DBL_EPSILON) {
      double L = crdTransf->getInitialLength();
      
      LI = Se(1)/(Se(1)+Se(2))*L;
    }

    return eleInfo.setDouble(LI);
  }

  // Tangent drift
  else if (responseID == 6) {
    double d2 = 0.0;
    double d3 = 0.0;

    double L = crdTransf->getInitialLength();
    
    // Location of inflection point from node I
    double LI = 0.0;
    if (fabs(Se(1)+Se(2)) > DBL_EPSILON)
      LI = Se(1)/(Se(1)+Se(2))*L;
      
    double wts[maxNumSections];
    beamIntegr->getSectionWeights(numSections, L, wts);
    
    double pts[maxNumSections];
    beamIntegr->getSectionLocations(numSections, L, pts);
    
    int i;
    for (i = 0; i < numSections; i++) {
      double x = pts[i]*L;
      if (x > LI)
	continue;
      const ID &type = sections[i]->getType();
      int order = sections[i]->getOrder();
      double kappa = 0.0;
      for (int j = 0; j < order; j++)
	if (type(j) == SECTION_RESPONSE_MZ)
	  kappa += vs[i](j);
      double b = -LI+x;
      d2 += (wts[i]*L)*kappa*b;
    }
    
    d2 += beamIntegr->getTangentDriftI(L, LI, Se(1), Se(2));
    
    for (i = numSections-1; i >= 0; i--) {
      double x = pts[i]*L;
      if (x < LI)
	continue;
      const ID &type = sections[i]->getType();
      int order = sections[i]->getOrder();
      double kappa = 0.0;
      for (int j = 0; j < order; j++)
	if (type(j) == SECTION_RESPONSE_MZ)
	  kappa += vs[i](j);
      double b = x-LI;
      d3 += (wts[i]*L)*kappa*b;
    }
    
    d3 += beamIntegr->getTangentDriftJ(L, LI, Se(1), Se(2));

    static Vector d(2);
    d(0) = d2;
    d(1) = d3;

    return eleInfo.setVector(d);
  }

  else if (responseID == 7)
    return eleInfo.setVector(Se);

  else if (responseID == 19)
    return eleInfo.setMatrix(kv);


  /*
  // Curvature sensitivity
  else if (responseID == 7) {
    Vector curv(numSections);
    for (int i = 0; i < numSections; i++) {
      int order = sections[i]->getOrder();
      const ID &type = sections[i]->getType();
      const Vector &dedh = sections[i]->getdedh();
      for (int j = 0; j < order; j++) {
	if (type(j) == SECTION_RESPONSE_MZ)
	  curv(i) = dedh(j);
      }
    }
    return eleInfo.setVector(curv);
  }
  */

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

  else
    return -1;
}

int 
ForceBeamColumn2d::getResponseSensitivity(int responseID, int gradNumber,
					  Information &eleInfo)
{
  // Basic deformation sensitivity
  if (responseID == 3) {  
    const Vector &dvdh = crdTransf->getBasicDisplSensitivity(gradNumber);
    return eleInfo.setVector(dvdh);
  }

  // Basic force sensitivity
  else if (responseID == 7) {
    static Vector dqdh(3);

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
    //opserr << "FBC2d::getRespSens dspdh: " << dsdh;
    static Vector dqdh(3);

    const Vector &dvdh = crdTransf->getBasicDisplSensitivity(gradNumber);

    dqdh.addMatrixVector(0.0, kv, dvdh, 1.0);

    dqdh.addVector(1.0, this->computedqdh(gradNumber), 1.0);

    //opserr << "FBC2d::getRespSens dqdh: " << dqdh;
 
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
      default:
	break;
      }
    }

    /*
    opserr << "FBC2d::getRespSens dsdh=b*dqdh+dspdh: " << dsdh;

    dsdh.Zero();
    if (numEleLoads > 0) {
      this->computeSectionForceSensitivity(dsdh, sectionNum-1, gradNumber);
    }
    const Matrix &ks = sections[sectionNum-1]->getSectionTangent();
    const Vector &dedh  = sections[sectionNum-1]->getSectionDeformationSensitivity(gradNumber);
    dsdh.addMatrixVector(1.0, ks, dedh, 1.0);
    dsdh.addVector(1.0, sections[sectionNum-1]->getStressResultantSensitivity(gradNumber, true), 1.0);

    opserr << "FBC2d::getRespSens dsdh=b*dqdh+dspdh: " << dsdh;
    */

    return eleInfo.setVector(dsdh);
  }

  // Plastic deformation sensitivity
  else if (responseID == 4) {
    static Vector dvpdh(3);

    const Vector &dvdh = crdTransf->getBasicDisplSensitivity(gradNumber);

    dvpdh = dvdh;
    //opserr << dvpdh;

    static Matrix fe(3,3);
    this->getInitialFlexibility(fe);

    const Vector &dqdh = this->computedqdh(gradNumber);

    dvpdh.addMatrixVector(1.0, fe, dqdh, -1.0);
    //opserr << dvpdh;

    static Matrix fek(3,3);
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
ForceBeamColumn2d::setParameter(const char **argv, int argc, Parameter &param)
{
  if (argc < 1)
    return -1;

  int result = -1;

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

  // If the parameter belongs to a section or lower
  if (strstr(argv[0],"section") != 0) {
    
    if (argc < 3)
      return -1;

    // Get section number: 1...Np
    int sectionNum = atoi(argv[1]);
    
    if (sectionNum > 0 && sectionNum <= numSections)
      return sections[sectionNum-1]->setParameter(&argv[2], argc-2, param);

    else
      return -1;

    /*
    // Find the right section and call its setParameter method
    int ok = 0;
    for (int i = 0; i < numSections; i++)
      if (paramSectionTag == sections[i]->getTag())
	ok += sections[i]->setParameter(&argv[2], argc-2, param);

    return ok;
    */
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
ForceBeamColumn2d::updateParameter (int parameterID, Information &info)
{
  if (parameterID == 1) {
    rho = info.theDouble;
    return 0;
  }
  else
    return -1;
}

int
ForceBeamColumn2d::activateParameter(int passedParameterID)
{
  parameterID = passedParameterID;

  return 0;  
}

const Matrix&
ForceBeamColumn2d::getInitialStiffSensitivity(int gradNumber)
{
  theMatrix.Zero();
  return theMatrix;
}

const Matrix&
ForceBeamColumn2d::getMassSensitivity(int gradNumber)
{
  theMatrix.Zero();
  
  double L = crdTransf->getInitialLength();
  if (rho != 0.0 && parameterID == 1)
    theMatrix(0,0) = theMatrix(1,1) = theMatrix(3,3) = theMatrix(4,4) = 0.5*L;
  
  return theMatrix;
}

const Vector&
ForceBeamColumn2d::getResistingForceSensitivity(int gradNumber)
{
  static Vector dqdh(3);
  dqdh = this->computedqdh(gradNumber);

  // Transform forces
  double dp0dh[3]; dp0dh[0] = 0.0; dp0dh[1] = 0.0; dp0dh[2] = 0.0;
  this->computeReactionSensitivity(dp0dh, gradNumber);
  Vector dp0dhVec(dp0dh, 3);

  static Vector P(6);
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
ForceBeamColumn2d::commitSensitivity(int gradNumber, int numGrads)
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

  static Vector dqdh(3);
  dqdh = this->computedqdh(gradNumber);

  // dvdh = A dudh + dAdh u
  const Vector &dvdh = crdTransf->getBasicDisplSensitivity(gradNumber);
  dqdh.addMatrixVector(1.0, kv, dvdh, 1.0);  // A dudh

  if (crdTransf->isShapeSensitivity()) {
    //const Vector &dAdh_u = crdTransf->getBasicTrialDispShapeSensitivity();
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
ForceBeamColumn2d::computedqdh(int gradNumber)
{
  //opserr << "FBC2d::computedqdh " << gradNumber << endln;

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

  static Vector dvdh(3);
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
      default:
	break;
      }
    }
  }

  static Matrix dfedh(3,3);
  dfedh.Zero();

  if (beamIntegr->addElasticFlexDeriv(L, dfedh, dLdh) < 0)
    dvdh.addMatrixVector(1.0, dfedh, Se, -1.0);
  
  //opserr << "dfedh: " << dfedh << endln;

  static Vector dqdh(3);
  dqdh.addMatrixVector(0.0, kv, dvdh, 1.0);
  
  //opserr << "dqdh: " << dqdh << endln;

  return dqdh;
}

const Matrix&
ForceBeamColumn2d::computedfedh(int gradNumber)
{
  static Matrix dfedh(3,3);

  dfedh.Zero();

  double L = crdTransf->getInitialLength();
  double oneOverL  = 1.0/L;  

  double dLdh = crdTransf->getdLdh();
  double d1oLdh = crdTransf->getd1overLdh();

  beamIntegr->addElasticFlexDeriv(L, dfedh, dLdh);

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
