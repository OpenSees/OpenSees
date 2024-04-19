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

// $Revision: 1.1 $
// $Date: 2010-05-04 17:14:45 $
// $Source: /scratch/slocal/chroot/cvsroot/openseescomp/CompositePackages/mixedBeamColumn/MixedBeamColumn2d.cpp,v $

#include <MixedBeamColumn2d.h>
#include <elementAPI.h>
#include <G3Globals.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <iomanip>

#include <Information.h>
#include <MatrixUtil.h>
#include <Domain.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Renderer.h>
#include <math.h>
#include <ElementResponse.h>
#include <ElementalLoad.h>
#include <iostream>
#include <fstream>
#include <Node.h>
#include <Message.h>

#include <LobattoBeamIntegration.h>
#include <LegendreBeamIntegration.h>
#include <RadauBeamIntegration.h>
#include <NewtonCotesBeamIntegration.h>
#include <TrapezoidalBeamIntegration.h>
#include <RegularizedHingeIntegration.h>

// Constants that define the dimensionality
#define  NDM   2                      // dimension of the problem (2d)
#define  NND   3                      // number of nodal dof's
#define  NEGD  6                      // number of element global dof's
#define  NDM_SECTION  2               // number of section dof's without torsion
#define  NDM_NATURAL  3               // number of element dof's in the basic system without torsion
#define  MAX_NUM_SECTIONS  10         // maximum number of sections allowed

using namespace std;


Matrix MixedBeamColumn2d::theMatrix(NEGD,NEGD);
Matrix MixedBeamColumn2d::theNaturalMatrix(NDM_NATURAL,NDM_NATURAL);
Matrix MixedBeamColumn2d::theSectionNaturalMatrix(NDM_SECTION,NDM_NATURAL);
Matrix MixedBeamColumn2d::G(NDM_NATURAL,NDM_NATURAL);
Matrix MixedBeamColumn2d::G2(NDM_NATURAL,NDM_NATURAL);
Matrix MixedBeamColumn2d::H(NDM_NATURAL,NDM_NATURAL);
Matrix MixedBeamColumn2d::H12(NDM_NATURAL,NDM_NATURAL);
Matrix MixedBeamColumn2d::H22(NDM_NATURAL,NDM_NATURAL);
Matrix MixedBeamColumn2d::Md(NDM_NATURAL,NDM_NATURAL);
Matrix MixedBeamColumn2d::Kg(NDM_NATURAL,NDM_NATURAL);
Matrix MixedBeamColumn2d::GT(NDM_NATURAL,NDM_NATURAL);
Matrix MixedBeamColumn2d::G2T(NDM_NATURAL,NDM_NATURAL);
Matrix MixedBeamColumn2d::GMHT(NDM_NATURAL,NDM_NATURAL);
Matrix MixedBeamColumn2d::ks(NDM_SECTION,NDM_SECTION);
Matrix MixedBeamColumn2d::fs(NDM_SECTION,NDM_SECTION);


Vector MixedBeamColumn2d::theVector(NEGD);
Vector MixedBeamColumn2d::theNaturalVector(NDM_NATURAL);
Vector MixedBeamColumn2d::naturalDisp(NDM_NATURAL);
Vector MixedBeamColumn2d::naturalIncrDeltaDisp(NDM_NATURAL);
Vector MixedBeamColumn2d::V2(NDM_NATURAL);

Vector *MixedBeamColumn2d::sectionDefShapeFcn = 0;
Vector *MixedBeamColumn2d::sectionForceShapeFcn = 0;
Matrix *MixedBeamColumn2d::nldhat = 0;
Matrix *MixedBeamColumn2d::nd1 = 0;
Matrix *MixedBeamColumn2d::nd2 = 0;
Matrix *MixedBeamColumn2d::nd1T = 0;
Matrix *MixedBeamColumn2d::nd2T = 0;

// Documentation: Two Dimensional Mixed Beam Column Element
// element MixedBeamColumn2d $tag $iNode $jNode $numIntgrPts $secTag $transfTag <-mass $massDens>
//   <-integration $intType> <-doRayleigh $rFlag> <-geomNonlinear>
//
// Required Input Parameters:
//   $tag					integer tag identifying the element
//   $iNode, $jNode         end nodes
//   $numIntgrPts 			number of integration points along the element length
//   $secTag 				identifier for previously-defined section object
//   $transfTag   			identifier for previously-defined coordinate-transformation (CrdTransf) object
//
//
// Optional Input:
//   -mass $massDens
//       $massDens          element mass density (per unit length), from which a lumped-mass matrix is formed (optional, default=0.0)
//   -integration $intType
//       $intType           numerical integration type, options are Lobotto, Legendre, Radau, NewtonCotes, Trapezoidal (optional, default= Lobotto)
//   -doRayleigh $rFlag
//       $rFlag             optional, default = 1
//                              rFlag = 0 no rayleigh damping
//                              rFlag = 1 include rayleigh damping (default)
//   -geomNonlinear            perform analysis with internal geometric nonlinearity
//
//
// References:
//   1. Bulent N. Alemdar and Donald W. White, “Displacement, Flexibility, and Mixed Beam-Column Finite
//      Element Formulations for Distributed Plasticity Analysis,” Journal of Structural Engineering 131,
//      no. 12 (December 2005): 1811-1819.
//   2. Cenk Tort and Jerome F. Hajjar, “Mixed Finite Element for Three-Dimensional Nonlinear Dynamic
//      Analysis of Rectangular Concrete-Filled Steel Tube Beam-Columns,” Journal of Engineering Mechanics
//      136, no. 11 (November 0, 2010): 1329-1339.
//   3. Denavit, M. D. and Hajjar, J. F. (2010). "Nonlinear Seismic Analysis of Circular Concrete-Filled
//      Steel Tube Members and Frames," Report No. NSEL-023, Newmark Structural Laboratory Report Series
//      (ISSN 1940-9826), Department of Civil and Environmental Engineering, University of Illinois at
//      Urbana-Champaign, Urbana, Illinois, March.
//

void * OPS_MixedBeamColumn2d() {
  // Variables to retrieve input
  int iData[10];
  double dData[10];
  int sDataLength = 40;
  //char *sData  = new char[sDataLength];
  //char *sData2 = new char[sDataLength];
  int numData;

  // Check the number of dimensions
  if (OPS_GetNDM() != NDM) {
   opserr << "ERROR: MixedBeamColumn2d: invalid number of dimensions\n";
   return 0;
  }

  // Check the number of degrees of freedom
  if (OPS_GetNDF() != NND) {
   opserr << "ERROR: MixedBeamColumn2d: invalid number of degrees of freedom\n";
   return 0;
  }

  // Check for minimum number of arguments
  if (OPS_GetNumRemainingInputArgs() < 5) {
    opserr << "ERROR: MixedBeamColumn2d, too few arguments: eleTag,ndI,ndJ,transfTag,integrationTag\n";
    return 0;
  }

  numData = 5;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid element data - MixedBeamColumn2d\n";
    return 0;
  }

  int eleTag = iData[0];
  int nodeI = iData[1];
  int nodeJ = iData[2];
  int transfTag = iData[3];
  int beamIntTag = iData[4];
  
  // Get the coordinate transformation
  CrdTransf *theTransf = OPS_getCrdTransf(transfTag);
  if (theTransf == 0) {
    opserr << "WARNING geometric transformation with tag " << transfTag << "not found for element " << eleTag << endln;
    return 0;
  }

  // Get beam integrataion
  BeamIntegrationRule* theRule = OPS_getBeamIntegrationRule(beamIntTag);
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
    
  // Set Default Values for Optional Input
  int doRayleigh = 1;
  double massDens = 0.0;
  bool geomLinear = true;

  // Loop through remaining arguments to get optional input
  while ( OPS_GetNumRemainingInputArgs() > 0 ) {
    const char *sData = OPS_GetString();
    if ( strcmp(sData,"-mass") == 0 ) {
      numData = 1;
      if (OPS_GetDoubleInput(&numData, dData) != 0) {
      opserr << "WARNING invalid input, want: -mass $massDens \n";
      return 0;
      }
      massDens = dData[0];

    } else if ( strcmp(sData,"-doRayleigh") == 0 ) {
        numData = 1;
        if (OPS_GetInt(&numData, &doRayleigh) != 0) {
          opserr << "WARNING: Invalid doRayleigh in element MixedBeamColumn2d " << eleTag;
          return 0;
        }

    } else if ( strcmp(sData,"-geomNonlinear") == 0 ) {
      geomLinear = false;

    } else {
      opserr << "WARNING unknown option " << sData << "\n";
    }
  }

  // now create the element and add it to the Domain
  Element *theElement = new MixedBeamColumn2d(eleTag, nodeI, nodeJ, secTags.Size(), sections, *bi, *theTransf, massDens, doRayleigh, geomLinear);

  if (theElement == 0) {
    opserr << "WARNING ran out of memory creating element with tag " << eleTag << endln;
    return 0;
  }

  delete [] sections;
  return theElement;
}


// constructor which takes the unique element tag, sections,
// and the node ID's of it's nodal end points.
// allocates the necessary space needed by each object
MixedBeamColumn2d::MixedBeamColumn2d (int tag, int nodeI, int nodeJ, int numSec,
                                      SectionForceDeformation **sec,
                                      BeamIntegration &bi,
                                      CrdTransf &coordTransf,
                                      double massDensPerUnitLength,
                                      int damp, bool geomLin):
  Element(tag,ELE_TAG_MixedBeamColumn2d),
  connectedExternalNodes(2), beamIntegr(0), numSections(0), sections(0),
  crdTransf(0), doRayleigh(damp), geomLinear(geomLin),
  rho(massDensPerUnitLength), initialLength(0.0),
  initialFlag(0), itr(0),
  V(NDM_NATURAL), committedV(NDM_NATURAL),
  naturalForce(NDM_NATURAL), committedNaturalForce(NDM_NATURAL),
  internalForce(NDM_NATURAL), committedInternalForce(NDM_NATURAL),
  lastNaturalDisp(NDM_NATURAL), committedLastNaturalDisp(NDM_NATURAL),
  sp(0),
  Hinv(NDM_NATURAL,NDM_NATURAL), committedHinv(NDM_NATURAL,NDM_NATURAL),
  GMH(NDM_NATURAL,NDM_NATURAL), committedGMH(NDM_NATURAL,NDM_NATURAL),
  kv(NDM_NATURAL,NDM_NATURAL),
  kvcommit(NDM_NATURAL,NDM_NATURAL),
  Ki(0),
  sectionForceFibers(0), committedSectionForceFibers(0),
  sectionDefFibers(0), committedSectionDefFibers(0),
  sectionFlexibility(0), committedSectionFlexibility(0)
{
   theNodes[0] = 0;
   theNodes[1] = 0;

   connectedExternalNodes(0) = nodeI;
   connectedExternalNodes(1) = nodeJ;

   // get copy of the beam integration object
   beamIntegr = bi.getCopy();
   if(beamIntegr == 0) {
     opserr<<"Error: MixedBeamColumn2d::MixedBeamColumn2d: could not create copy of beam integration object" << endln;
     exit(-1);
   }

   // get copy of the transformation object
   crdTransf = coordTransf.getCopy2d();
   if (crdTransf == 0) {
      opserr << "Error: MixedBeamColumn2d::MixedBeamColumn2d: could not create copy of coordinate transformation object" << endln;
      exit(-1);
   }


   //this->setSectionPointers(numSec,sec);
   if (numSec > MAX_NUM_SECTIONS) {
     opserr << "Error: MixedBeamColumn2d::setSectionPointers -- max number of sections exceeded";
   }

   numSections = numSec;

   if (sec == 0) {
     opserr << "Error: MixedBeamColumn2d::setSectionPointers -- invalid section pointer";
   }

   sections = new SectionForceDeformation *[numSections];
   if (sections == 0) {
     opserr << "Error: MixedBeamColumn2d::setSectionPointers -- could not allocate section pointers";
   }

   for (int i = 0; i < numSections; i++) {

     if (sec[i] == 0) {
       opserr << "Error: MixedBeamColumn2d::setSectionPointers -- null section pointer " << i << endln;
     }

     sections[i] = (SectionForceDeformation*) sec[i]->getCopy();

     if (sections[i] == 0) {
       opserr << "Error: MixedBeamColumn2d::setSectionPointers -- could not create copy of section " << i << endln;
     }

   }

   p0[0] = 0.0;
   p0[1] = 0.0;
   p0[2] = 0.0;

   // Allocate section vectors and matrices   
   this->setSectionPointers();
   

   V.Zero();
   naturalForce.Zero();
   internalForce.Zero();
   lastNaturalDisp.Zero();
   Hinv.Zero();
   GMH.Zero();
   kv.Zero();

   committedV.Zero();
   committedNaturalForce.Zero();
   committedInternalForce.Zero();
   committedLastNaturalDisp.Zero();
   committedHinv.Zero();
   committedGMH.Zero();
   kvcommit.Zero();

   if (sectionDefShapeFcn == 0)
     sectionDefShapeFcn = new Vector [MAX_NUM_SECTIONS];
   if (sectionForceShapeFcn == 0)
     sectionForceShapeFcn = new Vector [MAX_NUM_SECTIONS];
   if (nldhat == 0)
     nldhat  = new Matrix [MAX_NUM_SECTIONS];
   if (nd1 == 0)
     nd1  = new Matrix [MAX_NUM_SECTIONS];
   if (nd2 == 0)
     nd2  = new Matrix [MAX_NUM_SECTIONS];
   if (nd1T == 0)
     nd1T  = new Matrix [MAX_NUM_SECTIONS];
   if (nd2T == 0)
     nd2T  = new Matrix [MAX_NUM_SECTIONS];
   if (!sectionDefShapeFcn || !sectionForceShapeFcn || !nldhat || !nd1 || !nd2 || !nd1T || !nd2T ) {
     opserr << "MixedBeamColumn2d::MixedBeamColumn2d() -- failed to allocate static section arrays";
     exit(-1);
   }

   int i;
   for ( i=0; i<MAX_NUM_SECTIONS; i++ ){
     nd1T[i] = Matrix(NDM_NATURAL,NDM_SECTION);
     nd2T[i] = Matrix(NDM_NATURAL,NDM_SECTION);
   }
}

void
MixedBeamColumn2d::setSectionPointers(void)
{
  if (numSections < 1)
    return;
  
  // Element vectors and matrices
  if (sectionForceFibers != 0)
    delete [] sectionForceFibers;
  sectionForceFibers = new Vector [numSections];

  if (committedSectionForceFibers != 0)
    delete [] committedSectionForceFibers;
  committedSectionForceFibers = new Vector [numSections];

  if (sectionDefFibers != 0)
    delete [] sectionDefFibers;  
  sectionDefFibers = new Vector [numSections];

  if (committedSectionDefFibers != 0)
    delete [] committedSectionDefFibers;
  committedSectionDefFibers = new Vector [numSections];

  if (sectionFlexibility != 0)
    delete [] sectionFlexibility;
  sectionFlexibility = new Matrix [numSections];

  if (committedSectionFlexibility != 0)
    delete [] committedSectionFlexibility;
  committedSectionFlexibility = new Matrix [numSections];


  for (int i = 0; i < numSections; i++){
    sectionForceFibers[i] = Vector(NDM_SECTION);
    sectionForceFibers[i].Zero();
    committedSectionForceFibers[i] = Vector(NDM_SECTION);
    committedSectionForceFibers[i].Zero();
    sectionDefFibers[i] = Vector(NDM_SECTION);
    sectionDefFibers[i].Zero();
    committedSectionDefFibers[i] = Vector(NDM_SECTION);
    committedSectionDefFibers[i].Zero();
    sectionFlexibility[i] = Matrix(NDM_SECTION,NDM_SECTION);
    sectionFlexibility[i].Zero();
    committedSectionFlexibility[i] = Matrix(NDM_SECTION,NDM_SECTION);
    committedSectionFlexibility[i].Zero();
  }
}

// constructor:
// invoked by a FEM_ObjectBroker, recvSelf() needs to be invoked on this object.
// CONSTRUCTOR FOR PARALLEL PROCESSING
MixedBeamColumn2d::MixedBeamColumn2d():
  Element(0,ELE_TAG_MixedBeamColumn2d),
  connectedExternalNodes(2), beamIntegr(0), numSections(0), sections(0), crdTransf(0), doRayleigh(0), geomLinear(true),
  rho(0.0), initialLength(0.0),
  initialFlag(0), itr(0),
  V(NDM_NATURAL), committedV(NDM_NATURAL),
  naturalForce(NDM_NATURAL), committedNaturalForce(NDM_NATURAL),
  internalForce(NDM_NATURAL), committedInternalForce(NDM_NATURAL),
  lastNaturalDisp(NDM_NATURAL), committedLastNaturalDisp(NDM_NATURAL),
  sp(0),
  Hinv(NDM_NATURAL,NDM_NATURAL), committedHinv(NDM_NATURAL,NDM_NATURAL),
  GMH(NDM_NATURAL,NDM_NATURAL), committedGMH(NDM_NATURAL,NDM_NATURAL),
  kv(NDM_NATURAL,NDM_NATURAL), kvcommit(NDM_NATURAL,NDM_NATURAL),
  Ki(0),
  sectionForceFibers(0), committedSectionForceFibers(0), sectionDefFibers(0), committedSectionDefFibers(0),
  sectionFlexibility(0), committedSectionFlexibility(0)
{
  theNodes[0] = 0;
  theNodes[1] = 0;

  p0[0] = 0.0;
  p0[1] = 0.0;
  p0[2] = 0.0;

  V.Zero();
  naturalForce.Zero();
  internalForce.Zero();
  lastNaturalDisp.Zero();
  Hinv.Zero();
  GMH.Zero();
  kv.Zero();

  committedV.Zero();
  committedNaturalForce.Zero();
  committedInternalForce.Zero();
  committedLastNaturalDisp.Zero();
  committedHinv.Zero();
  committedGMH.Zero();
  kvcommit.Zero();

  if (sectionDefShapeFcn == 0)
     sectionDefShapeFcn = new Vector [MAX_NUM_SECTIONS];
  if (sectionForceShapeFcn == 0)
     sectionForceShapeFcn = new Vector [MAX_NUM_SECTIONS];
  if (nldhat == 0)
     nldhat = new Matrix [MAX_NUM_SECTIONS];
  if (nd1 == 0)
     nd1  = new Matrix [MAX_NUM_SECTIONS];
  if (nd2 == 0)
     nd2  = new Matrix [MAX_NUM_SECTIONS];
  if (nd1T == 0)
     nd1T  = new Matrix [MAX_NUM_SECTIONS];
  if (nd2T == 0)
     nd2T  = new Matrix [MAX_NUM_SECTIONS];
  if (!sectionDefShapeFcn || !sectionForceShapeFcn || !nldhat || !nd1 || !nd2 || !nd1T || !nd2T ) {
    opserr << "MixedBeamColumn2d::MixedBeamColumn2d() -- failed to allocate static section arrays";
    exit(-1);
  }

  int i;
  for ( i=0; i<MAX_NUM_SECTIONS; i++ ){
    nd1T[i] = Matrix(NDM_NATURAL,NDM_SECTION);
    nd2T[i] = Matrix(NDM_NATURAL,NDM_SECTION);
  }

}

MixedBeamColumn2d::~MixedBeamColumn2d() {

  if (sections) {
    for (int i=0; i < numSections; i++) {
      if (sections[i]) {
        delete sections[i];
      }
    }
    delete [] sections;
  }

  if (crdTransf)
   delete crdTransf;

  if(beamIntegr != 0)
   delete beamIntegr;

  if (sp != 0)
    delete sp;

  if (Ki != 0)
   delete Ki;

  if(sectionForceFibers != 0)
   delete [] sectionForceFibers;

  if(committedSectionForceFibers != 0)
   delete [] committedSectionForceFibers;

  if(sectionDefFibers != 0)
   delete [] sectionDefFibers;

  if(committedSectionDefFibers != 0)
   delete [] committedSectionDefFibers;

  if(sectionFlexibility != 0)
   delete [] sectionFlexibility;

  if(committedSectionFlexibility != 0)
   delete [] committedSectionFlexibility;
}


int MixedBeamColumn2d::getNumExternalNodes(void) const {
   return 2;
}


const ID & MixedBeamColumn2d::getExternalNodes(void) {
   return connectedExternalNodes;
}

Node ** MixedBeamColumn2d::getNodePtrs(void) {
   return theNodes;
}

int MixedBeamColumn2d::getNumDOF(void) {
   return NEGD;
}

void MixedBeamColumn2d::setDomain(Domain *theDomain)
{
  // check Domain is not null - invoked when object removed from a domain
  if (theDomain == 0) {
    theNodes[0] = 0;
    theNodes[1] = 0;

    opserr << "MixedBeamColumn2d::setDomain:  theDomain = 0 ";
    exit(0);
  }

  // get pointers to the nodes
  int Nd1 = connectedExternalNodes(0);
  int Nd2 = connectedExternalNodes(1);

  theNodes[0] = theDomain->getNode(Nd1);
  theNodes[1] = theDomain->getNode(Nd2);

  if (theNodes[0] == 0) {
    opserr << "MixedBeamColumn2d::setDomain: Nd1: ";
    opserr << Nd1 << "does not exist in model\n";
    exit(0);
  }

  if (theNodes[1] == 0) {
    opserr << "MixedBeamColumn2d::setDomain: Nd2: ";
    opserr << Nd2 << "does not exist in model\n";
    exit(0);
  }

  // call the DomainComponent class method
  this->DomainComponent::setDomain(theDomain);

  // ensure connected nodes have correct number of dof's
  int dofNode1 = theNodes[0]->getNumberDOF();
  int dofNode2 = theNodes[1]->getNumberDOF();

  if ((dofNode1 != NND) || (dofNode2 != NND)) {
    opserr << "MixedBeamColumn2d::setDomain(): Nd2 or Nd1 incorrect dof ";
    exit(0);
  }

  // initialize the transformation
  if (crdTransf->initialize(theNodes[0], theNodes[1])) {
    opserr << "MixedBeamColumn2d::setDomain(): Error initializing coordinate transformation";
    exit(0);
  }

  // Check element length
  if (crdTransf->getInitialLength() == 0.0) {
    opserr << "MixedBeamColumn2d::setDomain(): Zero element length:" << this->getTag();
    exit(0);
  }
}

int MixedBeamColumn2d::commitState() {
  int err = 0; // error flag
  int i = 0; // integers for loops

  // call element commitState to do any base class stuff
  if ((err = this->Element::commitState()) != 0) {
   opserr << "MixedBeamColumn2d::commitState () - failed in base class";
   return err;
  }

  // commit the sections
  do {
    err = sections[i++]->commitState();
  } while (err == 0 && i < numSections);

  if (err)
    return err;

  // commit the transformation between coord. systems
  if ((err = crdTransf->commitState()) != 0)
    return err;

  // commit the element variables state
  committedV = V;
  committedNaturalForce = naturalForce;
  committedInternalForce = internalForce;
  committedLastNaturalDisp = lastNaturalDisp;
  committedHinv = Hinv;
  committedGMH = GMH;
  kvcommit = kv;
  for( i = 0; i < numSections; i++){
    committedSectionForceFibers[i] = sectionForceFibers[i];
    committedSectionDefFibers[i] = sectionDefFibers[i];
    committedSectionFlexibility[i] = sectionFlexibility[i];
  }

  // Reset iteration counter
  itr = 0;

  return err;
}


int MixedBeamColumn2d::revertToLastCommit() {
  int err;
  int i = 0;

  do {
    err = sections[i]->revertToLastCommit();
    i++;
  } while (err == 0 && i < numSections);

  if (err)
    return err;

  // revert the transformation to last commit
  if ((err = crdTransf->revertToLastCommit()) != 0)
    return err;

  // revert the element state to last commit
  V = committedV;
  internalForce = committedInternalForce;
  naturalForce = committedNaturalForce;
  lastNaturalDisp = committedLastNaturalDisp;
  Hinv = committedHinv;
  GMH = committedGMH;
  kv   = kvcommit;
  for( i = 0; i < numSections; i++){
    sectionForceFibers[i] = committedSectionForceFibers[i];
    sectionDefFibers[i] = committedSectionDefFibers[i];
    sectionFlexibility[i] = committedSectionFlexibility[i];
  }

  // Reset iteration counter
  itr = 0;

  return err;
}


int MixedBeamColumn2d::revertToStart() {
  int err;
  int i,j,k; // for loops
  i = 0;

  // revert the sections state to start
  do {
     err = sections[i++]->revertToStart();

  }while (err == 0 && i < numSections);

  if (err)
    return err;

  // revert the transformation to start
  if ((err = crdTransf->revertToStart()) != 0)
    return err;

  // revert the element state to start

  // Set initial length
  initialLength = crdTransf->getInitialLength();

  // Get the numerical integration weights
  double wt[MAX_NUM_SECTIONS]; // weights of sections or gauss points of integration points
  beamIntegr->getSectionWeights(numSections, initialLength, wt);

  // Vector of zeros to use at initial natural displacements
  theNaturalVector.Zero();

  // Set initial shape functions
  for ( i = 0; i < numSections; i++ ){
    nldhat[i] = this->getNld_hat(i, theNaturalVector, initialLength, geomLinear);
    nd1[i] = this->getNd1(i, theNaturalVector, initialLength, geomLinear);
    nd2[i] = this->getNd2(i, 0, initialLength);

    for( j = 0; j < NDM_SECTION; j++ ){
      for( k = 0; k < NDM_NATURAL; k++ ){
        nd1T[i](k,j) = nd1[i](j,k);
        nd2T[i](k,j) = nd2[i](j,k);
      }
    }
  }

  // Set initial and committed section flexibility and GJ
  for ( i = 0; i < numSections; i++ ){
    getSectionTangent(i,2,ks);
    invertMatrix(NDM_SECTION,ks,sectionFlexibility[i]);
    committedSectionFlexibility[i] = sectionFlexibility[i];
  }

  // Set initial and committed section forces and deformations
  for ( i = 0; i < numSections; i++ ){
    sectionForceFibers[i].Zero();
    committedSectionForceFibers[i].Zero();
    sectionDefFibers[i].Zero();
    committedSectionDefFibers[i].Zero();
  }

  // Compute the following matrices: G, G2, H, H12, H22, Md, Kg
  G.Zero();
  G2.Zero();
  H.Zero();
  H12.Zero();
  H22.Zero();
  Md.Zero();
  Kg.Zero();
  for( i = 0; i < numSections; i++ ){
    G   = G   + initialLength * wt[i] * nd1T[i] * nldhat[i];
    G2  = G2  + initialLength * wt[i] * nd2T[i] * nldhat[i];
    H   = H   + initialLength * wt[i] * nd1T[i] * sectionFlexibility[i] * nd1[i];
    H12 = H12 + initialLength * wt[i] * nd1T[i] * sectionFlexibility[i] * nd2[i];
    H22 = H22 + initialLength * wt[i] * nd2T[i] * sectionFlexibility[i] * nd2[i];
    // Md is zero since deformations are zero
    Kg  = Kg  + initialLength * wt[i] * this->getKg(i, 0.0, initialLength);
  }

  // Compute the inverse of the H matrix
  invertMatrix(NDM_NATURAL, H, Hinv);
  committedHinv = Hinv;

  // Compute the GMH matrix ( G + Md - H12 ) and its transpose
  GMH = G + Md - H12;
  //GMH = G; // Omit P-small delta
  committedGMH = GMH;

  // Compute the transposes of the following matrices: G2, GMH
  for( i = 0; i < NDM_NATURAL; i++ ){
    for( j = 0; j < NDM_NATURAL; j++ ){
      G2T(i,j) = G2(j,i);
      GMHT(i,j) = GMH(j,i);
    }
  }

  // Compute the stiffness matrix
  kv.Zero();
  kv = ( Kg + G2 + G2T - H22 ) + GMHT * Hinv * GMH;
  kvcommit = kv;
  Ki = new Matrix(crdTransf->getInitialGlobalStiffMatrix(kv));

  // Vector V is zero at initial state
  V.Zero();
  committedV.Zero();

  // Internal force is zero at initial state
  internalForce.Zero();
  committedInternalForce.Zero();
  naturalForce.Zero();
  committedNaturalForce.Zero();

  // Last natural displacement is zero at initial state
  lastNaturalDisp.Zero();
  committedLastNaturalDisp.Zero();

  // Reset iteration counter
  itr = 0;

  // Set initialFlag to 1 so update doesn't call again
  initialFlag = 1;

  return err;
}

const Matrix & MixedBeamColumn2d::getInitialStiff(void) {
  // If things haven't be initialized, then do so
  if (initialFlag == 0) {
    this->revertToStart();
  }
  return *Ki;
}

const Matrix & MixedBeamColumn2d::getTangentStiff(void) {
  // If things haven't be initialized, then do so
  if (initialFlag == 0) {
    this->revertToStart();
  }
  crdTransf->update();  // Will remove once we clean up the corotational 3d transformation -- MHS
  return crdTransf->getGlobalStiffMatrix(kv,internalForce);
}

const Vector & MixedBeamColumn2d::getResistingForce(void) {
  crdTransf->update();  // Will remove once we clean up the corotational 3d transformation -- MHS
  Vector p0Vec(p0, NDM_NATURAL);
  return crdTransf->getGlobalResistingForce(internalForce, p0Vec);
}

int MixedBeamColumn2d::update() {

  // If things haven't be initialized, then do so
  if (initialFlag == 0) {
    this->revertToStart();
  }

  int i,j,k; // integers for loops

  // Update iteration counter
  // says how many times update has been called since the last commit state
  itr++;

  // Update Coordinate Transformation
  crdTransf->update();

  // Current Length
  double currentLength;
  if (geomLinear) {
    currentLength = initialLength;
  } else {
    currentLength = crdTransf->getDeformedLength();
  }

  // Compute the natural displacements
  naturalDisp.Zero();
  naturalDisp = crdTransf->getBasicTrialDisp();

  naturalIncrDeltaDisp.Zero();
  naturalIncrDeltaDisp = naturalDisp - lastNaturalDisp;
  lastNaturalDisp = naturalDisp;

  // Get the numerical integration weights
  double wt[MAX_NUM_SECTIONS]; // weights of sections or gauss points of integration points
  beamIntegr->getSectionWeights(numSections, initialLength, wt);

  // Compute shape functions and their transposes
  for ( i = 0; i < numSections; i++ ){
    nldhat[i] = this->getNld_hat(i, naturalDisp, currentLength, geomLinear);
    sectionDefShapeFcn[i] = this->getd_hat(i, naturalDisp, currentLength, geomLinear);
    nd1[i] = this->getNd1(i, naturalDisp, currentLength, geomLinear);
    if (geomLinear) {
      nd2[i].Zero();
    } else {
      nd2[i] = this->getNd2(i, internalForce(0), currentLength);
    }

    for( j = 0; j < NDM_SECTION; j++ ){
      for( k = 0; k < NDM_NATURAL; k++ ){
        nd1T[i](k,j) = nd1[i](j,k);
        nd2T[i](k,j) = nd2[i](j,k);
      }
    }
  }

  // Update natural force
  if (geomLinear) {
    naturalForce = naturalForce + Hinv * ( GMH * naturalIncrDeltaDisp + V );
  } else {
    naturalForce = naturalForce + Hinv * ( GMH * naturalIncrDeltaDisp + V );
  }

  // Update sections
  for ( i = 0; i < numSections; i++){
    // Compute section deformations
    sectionForceShapeFcn[i] = nd1[i] * naturalForce;
    if (sp != 0) {
      const Matrix &s_p = *sp;
      for ( j = 0; j < NDM_SECTION; j++ ) {
        sectionForceShapeFcn[i](j) += s_p(j,i);
      }
    }
    sectionDefFibers[i] = sectionDefFibers[i] + sectionFlexibility[i] * ( sectionForceShapeFcn[i] - sectionForceFibers[i] );

    // Send section deformation to section object
    setSectionDeformation(i,sectionDefFibers[i]);

    // Get section force vector
    getSectionStress(i,sectionForceFibers[i]);

    // Get section tangent matrix
    getSectionTangent(i,1,ks);

    // Compute section flexibility matrix
    invertMatrix(NDM_SECTION,ks,sectionFlexibility[i]);

  }

  // Compute the following matrices: V, V2, G, G2, H, H12, H22, Md, Kg
  V.Zero();
  V2.Zero();
  G.Zero();
  G2.Zero();
  H.Zero();
  H12.Zero();
  H22.Zero();
  Md.Zero();
  Kg.Zero();

  for( i = 0; i < numSections; i++ ){
    V   = V   + initialLength * wt[i] * nd1T[i] * (sectionDefShapeFcn[i] - sectionDefFibers[i] - sectionFlexibility[i] * ( sectionForceShapeFcn[i] - sectionForceFibers[i] ) );
    V2  = V2  + initialLength * wt[i] * nd2T[i] * (sectionDefShapeFcn[i] - sectionDefFibers[i]);
    G   = G   + initialLength * wt[i] * nd1T[i] * nldhat[i];
    G2  = G2  + initialLength * wt[i] * nd2T[i] * nldhat[i];
    H   = H   + initialLength * wt[i] * nd1T[i] * sectionFlexibility[i] * nd1[i];
    H12 = H12 + initialLength * wt[i] * nd1T[i] * sectionFlexibility[i] * nd2[i];
    H22 = H22 + initialLength * wt[i] * nd2T[i] * sectionFlexibility[i] * nd2[i];
    if (!geomLinear) {
      Kg = Kg  + initialLength * wt[i] * this->getKg(i, sectionForceFibers[i](0), currentLength);
        // sectionForceFibers[i](0) is the axial load, P
      Md = Md  + initialLength * wt[i] * this->getMd(i, sectionDefShapeFcn[i], sectionDefFibers[i], currentLength);
    }
  }

  // Compute the inverse of the H matrix
  invertMatrix(NDM_NATURAL, H, Hinv);

  // Compute the GMH matrix ( G + Md - H12 ) and its transpose
  GMH = G + Md - H12;
  //GMH = G; // Omit P-small delta

  // Compute the transposes of the following matrices: G, G2, GMH
  for( i = 0; i < NDM_NATURAL; i++ ) {
    for( j = 0; j < NDM_NATURAL; j++ ) {
      GT(i,j) = G(j,i);
      G2T(i,j) = G2(j,i);
      GMHT(i,j) = GMH(j,i);
    }
  }


  // Define the internal force
  if (geomLinear) {
    internalForce = GT * naturalForce + V2 + GMHT * Hinv * V;
  } else {
    internalForce = GT * naturalForce + V2 + GMHT * Hinv * V;
  }


  // Compute the stiffness matrix without the torsion term
  if (geomLinear) {
    kv = ( Kg + G2 + G2T - H22 ) + GMHT * Hinv * GMH;
  } else {
    kv = ( Kg + G2 + G2T - H22 ) + GMHT * Hinv * GMH;
  }

  return 0;
}

const Matrix & MixedBeamColumn2d::getMass(void) {
  theMatrix.Zero();

  if (rho != 0.0) {
    theMatrix(0,0) = theMatrix(1,1) =
    theMatrix(3,3) = theMatrix(4,4) = 0.5*initialLength*rho;
  }

  return theMatrix;
}

const Matrix & MixedBeamColumn2d::getDamp(void) {
  theMatrix.Zero();

  // Add the damping forces
  if ( doRayleigh == 1 && (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0) ) {
    theMatrix = this->Element::getDamp();
  }

  return theMatrix;
}

void MixedBeamColumn2d::zeroLoad(void) {
  if (sp != 0)
    sp->Zero();

  p0[0] = 0.0;
  p0[1] = 0.0;
  p0[2] = 0.0;
}

int MixedBeamColumn2d::addLoad(ElementalLoad *theLoad, double loadFactor) {

  int type;
  const Vector &data = theLoad->getData(type, loadFactor);

  if (sp == 0) {
    sp = new Matrix(NDM_SECTION,numSections);
    if (sp == 0) {
      opserr << "MixedBeamColumn2d::addLoad -- out of memory\n";
      exit(-1);
    }
  }

  double L = crdTransf->getInitialLength();

  double xi[MAX_NUM_SECTIONS];
  beamIntegr->getSectionLocations(numSections, L, xi);

  if (type == LOAD_TAG_Beam2dUniformLoad) {
    double wy = data(0)*loadFactor;  // Transverse
    double wx = data(1)*loadFactor;  // Axial

    Matrix &s_p = *sp;

    // Accumulate applied section forces due to element loads
    for (int i = 0; i < numSections; i++) {
      double x = xi[i]*L;
      // Axial
      s_p(0,i) += wx*(L-x);
      // Moment
      s_p(1,i) += wy*0.5*x*(x-L);
    }

    // Accumulate reactions in basic system
    p0[0] -= wx*L;
    double V;
    V = 0.5*wy*L;
    p0[1] -= V;
    p0[2] -= V;

  } else if (type == LOAD_TAG_Beam3dPointLoad) {
    double Py = data(0)*loadFactor;
    double N  = data(1)*loadFactor;
    double aOverL = data(2);

    if (aOverL < 0.0 || aOverL > 1.0)
      return 0;

    double a = aOverL*L;

    double Vy2 = Py*aOverL;
    double Vy1 = Py-Vy2;

    Matrix &s_p = *sp;

    // Accumulate applied section forces due to element loads
    for (int i = 0; i < numSections; i++) {
      double x = xi[i]*L;
      if (x <= a) {
        s_p(0,i) += N;
        s_p(1,i) -= x*Vy1;
      }
      else {
        s_p(1,i) -= (L-x)*Vy2;
      }
    }

    // Accumulate reactions in basic system
    p0[0] -= N;
    p0[1] -= Vy1;
    p0[2] -= Vy2;

  } else {
    opserr << "MixedBeamColumn2d::addLoad() -- load type unknown for element with tag: " <<
        this->getTag() << endln;

    return -1;
  }

  return 0;
}

const Vector & MixedBeamColumn2d::getResistingForceIncInertia() {

  // Compute the current resisting force
  theVector = this->getResistingForce();

  // Add the inertial forces
  if (rho != 0.0) {
    const Vector &accel1 = theNodes[0]->getTrialAccel();
    const Vector &accel2 = theNodes[1]->getTrialAccel();

    double L = crdTransf->getInitialLength();
    double m = 0.5*rho*L;

    theVector(0) += m*accel1(0);
    theVector(1) += m*accel1(1);
    theVector(4) += m*accel2(0);
    theVector(5) += m*accel2(1);
  }

  // Add the damping forces
  if ( doRayleigh == 1 && (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0) ) {
    theVector += this->getRayleighDampingForces();
  }

  return theVector;
}

void MixedBeamColumn2d::Print(OPS_Stream &s, int flag) {

  if (flag == 1) {
    s << "\nElement: " << this->getTag() << " Type: MixedBeamColumn2d ";
    s << "\tConnected Nodes: " << connectedExternalNodes ;
    s << "\tNumber of Sections: " << numSections;
    s << "\tMass density: " << rho;
    for (int i = 0; i < numSections; i++)
       s << "\nSection "<<i<<" :" << *sections[i];

  } else if (flag == 33) {
    s << "\nElement: " << this->getTag() << " Type: MixedBeamColumn2d ";
    double xi[MAX_NUM_SECTIONS]; // location of sections or gauss points or integration points
    beamIntegr->getSectionLocations(numSections, initialLength, xi);
    double wt[MAX_NUM_SECTIONS]; // weights of sections or gauss points of integration points
    beamIntegr->getSectionWeights(numSections, initialLength, wt);
    s << "\n section xi wt";
    for (int i = 0; i < numSections; i++)
      s << "\n"<<i<<" "<<xi[i]<<" "<<wt[i];

  } else if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << "\t\t\t{";
    s << "\"name\": " << this->getTag() << ", ";
    s << "\"type\": \"MixedBeamColumn2d\", ";
    s << "\"nodes\": [" << connectedExternalNodes(0) << ", " << connectedExternalNodes(1) << "], ";
    s << "\"sections\": [";
    for (int i = 0; i < numSections - 1; i++)
      s << "\"" << sections[i]->getTag() << "\", ";
    s << "\"" << sections[numSections - 1]->getTag() << "\"], ";
    s << "\"integration\": ";
    beamIntegr->Print(s, flag);
    s << ", \"massperlength\": " << rho << ", ";
    s << "\"crdTransformation\": \"" << crdTransf->getTag() << "\"";
    if (!doRayleigh)
      s << ", \"doRayleigh\": false";
    if (geomLinear)
      s << ", \"geomLinear\": true";
    s << "}";

  } else {
    s << "\nElement: " << this->getTag() << " Type: MixedBeamColumn2d ";
    s << "\tConnected Nodes: " << connectedExternalNodes ;
    s << "\tNumber of Sections: " << numSections;
    s << "\tMass density: " << rho << endln;
  }
}

OPS_Stream &operator<<(OPS_Stream &s, MixedBeamColumn2d &E) {
    E.Print(s);
    return s;
}


Response* MixedBeamColumn2d::setResponse(const char **argv, int argc,
                                         OPS_Stream &output) {

  Response *theResponse = 0;

  output.tag("ElementOutput");
  output.attr("eleType","MixedBeamColumn2d");
  output.attr("eleTag",this->getTag());
  output.attr("node1",connectedExternalNodes[0]);
  output.attr("node2",connectedExternalNodes[1]);

  //
  // we compare argv[0] for known response types
  //

  // global force -
  if (strcmp(argv[0],"forces") == 0 ||
      strcmp(argv[0],"force") == 0 ||
      strcmp(argv[0],"globalForce") == 0 ||
      strcmp(argv[0],"globalForces") == 0) {

    output.tag("ResponseType","Px_1");
    output.tag("ResponseType","Py_1");
    output.tag("ResponseType","Mz_1");
    output.tag("ResponseType","Px_2");
    output.tag("ResponseType","Py_2");
    output.tag("ResponseType","Mz_2");

    theResponse = new ElementResponse(this, 1, theVector);

  // local force
  }  else if (strcmp(argv[0],"localForce") == 0 ||
              strcmp(argv[0],"localForces") == 0) {

    output.tag("ResponseType","N_1");
    output.tag("ResponseType","V_1");
    output.tag("ResponseType","M_1");
    output.tag("ResponseType","N_2");
    output.tag("ResponseType","V_2");
    output.tag("ResponseType","M_2");

    theResponse = new ElementResponse(this, 2, theVector);

  // basic or natural forces
  } else if (strcmp(argv[0],"basicForce") == 0 ||
             strcmp(argv[0],"basicForces") == 0) {

    output.tag("ResponseType","N");
    output.tag("ResponseType","M_1");
    output.tag("ResponseType","M_2");

    theResponse = new ElementResponse(this, 3, Vector(3));

  } else if (strcmp(argv[0],"sectionDeformation_Force") == 0) {

    int i;
    char *q  = new char[80];
    for ( i = 0; i < numSections; i++ ){
      sprintf(q,"axialStrain_%d",i+1);
      output.tag("ResponseType",q);
      sprintf(q,"curvature_%d",i+1);
      output.tag("ResponseType",q);
    }
    delete [] q;

    theResponse =  new ElementResponse(this, 4, Vector(2*numSections));

  } else if (strcmp(argv[0],"plasticSectionDeformation_Force") == 0) {

    int i;
    char *q  = new char[80];
    for ( i = 0; i < numSections; i++ ){
      sprintf(q,"plasticAxialStrain_%d",i+1);
      output.tag("ResponseType",q);
      sprintf(q,"plasticCurvature_%d",i+1);
      output.tag("ResponseType",q);
    }
    delete [] q;

    theResponse =  new ElementResponse(this, 5, Vector(2*numSections));

  } else if (strcmp(argv[0],"sectionStiffness") == 0) {

    int i;
    char *q  = new char[80];
    for ( i = 0; i < numSections; i++ ){
      sprintf(q,"sectionStiffness_EA_%d",i+1);
      output.tag("ResponseType",q);
      sprintf(q,"sectionStiffness_EI_%d",i+1);
      output.tag("ResponseType",q);
    }
    delete [] q;

    theResponse =  new ElementResponse(this, 6, Vector(2*numSections));

  } else if (strcmp(argv[0],"basicStiffness") == 0) {

    output.tag("ResponseType","N");
    output.tag("ResponseType","M1");
    output.tag("ResponseType","M2");

    theResponse =  new ElementResponse(this, 19, Matrix(3,3));
    
  } else if (strcmp(argv[0],"integrationPoints") == 0) {
    theResponse =  new ElementResponse(this, 100, Vector(numSections));

  } else if (strcmp(argv[0],"integrationWeights") == 0) {
    theResponse =  new ElementResponse(this, 101, Vector(numSections));

  } else if (strcmp(argv[0],"sectionTags") == 0) {
    theResponse = new ElementResponse(this, 110, ID(numSections));
    
  } else if (strcmp(argv[0],"connectedNodes") == 0) {
    theResponse =  new ElementResponse(this, 102, ID(2));

  } else if (strcmp(argv[0],"numSections") == 0 ||
             strcmp(argv[0],"numberOfSections") == 0 ) {
    theResponse =  new ElementResponse(this, 103, ID(1));
  }

  else if (strcmp(argv[0],"sectionX") ==0) {
    if (argc > 2) {
      
      float sectionLoc = atof(argv[1]);
      
      double xi[MAX_NUM_SECTIONS];
      double L = crdTransf->getInitialLength();
      beamIntegr->getSectionLocations(numSections, L, xi);
      
      sectionLoc /= L;
      
      float minDistance = fabs(xi[0]-sectionLoc);
      int sectionNum = 1;
      for (int i = 1; i < numSections; i++) {
	if (fabs(xi[i]-sectionLoc) < minDistance) {
	  minDistance = fabs(xi[i]-sectionLoc);
	  sectionNum = i+1;
	}
      }
      
      output.tag("GaussPointOutput");
      output.attr("number",sectionNum);
      output.attr("eta",xi[sectionNum-1]*L);
      
      // A kinda hacky thing to record section shear even though this
      // element doesn't include shear effects
      bool thisSectionHasShear = false;
      int order = sections[sectionNum-1]->getOrder();
      const ID &type = sections[sectionNum-1]->getType();
      for (int i = 0; i < order; i++) {
	if (type(i) == SECTION_RESPONSE_VY)
	  thisSectionHasShear = true;
      }
      
      if (!thisSectionHasShear || strcmp(argv[2],"force") != 0)
	theResponse =  sections[sectionNum-1]->setResponse(&argv[2], argc-2, output);	  
      else
	theResponse = new ElementResponse(this, 500 + sectionNum, Vector(order));
      
      output.endTag();
    }
  }
  
  else if (strcmp(argv[0],"section") ==0) {
    if (argc > 2) {

      int sectionNum = atoi(argv[1]);
      if (sectionNum > 0 && sectionNum <= numSections) {

        double xi[MAX_NUM_SECTIONS];
        double L = crdTransf->getInitialLength();
        beamIntegr->getSectionLocations(numSections, L, xi);

        output.tag("GaussPointOutput");
        output.attr("number",sectionNum);
        output.attr("eta",xi[sectionNum-1]*L);

	// A kinda hacky thing to record section shear even though this
	// element doesn't include shear effects
	bool thisSectionHasShear = false;
	int order = sections[sectionNum-1]->getOrder();
	const ID &type = sections[sectionNum-1]->getType();
	for (int i = 0; i < order; i++) {
	  if (type(i) == SECTION_RESPONSE_VY)
	    thisSectionHasShear = true;
	}

	if (!thisSectionHasShear || strcmp(argv[2],"force") != 0)
	  theResponse =  sections[sectionNum-1]->setResponse(&argv[2], argc-2, output);	  
	else
	  theResponse = new ElementResponse(this, 500 + sectionNum, Vector(order));
	
        output.endTag();
      }
    }
  }

  if (theResponse == 0)
    theResponse = crdTransf->setResponse(argv, argc, output);
  
  output.endTag();
  return theResponse;
}


int MixedBeamColumn2d::getResponse(int responseID, Information &eleInfo) {
  if (responseID == 1) { // global forces
    return eleInfo.setVector(this->getResistingForce());

  } else if (responseID == 2) { // local forces
    // Axial
    double N = internalForce(0);
    theVector(3) =  N;
    theVector(0) = -N+p0[0];

    // Moments about z and shears along y
    double M1 = internalForce(1);
    double M2 = internalForce(2);
    theVector(2)  = M1;
    theVector(5) = M2;
    double L = crdTransf->getInitialLength();
    double V = (M1+M2)/L;
    theVector(1) =  V+p0[1];
    theVector(4) = -V+p0[2];

    return eleInfo.setVector(theVector);

  } else if (responseID == 3) { // basic forces
    return eleInfo.setVector(internalForce);

  } else if (responseID == 4) { // section deformation (from forces)

    int i;
    Vector tempVector(2*numSections);
    tempVector.Zero();
    for ( i = 0; i < numSections; i++ ){
      tempVector(2*i)   = sectionDefFibers[i](0);
      tempVector(2*i+1) = sectionDefFibers[i](1);
    }

    return eleInfo.setVector(tempVector);

  } else if (responseID == 5) { // plastic section deformation (from forces)

    int i;
    Vector tempVector(2*numSections);
    Vector sectionForce(NDM_SECTION);
    Vector plasticSectionDef(NDM_SECTION);
    tempVector.Zero();
    for ( i = 0; i < numSections; i++ ){

      getSectionStress(i,sectionForce);
      getSectionTangent(i,2,ks);
      invertMatrix(NDM_SECTION,ks,fs);

      plasticSectionDef = sectionDefFibers[i] - fs*sectionForce;

      tempVector(2*i)   = plasticSectionDef(0);
      tempVector(2*i+1) = plasticSectionDef(1);
    }

    return eleInfo.setVector(tempVector);

  } else if (responseID == 6) { // section stiffness (EA and EI)

    int i;
    Vector tempVector(2*numSections);
    tempVector.Zero();
    for ( i = 0; i < numSections; i++ ){
      getSectionTangent(i,1,ks);
      tempVector(2*i)   = ks(0,0); // EA
      tempVector(2*i+1) = ks(1,1); // EI
    }

    return eleInfo.setVector(tempVector);

  } else if (responseID == 19)
    return eleInfo.setMatrix(kv);
  
  else if (responseID == 100) { // integration points

    double L = crdTransf->getInitialLength();
    double pts[MAX_NUM_SECTIONS];
    beamIntegr->getSectionLocations(numSections, L, pts);
    Vector locs(numSections);
    for (int i = 0; i < numSections; i++)
      locs(i) = pts[i]*L;
    return eleInfo.setVector(locs);

  } else if (responseID == 101) { // integration weights
      double L = crdTransf->getInitialLength();
      double wts[MAX_NUM_SECTIONS];
      beamIntegr->getSectionWeights(numSections, L, wts);
      Vector weights(numSections);
      for (int i = 0; i < numSections; i++)
        weights(i) = wts[i]*L;
      return eleInfo.setVector(weights);

  } else if (responseID == 110) {
    ID tags(numSections);
    for (int i = 0; i < numSections; i++)
      tags(i) = sections[i]->getTag();
    return eleInfo.setID(tags);
      
  } else if (responseID == 102) { // connected nodes
    ID tempVector(2);
    tempVector(0) = connectedExternalNodes(0);
    tempVector(1) = connectedExternalNodes(1);
    return eleInfo.setID(tempVector);

  } else if (responseID == 103) { // number of sections
    ID tempVector(1);
    tempVector(0) = numSections;
    return eleInfo.setID(tempVector);

  } else if (responseID > 500 && responseID <= 550) {
    int sectionNum = responseID % 500;
    double L = crdTransf->getInitialLength();
    
    double M1 = internalForce(1);
    double M2 = internalForce(2);
    double Vy = (M1+M2)/L;

    int order = sections[sectionNum-1]->getOrder();
    Vector s(sections[sectionNum-1]->getStressResultant());
    const ID &type = sections[sectionNum-1]->getType();
    for (int i = 0; i < order; i++) {
      if (type(i) == SECTION_RESPONSE_VY) 
	s(i) = Vy;
    }
    return eleInfo.setVector(s);
  }

  else {
    return -1;

  }
}

Vector
MixedBeamColumn2d::getd_hat(int sec, const Vector &v, double L, bool geomLinear){
  double xi[MAX_NUM_SECTIONS];
  beamIntegr->getSectionLocations(numSections, L, xi);

  double x, C, E, F;
  Vector D_hat(NDM_SECTION);
  D_hat.Zero();

  x = L*xi[sec];
  C = 1/L;
  E = -4/L + 6*x/(L*L);
  F = -2/L + 6*x/(L*L);

  if (geomLinear) {

    D_hat(0) =  C*v(0);
    D_hat(1) =  E*v(1) + F*v(2);

  } else {

    double A,B;
    A = 1 - 4*(x/L) + 3*pow(x/L,2);
    B =   - 2*(x/L) + 3*pow(x/L,2);

    D_hat(0) =  C*v(0) +
                0.5*A*A*v(1)*v(1) +
                    A*B*v(1)*v(2) +
                0.5*B*B*v(2)*v(2);
    D_hat(1) =  E*v(1) + F*v(2);

  }

  return D_hat;
}

Matrix
MixedBeamColumn2d::getKg(int sec, double P, double L) {
  double xi[MAX_NUM_SECTIONS];
  beamIntegr->getSectionLocations(numSections, L, xi);

  double x, A, B;
  x = L*xi[sec];
  A = 1 - 4*(x/L) + 3*pow(x/L,2);
  B =   - 2*(x/L) + 3*pow(x/L,2);

  theNaturalMatrix.Zero();
  theNaturalMatrix(0,0) = P / ( L * L );
  theNaturalMatrix(1,1) = P*A*A;
  theNaturalMatrix(2,2) = P*B*B;
  theNaturalMatrix(1,2) = P*A*B;
  theNaturalMatrix(2,1) = P*A*B;

  return theNaturalMatrix;
}

Matrix MixedBeamColumn2d::getMd(int sec, Vector dShapeFcn, Vector dFibers, double L) {
  double xi[MAX_NUM_SECTIONS];
  beamIntegr->getSectionLocations(numSections, L, xi);

  double x, A, B;
  x = L*xi[sec];
  A =  ( x/L - 2*pow(x/L,2) + pow(x/L,3) )*L;
  B =          (-pow(x/L,2) + pow(x/L,3) )*L;

  theNaturalMatrix.Zero();
  theNaturalMatrix(0,1) = A * ( dShapeFcn(1) - dFibers(1) );
  theNaturalMatrix(0,2) = B * ( dShapeFcn(1) - dFibers(1) );

  return theNaturalMatrix;
}

Matrix MixedBeamColumn2d::getNld_hat(int sec, const Vector &v, double L, bool geomLinear) {
  double xi[MAX_NUM_SECTIONS];
  beamIntegr->getSectionLocations(numSections, L, xi);

  double x, C, E, F;

  x = L*xi[sec];
  C =  1/L;
  E = -4/L + 6*x/(L*L);
  F = -2/L + 6*x/(L*L);

  theSectionNaturalMatrix.Zero();
  if (geomLinear) {
    theSectionNaturalMatrix(0,0) = C;
    theSectionNaturalMatrix(1,1) = E;
    theSectionNaturalMatrix(1,2) = F;
  } else {
    double A,B;
    A = 1 - 4 * ( x / L ) + 3 * pow ( ( x / L ) , 2 );
    B = - 2 * ( x / L ) + 3 * pow ( ( x / L ) , 2 );
    theSectionNaturalMatrix(0,0) = C + C*C*v(0);
    theSectionNaturalMatrix(0,1) = A*A*v(1) + A*B*v(2);
    theSectionNaturalMatrix(0,2) = A*B*v(1) + B*B*v(2);
    theSectionNaturalMatrix(1,1) = E;
    theSectionNaturalMatrix(1,2) = F;
  }

  return theSectionNaturalMatrix;
}

Matrix
MixedBeamColumn2d::getNd2(int sec, double P, double L){
   double xi[MAX_NUM_SECTIONS];
   beamIntegr->getSectionLocations(numSections, L, xi);

   double x, A, B;

   x = L * xi[sec];
   A = L * ( x / L - 2 * pow( x / L, 2 ) + pow( x / L, 3 ) );
   B = L * ( -pow( x / L, 2 ) + pow( x / L, 3 ) );

   theSectionNaturalMatrix.Zero();
   theSectionNaturalMatrix(1,1) = P * A;
   theSectionNaturalMatrix(1,2) = P * B;

   return theSectionNaturalMatrix;
}


Matrix
MixedBeamColumn2d::getNd1(int sec, const Vector &v, double L, bool geomLinear){
   double xi[MAX_NUM_SECTIONS];
   beamIntegr->getSectionLocations(numSections, L, xi);

   double x = L * xi[sec];

   theSectionNaturalMatrix.Zero();
   if (geomLinear) {
     theSectionNaturalMatrix(0,0) = 1.0;
     theSectionNaturalMatrix(1,1) = -x/L + 1.0;
     theSectionNaturalMatrix(1,2) =  x/L;
   } else {
     double A;
     A = L * ( x/L - 2*pow(x/L,2) + pow(x/L,3) ) * v[1]
              + L * ( -pow(x/L,2) + pow(x/L,3) ) * v[2];

     theSectionNaturalMatrix(0,0) = 1.0;
     theSectionNaturalMatrix(1,0) = A;
     theSectionNaturalMatrix(1,1) = -x/L + 1.0;
     theSectionNaturalMatrix(1,2) =  x/L;
   }

   return theSectionNaturalMatrix;
}


void MixedBeamColumn2d::getSectionTangent(int sec,int type,Matrix &kSection) {
  int order = sections[sec]->getOrder();
  const ID &code = sections[sec]->getType();

  // Initialize formulation friendly variables
  kSection.Zero();

  // Get the stress resultant from section
  Matrix sectionTangent(order,order);
  if ( type == 1 ) {
    sectionTangent = sections[sec]->getSectionTangent();
  } else if ( type == 2 ) {
    sectionTangent = sections[sec]->getInitialTangent();
  } else {
    sectionTangent.Zero();
  }

  // Set Components of Section Tangent
  int i,j;
  for (i = 0; i < order; i++) {
    for (j = 0; j < order; j++) {
      switch(code(i)) {
        case SECTION_RESPONSE_P:
          switch(code(j)) {
            case SECTION_RESPONSE_P:
              kSection(0,0) = sectionTangent(i,j);
              break;
            case SECTION_RESPONSE_MZ:
              kSection(0,1) = sectionTangent(i,j);
              break;
            default:
              break;
          }
          break;
        case SECTION_RESPONSE_MZ:
          switch(code(j)) {
            case SECTION_RESPONSE_P:
              kSection(1,0) = sectionTangent(i,j);
              break;
            case SECTION_RESPONSE_MZ:
              kSection(1,1) = sectionTangent(i,j);
              break;
            default:
              break;
          }
          break;
        default:
          break;
      }
    }
  }
}

void MixedBeamColumn2d::getSectionStress(int sec,Vector &fSection) {
  int order = sections[sec]->getOrder();
  const ID &code = sections[sec]->getType();

  // Get the stress resultant from section
  Vector stressResultant = sections[sec]->getStressResultant();

  // Initialize formulation friendly variables
  fSection.Zero();

  // Set Components of Section Stress Resultant
  int j;
  for (j = 0; j < order; j++) {
    switch(code(j)) {
      case SECTION_RESPONSE_P:
        fSection(0) = stressResultant(j);
        break;
      case SECTION_RESPONSE_MZ:
        fSection(1) = stressResultant(j);
        break;
      default:
        break;
    }
  }
}

void
MixedBeamColumn2d::setSectionDeformation(int sec,Vector &defSection) {
  int order = sections[sec]->getOrder();
  const ID &code = sections[sec]->getType();

  // Initialize Section Deformation Vector
  Vector sectionDeformation(order);
  sectionDeformation.Zero();

  // Set Components of Section Deformations
  int j;
  for (j = 0; j < order; j++) {
    switch(code(j)) {
      case SECTION_RESPONSE_P:
        sectionDeformation(j) = defSection(0);
        break;
      case SECTION_RESPONSE_MZ:
        sectionDeformation(j) = defSection(1);
        break;
      default:
        break;
    }
  }

  // Set the section deformations
  int res = sections[sec]->setTrialSectionDeformation(sectionDeformation);
}


int MixedBeamColumn2d::sendSelf(int commitTag, Channel &theChannel)
{
  int dbTag = this->getDbTag();

  static ID idData(11); // Make sure cannot be 2*numSections
  idData(0) = this->getTag();
  idData(1) = connectedExternalNodes(0);
  idData(2) = connectedExternalNodes(1);  
  idData(3) = numSections;
  idData(4) = crdTransf->getClassTag();
  int crdTransfDbTag  = crdTransf->getDbTag();
  if (crdTransfDbTag  == 0) {
    crdTransfDbTag = theChannel.getDbTag();
    if (crdTransfDbTag  != 0) 
      crdTransf->setDbTag(crdTransfDbTag);
  }
  idData(5) = crdTransfDbTag;
  idData(6) = beamIntegr->getClassTag();
  int beamIntDbTag  = beamIntegr->getDbTag();
  if (beamIntDbTag  == 0) {
    beamIntDbTag = theChannel.getDbTag();
    if (beamIntDbTag  != 0) 
      beamIntegr->setDbTag(beamIntDbTag);
  }
  idData(7) = beamIntDbTag;
  idData(8) = geomLinear ? 1 : 0;
  idData(9) = initialFlag;
  idData(10) = doRayleigh;

  if (theChannel.sendID(dbTag, commitTag, idData) < 0) {
    opserr << "MixedBeamColumn2d::sendSelf() - failed to send data ID" << endln;
    return -1;
  }
  
  static Vector data(2);
  data(0) = rho;
  data(1) = initialLength;

  if (theChannel.sendVector(dbTag, commitTag, data) < 0) {
    opserr << "MixedBeamColumn2d::sendSelf() - failed to send data Vector" << endln;
    return -2;
  }

  // send the coordinate transformation
  if (crdTransf->sendSelf(commitTag, theChannel) < 0) {
    opserr << "MixedBeamColumn2d::sendSelf() - failed to send crdTranf" << endln;
    return -3;
  }      

  // send the beam integration
  if (beamIntegr->sendSelf(commitTag, theChannel) < 0) {
    opserr << "MixedBeamColumn2d::sendSelf() - failed to send beamInt" << endln;
    return -4;
  }

  //
  // send an ID for the sections containing each sections dbTag and classTag
  // if section ha no dbTag get one and assign it
  //
  ID idSections(2*numSections);
  int loc = 0;
  for (int i = 0; i<numSections; i++) {
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
    opserr << "MixedBeamColumn2d::sendSelf() - failed to send ID data" << endln;
    return -5;
  }    

  //
  // send the sections
  //
  
  for (int i = 0; i<numSections; i++) {
    if (sections[i]->sendSelf(commitTag, theChannel) < 0) {
      opserr << "MixedBeamColumn2d::sendSelf() - section " << 
	i << " failed to send itself" << endln;
      return -6;
    }
  }

  /*
  Vector committedV; // NDM_NATURAL
  Vector committedInternalForce;  // NDM_NATURAL
  Vector committedNaturalForce; // NDM_NATURAL
  Vector committedLastNaturalDisp; // NDM_NATURAL
  Matrix committedHinv; // NDM_NATURAL x NDM_NATURAL
  Matrix committedGMH; // NDM_NATURAL x NDM_NATURAL
  Matrix kvcommit; // NDM_NATURAL x NDM_NATURAL
  */
  
  // 4*NDM_NATURAL + 3*NDM_NATURAL**2 = 4*3 + 3*3**2 = 12 + 27 = 39
  int lenElementData = 4*NDM_NATURAL + 3*NDM_NATURAL*NDM_NATURAL;
  static Vector elementData(lenElementData);
  for (int i = 0; i < NDM_NATURAL; i++) {
    elementData(i              ) = committedV(i);
    elementData(i+  NDM_NATURAL) = committedInternalForce(i);
    elementData(i+2*NDM_NATURAL) = committedNaturalForce(i);
    elementData(i+3*NDM_NATURAL) = committedLastNaturalDisp(i);    
  }
  loc = 4*NDM_NATURAL;
  for (int i = 0; i < NDM_NATURAL; i++) {
    for (int j = 0; j < NDM_NATURAL; j++)
      elementData(loc++) = committedHinv(i,j);
  }
  loc = 4*NDM_NATURAL + NDM_NATURAL*NDM_NATURAL;
  for (int i = 0; i < NDM_NATURAL; i++) {
    for (int j = 0; j < NDM_NATURAL; j++)
      elementData(loc++) = committedGMH(i,j);
  }
  loc = 4*NDM_NATURAL + 2*NDM_NATURAL*NDM_NATURAL;
  for (int i = 0; i < NDM_NATURAL; i++) {
    for (int j = 0; j < NDM_NATURAL; j++)
      elementData(loc++) = kvcommit(i,j);
  }

  if (theChannel.sendVector(dbTag, commitTag, elementData) < 0) {
    opserr << "MixedBeamColumn2d::sendSelf() - failed to send elementData Vector" << endln;
    return -7;
  }
 
  /*
  Vector *committedSectionForceFibers; // numSections
  Vector *committedSectionDefFibers; // numSections
  Matrix *committedSectionFlexibility; // numSections
  */
  
  // 2*numSections*order + numSections*order**2
  int order = sections[0]->getOrder(); // Assume all sections have same order
  order = NDM_SECTION;
  int lenSectionData = 2*numSections*order + numSections*order*order;
  if (lenSectionData == lenElementData) {
    lenSectionData++;
  }
  Vector sectionData(lenSectionData);
  for (int i = 0; i < numSections; i++) {
    for (int j = 0; j < order; j++) {
      sectionData(                    i*order + j) = committedSectionForceFibers[i](j);
      sectionData(numSections*order + i*order + j) = committedSectionDefFibers[i](j);
    }
  }
  loc = 2*numSections*order;
  for (int i = 0; i < numSections; i++) {
    for (int j = 0; j < order; j++)
      for (int k = 0; k < order; k++) {
	sectionData(loc++) = committedSectionFlexibility[i](j,k);
      }
  }

  if (theChannel.sendVector(dbTag, commitTag, sectionData) < 0) {
    opserr << "MixedBeamColumn2d::sendSelf() - failed to send sectionData Vector" << endln;
    return -8;
  }  

  return 0;
}

int MixedBeamColumn2d::recvSelf(int commitTag, Channel &theChannel,
                                FEM_ObjectBroker &theBroker)
{
  int dbTag = this->getDbTag();

  static ID idData(11);
  if (theChannel.recvID(dbTag, commitTag, idData) < 0) {
    opserr << "MixedBeamColumn2d::recvSelf() - failed to receive data ID" << endln;
    return -1;
  }

  this->setTag(idData(0));
  connectedExternalNodes(0) = idData(1);
  connectedExternalNodes(1) = idData(2);  
  int nSect = idData(3);
  int crdTransfClassTag = idData(4);
  int crdTransfDbTag = idData(5);
  int beamIntClassTag = idData(6);
  int beamIntDbTag = idData(7);  
  geomLinear = (idData(8) == 1) ? true : false;
  initialFlag = idData(9);
  doRayleigh = idData(10);

  
  static Vector data(2);
  if (theChannel.recvVector(dbTag, commitTag, data) < 0) {
    opserr << "MixedBeamColumn2d::recvSelf() - failed to receive data Vector" << endln;
    return -2;
  }

  rho = data(0);
  initialLength = data(1);

  // create a new crdTransf object if one needed
  if (crdTransf == 0 || crdTransf->getClassTag() != crdTransfClassTag) {
      if (crdTransf != 0)
	  delete crdTransf;

      crdTransf = theBroker.getNewCrdTransf(crdTransfClassTag);

      if (crdTransf == 0) {
	opserr << "MixedBeamColumn2d::recvSelf() - failed to obtain a CrdTransf object with classTag " <<
	  crdTransfClassTag << endln;
	exit(-1);
      }
  }
  crdTransf->setDbTag(crdTransfDbTag);

  // invoke recvSelf on the crdTransf object
  if (crdTransf->recvSelf(commitTag, theChannel, theBroker) < 0) {
    opserr << "MixedBeamColumn2d::sendSelf() - failed to recv crdTranf" << endln;
    return -3;
  }      

  // create a new beamInt object if one needed
  if (beamIntegr == 0 || beamIntegr->getClassTag() != beamIntClassTag) {
      if (beamIntegr != 0)
	  delete beamIntegr;

      beamIntegr = theBroker.getNewBeamIntegration(beamIntClassTag);

      if (beamIntegr == 0) {
	opserr << "MixedBeamColumn2d::recvSelf() - failed to obtain the beam integration object with classTag" <<
	  beamIntClassTag << endln;
	exit(-1);
      }
  }

  beamIntegr->setDbTag(beamIntDbTag);

  // invoke recvSelf on the beamInt object
  if (beamIntegr->recvSelf(commitTag, theChannel, theBroker) < 0)  
  {
    opserr << "MixedBeamColumn2d::sendSelf() - failed to recv beam integration" << endln;
    return -4;
  }      


  ID idSections(2*nSect);

  if (theChannel.recvID(dbTag, commitTag, idSections) < 0)  {
    opserr << "DispBeamColumn2d::recvSelf() - failed to recv ID data\n";
    return -5;
  }    

  //
  // now receive the sections
  //
  
  if (numSections != nSect) {

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

    // create a new array to hold pointers
    sections = new SectionForceDeformation *[nSect];
    if (sections == 0) {
      opserr << "MixedBeamColumn2d::recvSelf() - out of memory creating sections array of size " <<
	nSect << endln;
      return -5;
    }    

    // create a section and recvSelf on it
    numSections = nSect;
    int loc = 0;
    for (int i=0; i<numSections; i++) {
      int sectClassTag = idSections(loc);
      int sectDbTag = idSections(loc+1);
      loc += 2;
      sections[i] = theBroker.getNewSection(sectClassTag);
      if (sections[i] == 0) {
	opserr << "MixedpBeamColumn2d::recvSelf() - Broker could not create Section of class type " <<
	  sectClassTag << endln;
	exit(-1);
      }
      sections[i]->setDbTag(sectDbTag);
      if (sections[i]->recvSelf(commitTag, theChannel, theBroker) < 0) {
	opserr << "MixedBeamColumn2d::recvSelf() - section " << i << " failed to recv itself" << endln;
	return -5;
      }     
    }

  } else {

    // 
    // for each existing section, check it is of correct type
    // (if not delete old & create a new one) then recvSelf on it
    //
    
    int loc = 0;
    for (int i=0; i<numSections; i++) {
      int sectClassTag = idSections(loc);
      int sectDbTag = idSections(loc+1);
      loc += 2;

      // check of correct type
      if (sections[i]->getClassTag() !=  sectClassTag) {
	// delete the old section[i] and create a new one
	delete sections[i];
	sections[i] = theBroker.getNewSection(sectClassTag);
	if (sections[i] == 0) {
	  opserr << "MixedBeamColumn2d::recvSelf() - Broker could not create Section of class type " <<
	    sectClassTag << endln;
	  exit(-1);
	}
      }

      // recvSelf on it
      sections[i]->setDbTag(sectDbTag);
      if (sections[i]->recvSelf(commitTag, theChannel, theBroker) < 0) {
	opserr << "MixedBeamColumn2d::recvSelf() - section " << i << " failed to recv itself" << endln;
	return -5;
      }     
    }
  }

  
  int lenElementData = 4*NDM_NATURAL + 3*NDM_NATURAL*NDM_NATURAL;
  static Vector elementData(lenElementData);
  if (theChannel.recvVector(dbTag, commitTag, elementData) < 0) {
    opserr << "MixedBeamColumn2d::recvSelf() - failed to receive elementData Vector" << endln;
    return -6;
  }
  
  for (int i = 0; i < NDM_NATURAL; i++) {
    committedV(i)               = elementData(i              );
    committedInternalForce(i)   = elementData(i+  NDM_NATURAL);
    committedNaturalForce(i)    = elementData(i+2*NDM_NATURAL);
    committedLastNaturalDisp(i) = elementData(i+3*NDM_NATURAL);
  }
  int loc;
  loc = 4*NDM_NATURAL;
  for (int i = 0; i < NDM_NATURAL; i++) {
    for (int j = 0; j < NDM_NATURAL; j++)
      committedHinv(i,j) = elementData(loc++);
  }
  loc = 4*NDM_NATURAL + NDM_NATURAL*NDM_NATURAL;
  for (int i = 0; i < NDM_NATURAL; i++) {
    for (int j = 0; j < NDM_NATURAL; j++)
      committedGMH(i,j) = elementData(loc++);
  }
  loc = 4*NDM_NATURAL + 2*NDM_NATURAL*NDM_NATURAL;
  for (int i = 0; i < NDM_NATURAL; i++) {
    for (int j = 0; j < NDM_NATURAL; j++)
      kvcommit(i,j) = elementData(loc++);
  }  


  // Allocate section vectors and matrices
  this->setSectionPointers();
  
  int order = sections[0]->getOrder(); // Assume all sections have same order
  order = NDM_SECTION;
  int lenSectionData = 2*numSections*order + numSections*order*order;
  if (lenSectionData == lenElementData) {
    lenSectionData++;
  }
  Vector sectionData(lenSectionData);
  if (theChannel.recvVector(dbTag, commitTag, sectionData) < 0) {
    opserr << "MixedBeamColumn2d::recvSelf() - failed to receive sectionData Vector" << endln;
    return -7;
  }

  for (int i = 0; i < numSections; i++) {
    for (int j = 0; j < order; j++) {
      committedSectionForceFibers[i](j) = sectionData(                    i*order + j);
      committedSectionDefFibers[i](j)   = sectionData(numSections*order + i*order + j);
    }
  }
  loc = 2*numSections*order;
  for (int i = 0; i < numSections; i++) {
    for (int j = 0; j < order; j++)
      for (int k = 0; k < order; k++) {
	committedSectionFlexibility[i](j,k) = sectionData(loc++);
      }
  }  

  return 0;
}


int MixedBeamColumn2d::displaySelf(Renderer& theViewer, int displayMode, float fact, const char** modes, int numMode)
{
    static Vector v1(3);
    static Vector v2(3);

    theNodes[0]->getDisplayCrds(v1, fact, displayMode);
    theNodes[1]->getDisplayCrds(v2, fact, displayMode);

    return theViewer.drawLine(v1, v2, 1.0, 1.0, this->getTag());
}
