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
// $Source: /scratch/slocal/chroot/cvsroot/openseescomp/CompositePackages/mixedBeamColumn/MixedBeamColumn3d.cpp,v $

#include <MixedBeamColumn3d.h>
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
#define  NDM   3                       // dimension of the problem (3d)
#define  NND   6                       // number of nodal dof's
#define  NEGD  12                      // number of element global dof's
#define  NDM_SECTION  3                // number of section dof's without torsion
#define  NDM_NATURAL  5                // number of element dof's in the basic system without torsion
#define  NDM_NATURAL_WITH_TORSION  6   // number of element dof's in the basic system with torsion
#define  MAX_NUM_SECTIONS  10          // maximum number of sections allowed

using namespace std;


Matrix MixedBeamColumn3d::theMatrix(NEGD,NEGD);
Vector MixedBeamColumn3d::theVector(NEGD);
double MixedBeamColumn3d::workArea[400];
Matrix MixedBeamColumn3d::transformNaturalCoords(NDM_NATURAL_WITH_TORSION,NDM_NATURAL_WITH_TORSION);
Matrix MixedBeamColumn3d::transformNaturalCoordsT(NDM_NATURAL_WITH_TORSION,NDM_NATURAL_WITH_TORSION);

Vector *MixedBeamColumn3d::sectionDefShapeFcn = 0;
Vector* MixedBeamColumn3d::sectionForceShapeFcn = 0;
Matrix *MixedBeamColumn3d::nldhat = 0;
Matrix *MixedBeamColumn3d::nd1 = 0;
Matrix *MixedBeamColumn3d::nd2 = 0;
Matrix *MixedBeamColumn3d::nd1T = 0;
Matrix *MixedBeamColumn3d::nd2T = 0;

// Documentation: Three Dimensional Mixed Beam Column Element
// element MixedBeamColumn3d $tag $iNode $jNode $numIntgrPts $secTag $transfTag <-mass $massDens>
//   <-integration $intType> <-doRayleigh $rFlag> <-geomNonlinear>
//
// Required Input Parameters:
//   $tag                   integer tag identifying the element
//   $iNode, $jNode         end nodes
//   $numIntgrPts           number of integration points along the element length
//   $secTag                identifier for previously-defined section object
//   $transfTag             identifier for previously-defined coordinate-transformation (CrdTransf) object
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
//   1. Bulent N. Alemdar and Donald W. White, "Displacement, Flexibility, and Mixed Beam-Column Finite
//      Element Formulations for Distributed Plasticity Analysis," Journal of Structural Engineering 131,
//      no. 12 (December 2005): 1811-1819.
//   2. Cenk Tort and Jerome F. Hajjar, "Mixed Finite Element for Three-Dimensional Nonlinear Dynamic
//      Analysis of Rectangular Concrete-Filled Steel Tube Beam-Columns," Journal of Engineering Mechanics
//      136, no. 11 (November 0, 2010): 1329-1339.
//   3. Denavit, M. D. and Hajjar, J. F. (2010). "Nonlinear Seismic Analysis of Circular Concrete-Filled
//      Steel Tube Members and Frames," Report No. NSEL-023, Newmark Structural Laboratory Report Series
//      (ISSN 1940-9826), Department of Civil and Environmental Engineering, University of Illinois at
//      Urbana-Champaign, Urbana, Illinois, March.
//

void * OPS_MixedBeamColumn3d() {
  // Variables to retrieve input
  int iData[10];
  double dData[10];
  int sDataLength = 40;
  //char *sData  = new char[sDataLength];
  //char *sData2 = new char[sDataLength];
  int numData;

  // Check the number of dimensions
  if (OPS_GetNDM() != NDM) {
     opserr << "ERROR: MixedBeamColumn3d: invalid number of dimensions\n";
   return 0;
  }

  // Check the number of degrees of freedom
  if (OPS_GetNDF() != NND) {
     opserr << "ERROR: MixedBeamColumn3d: invalid number of degrees of freedom\n";
   return 0;
  }

  // Check for minimum number of arguments
  if (OPS_GetNumRemainingInputArgs() < 5) {
    opserr << "ERROR: MixedBeamColumn3d, too few arguments: eleTag,ndI,ndJ,transfTag,integrationTag\n";
    return 0;
  }

  // Get required input data
  numData = 5;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid element data - MixedBeamColumn3d\n";
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
          opserr << "WARNING: Invalid doRayleigh in element MixedBeamColumn3d " << eleTag;
          return 0;
        }

    } else if ( strcmp(sData,"-geomNonlinear") == 0 ) {
      geomLinear = false;

    } else {
      opserr << "WARNING unknown option " << sData << "\n";
    }
  }

  // now create the element and add it to the Domain
  Element *theElement = new MixedBeamColumn3d(eleTag, nodeI, nodeJ, secTags.Size(), sections, *bi, *theTransf, massDens, doRayleigh, geomLinear);

  if (theElement == 0) {
    opserr << "WARNING ran out of memory creating element with tag " << eleTag << endln;
    return 0;
  }

  return theElement;
}


// constructor which takes the unique element tag, sections,
// and the node ID's of it's nodal end points.
// allocates the necessary space needed by each object
MixedBeamColumn3d::MixedBeamColumn3d (int tag, int nodeI, int nodeJ, int numSec,
                                      SectionForceDeformation **sec,
                                      BeamIntegration &bi,
                                      CrdTransf &coordTransf,
                                      double massDensPerUnitLength,
                                      int damp, bool geomLin):
  Element(tag,ELE_TAG_MixedBeamColumn3d),
  connectedExternalNodes(2), beamIntegr(0), numSections(0), sections(0),
  crdTransf(0), doRayleigh(damp), geomLinear(geomLin),
  rho(massDensPerUnitLength), initialLength(0.0),
  itr(0), initialFlag(0),
  V(NDM_NATURAL), committedV(NDM_NATURAL),
  internalForceOpenSees(NDM_NATURAL_WITH_TORSION),
  committedInternalForceOpenSees(NDM_NATURAL_WITH_TORSION),
  naturalForce(NDM_NATURAL), commitedNaturalForce(NDM_NATURAL),
  lastNaturalDisp(NDM_NATURAL), commitedLastNaturalDisp(NDM_NATURAL),
  sp(0),
  Hinv(NDM_NATURAL,NDM_NATURAL), commitedHinv(NDM_NATURAL,NDM_NATURAL),
  GMH(NDM_NATURAL,NDM_NATURAL), commitedGMH(NDM_NATURAL,NDM_NATURAL),
  kv(NDM_NATURAL_WITH_TORSION,NDM_NATURAL_WITH_TORSION),
  kvcommit(NDM_NATURAL_WITH_TORSION,NDM_NATURAL_WITH_TORSION),
  Ki(0),
  sectionForceFibers(0), commitedSectionForceFibers(0),
  sectionDefFibers(0), commitedSectionDefFibers(0),
  sectionFlexibility(0), commitedSectionFlexibility(0)
{
  theNodes[0] = 0;
  theNodes[1] = 0;

  connectedExternalNodes(0) = nodeI;
  connectedExternalNodes(1) = nodeJ;

  // get copy of the beam integration object
  beamIntegr = bi.getCopy();
  if (beamIntegr == 0) {
    opserr<<"Error: MixedBeamColumn3d::MixedBeamColumn3d: could not create copy of beam integration object" << endln;
    exit(-1);
  }

  // get copy of the transformation object
  crdTransf = coordTransf.getCopy3d();
  if (crdTransf == 0) {
    opserr << "Error: MixedBeamColumn3d::MixedBeamColumn3d: could not create copy of coordinate transformation object" << endln;
    exit(-1);
  }


  //this->setSectionPointers(numSec,sec);
  if (numSec > MAX_NUM_SECTIONS) {
    opserr << "Error: MixedBeamColumn3d::setSectionPointers -- max number of sections exceeded";
  }

  numSections = numSec;

  if (sec == 0) {
    opserr << "Error: MixedBeamColumn3d::setSectionPointers -- invalid section pointer";
  }

  sections = new SectionForceDeformation *[numSections];
  if (sections == 0) {
    opserr << "Error: MixedBeamColumn3d::setSectionPointers -- could not allocate section pointers";
  }

  for (int i = 0; i < numSections; i++) {
    if (sec[i] == 0) {
      opserr << "Error: MixedBeamColumn3d::setSectionPointers -- null section pointer " << i << endln;
    }

    sections[i] = (SectionForceDeformation*) sec[i]->getCopy();

    if (sections[i] == 0) {
      opserr << "Error: MixedBeamColumn3d::setSectionPointers -- could not create copy of section " << i << endln;
    }
  }

  p0[0] = 0.0;
  p0[1] = 0.0;
  p0[2] = 0.0;
  p0[3] = 0.0;
  p0[4] = 0.0;

  // Element vectors and matrices
  sectionForceFibers = new Vector [numSections];
  commitedSectionForceFibers = new Vector [numSections];
  sectionDefFibers = new Vector [numSections];
  commitedSectionDefFibers = new Vector [numSections];
  sectionFlexibility = new Matrix [numSections];
  commitedSectionFlexibility = new Matrix [numSections];

  for (int i = 0; i < numSections; i++){
    sectionForceFibers[i] = Vector(NDM_SECTION);
    sectionForceFibers[i].Zero();
    commitedSectionForceFibers[i] = Vector(NDM_SECTION);
    commitedSectionForceFibers[i].Zero();
    sectionDefFibers[i] = Vector(NDM_SECTION);
    sectionDefFibers[i].Zero();
    commitedSectionDefFibers[i] = Vector(NDM_SECTION);
    commitedSectionDefFibers[i].Zero();
    sectionFlexibility[i] = Matrix(NDM_SECTION,NDM_SECTION);
    sectionFlexibility[i].Zero();
    commitedSectionFlexibility[i] = Matrix(NDM_SECTION,NDM_SECTION);
    commitedSectionFlexibility[i].Zero();
  }

  V.Zero();
  internalForceOpenSees.Zero();
  naturalForce.Zero();
  lastNaturalDisp.Zero();
  Hinv.Zero();
  GMH.Zero();
  kv.Zero();

  committedV.Zero();
  committedInternalForceOpenSees.Zero();
  commitedNaturalForce.Zero();
  commitedLastNaturalDisp.Zero();
  commitedHinv.Zero();
  commitedGMH.Zero();
  kvcommit.Zero();

  if ( transformNaturalCoords(1,1) != 1 ) {
    // if transformNaturalCoords hasn't been set yet then set it
    // This is needed because the way the state determination algorithm was formulated
    // and OpenSees assume a different order of natural degrees of freedom
    // Formulation --> { e, theta_zi, theta_yi, theta_zj, theta_yj, twist }
    // OpenSees    --> { e, theta_zi, theta_zj, theta_yi, theta_yj, twist }
    transformNaturalCoords.Zero();
    transformNaturalCoords(0,0) = 1;
    transformNaturalCoords(1,1) = 1;
    transformNaturalCoords(2,3) = 1;
    transformNaturalCoords(3,2) = 1;
    transformNaturalCoords(4,4) = 1;
    transformNaturalCoords(5,5) = 1;
    transformNaturalCoordsT.Zero();
    transformNaturalCoordsT(0,0) = 1;
    transformNaturalCoordsT(1,1) = 1;
    transformNaturalCoordsT(3,2) = 1;
    transformNaturalCoordsT(2,3) = 1;
    transformNaturalCoordsT(4,4) = 1;
    transformNaturalCoordsT(5,5) = 1;
  }

  if (sectionDefShapeFcn == 0)
    sectionDefShapeFcn  = new Vector [MAX_NUM_SECTIONS];
  if (sectionForceShapeFcn == 0)
      sectionForceShapeFcn = new Vector[MAX_NUM_SECTIONS];
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
    opserr << "MixedBeamColumn3d::MixedBeamColumn3d() -- failed to allocate static section arrays";
    exit(-1);
  }

  int i;
  for ( i=0; i< MAX_NUM_SECTIONS; i++ ){
    nd1T[i] = Matrix(NDM_NATURAL,NDM_SECTION);
    nd2T[i] = Matrix(NDM_NATURAL,NDM_SECTION);
  }

}

// constructor:
// invoked by a FEM_ObjectBroker, recvSelf() needs to be invoked on this object.
// CONSTRUCTOR FOR PARALLEL PROCESSING
MixedBeamColumn3d::MixedBeamColumn3d():
  Element(0,ELE_TAG_MixedBeamColumn3d),
  connectedExternalNodes(2), beamIntegr(0), numSections(0), sections(0), crdTransf(0), doRayleigh(0), geomLinear(true),
  rho(0.0), initialLength(0.0),
  itr(0), initialFlag(0),
  V(NDM_NATURAL), committedV(NDM_NATURAL),
  internalForceOpenSees(NDM_NATURAL_WITH_TORSION), committedInternalForceOpenSees(NDM_NATURAL_WITH_TORSION),
  naturalForce(NDM_NATURAL), commitedNaturalForce(NDM_NATURAL),
  lastNaturalDisp(NDM_NATURAL), commitedLastNaturalDisp(NDM_NATURAL),
  sp(0),
  Hinv(NDM_NATURAL,NDM_NATURAL), commitedHinv(NDM_NATURAL,NDM_NATURAL),
  GMH(NDM_NATURAL,NDM_NATURAL), commitedGMH(NDM_NATURAL,NDM_NATURAL),
  kv(NDM_NATURAL_WITH_TORSION,NDM_NATURAL_WITH_TORSION), kvcommit(NDM_NATURAL_WITH_TORSION,NDM_NATURAL_WITH_TORSION),
  Ki(0),
  sectionForceFibers(0), commitedSectionForceFibers(0), sectionDefFibers(0), commitedSectionDefFibers(0),
  sectionFlexibility(0), commitedSectionFlexibility(0)
{
  theNodes[0] = 0;
  theNodes[1] = 0;

  p0[0] = 0.0;
  p0[1] = 0.0;
  p0[2] = 0.0;
  p0[3] = 0.0;
  p0[4] = 0.0;

  // Element vectors and matrices
  sectionForceFibers = new Vector [numSections];
  commitedSectionForceFibers = new Vector [numSections];
  sectionDefFibers = new Vector [numSections];
  commitedSectionDefFibers = new Vector [numSections];
  sectionFlexibility = new Matrix [numSections];
  commitedSectionFlexibility = new Matrix [numSections];

  for (int i = 0; i < numSections; i++){
    sectionForceFibers[i] = Vector(NDM_SECTION);
    sectionForceFibers[i].Zero();
    commitedSectionForceFibers[i] = Vector(NDM_SECTION);
    commitedSectionForceFibers[i].Zero();
    sectionDefFibers[i] = Vector(NDM_SECTION);
    sectionDefFibers[i].Zero();
    commitedSectionDefFibers[i] = Vector(NDM_SECTION);
    commitedSectionDefFibers[i].Zero();
    sectionFlexibility[i] = Matrix(NDM_SECTION,NDM_SECTION);
    sectionFlexibility[i].Zero();
    commitedSectionFlexibility[i] = Matrix(NDM_SECTION,NDM_SECTION);
    commitedSectionFlexibility[i].Zero();
  }

  V.Zero();
  internalForceOpenSees.Zero();
  naturalForce.Zero();
  lastNaturalDisp.Zero();
  Hinv.Zero();
  GMH.Zero();
  kv.Zero();

  committedV.Zero();
  committedInternalForceOpenSees.Zero();
  commitedNaturalForce.Zero();
  commitedLastNaturalDisp.Zero();
  commitedHinv.Zero();
  commitedGMH.Zero();
  kvcommit.Zero();


  if ( transformNaturalCoords(1,1) != 1 ) {
    // if transformNaturalCoords hasn't been set yet then set it
    transformNaturalCoords.Zero();
    transformNaturalCoords(0,0) = 1;
    transformNaturalCoords(1,1) = 1;
    transformNaturalCoords(2,3) = -1;
    transformNaturalCoords(3,2) = 1;
    transformNaturalCoords(4,4) = -1;
    transformNaturalCoords(5,5) = 1;
    transformNaturalCoordsT.Zero();
    transformNaturalCoordsT(0,0) = 1;
    transformNaturalCoordsT(1,1) = 1;
    transformNaturalCoordsT(3,2) = -1;
    transformNaturalCoordsT(2,3) = 1;
    transformNaturalCoordsT(4,4) = -1;
    transformNaturalCoordsT(5,5) = 1;
  }

  if (sectionDefShapeFcn == 0)
    sectionDefShapeFcn  = new Vector [MAX_NUM_SECTIONS];
  if (sectionForceShapeFcn == 0)
      sectionForceShapeFcn = new Vector[MAX_NUM_SECTIONS];
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
    opserr << "MixedBeamColumn3d::MixedBeamColumn3d() -- failed to allocate static section arrays";
    exit(-1);
  }

  int i;
  for ( i=0; i< MAX_NUM_SECTIONS; i++ ){
    nd1T[i] = Matrix(NDM_NATURAL,NDM_SECTION);
    nd2T[i] = Matrix(NDM_NATURAL,NDM_SECTION);
  }

}

MixedBeamColumn3d::~MixedBeamColumn3d() {

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

  if (beamIntegr != 0)
    delete beamIntegr;

  if (sp != 0)
    delete sp;

  if (Ki != 0)
    delete Ki;

  if (sectionForceFibers != 0)
    delete [] sectionForceFibers;

  if (commitedSectionForceFibers != 0)
    delete [] commitedSectionForceFibers;

  if (sectionDefFibers != 0)
    delete [] sectionDefFibers;

  if (commitedSectionDefFibers != 0)
    delete [] commitedSectionDefFibers;

  if (sectionFlexibility != 0)
    delete [] sectionFlexibility;

  if (commitedSectionFlexibility != 0)
    delete [] commitedSectionFlexibility;
}

int MixedBeamColumn3d::getNumExternalNodes(void) const {
   return 2;
}

const ID & MixedBeamColumn3d::getExternalNodes(void) {
   return connectedExternalNodes;
}

Node ** MixedBeamColumn3d::getNodePtrs(void) {
   return theNodes;
}

int MixedBeamColumn3d::getNumDOF(void) {
   return NEGD;
}

void MixedBeamColumn3d::setDomain(Domain *theDomain) {

  // check Domain is not null - invoked when object removed from a domain
  if (theDomain == 0) {
    theNodes[0] = 0;
    theNodes[1] = 0;

    opserr << "MixedBeamColumn3d::setDomain:  theDomain = 0 ";
    exit(0);
  }

  // get pointers to the nodes
  int Nd1 = connectedExternalNodes(0);
  int Nd2 = connectedExternalNodes(1);

  theNodes[0] = theDomain->getNode(Nd1);
  theNodes[1] = theDomain->getNode(Nd2);

  if (theNodes[0] == 0) {
    opserr << "MixedBeamColumn3d::setDomain: Nd1: ";
    opserr << Nd1 << "does not exist in model\n";
    exit(0);
  }

  if (theNodes[1] == 0) {
    opserr << "MixedBeamColumn3d::setDomain: Nd2: ";
    opserr << Nd2 << "does not exist in model\n";
    exit(0);
  }

  // call the DomainComponent class method
  this->DomainComponent::setDomain(theDomain);

  // ensure connected nodes have correct number of dof's
  int dofNode1 = theNodes[0]->getNumberDOF();
  int dofNode2 = theNodes[1]->getNumberDOF();

  if ((dofNode1 != NND) || (dofNode2 != NND)) {
    opserr << "MixedBeamColumn3d::setDomain(): Nd2 or Nd1 incorrect dof ";
    exit(0);
  }

  // initialize the transformation
  if (crdTransf->initialize(theNodes[0], theNodes[1])) {
    opserr << "MixedBeamColumn3d::setDomain(): Error initializing coordinate transformation";
    exit(0);
  }

  // Check element length
  if (crdTransf->getInitialLength() == 0.0) {
    opserr << "MixedBeamColumn3d::setDomain(): Zero element length:" << this->getTag();
    exit(0);
  }

}

int MixedBeamColumn3d::commitState() {
  int err = 0; // error flag
  int i = 0; // integer for loops

  // call element commitState to do any base class stuff
  if ((err = this->Element::commitState()) != 0) {
    opserr << "MixedBeamColumn3d::commitState () - failed in base class";
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
  committedInternalForceOpenSees = internalForceOpenSees;
  commitedNaturalForce = naturalForce;
  commitedLastNaturalDisp = lastNaturalDisp;
  commitedHinv = Hinv;
  commitedGMH = GMH;
  kvcommit = kv;
  for( i = 0; i < numSections; i++){
    commitedSectionForceFibers[i] = sectionForceFibers[i];
    commitedSectionDefFibers[i] = sectionDefFibers[i];
    commitedSectionFlexibility[i] = sectionFlexibility[i];
  }

  // Reset iteration counter
  itr = 0;

  return err;
}


int MixedBeamColumn3d::revertToLastCommit() {
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
  internalForceOpenSees = committedInternalForceOpenSees;
  naturalForce = commitedNaturalForce;
  lastNaturalDisp = commitedLastNaturalDisp;
  Hinv = commitedHinv;
  GMH = commitedGMH;
  kv   = kvcommit;
  for( i = 0; i < numSections; i++){
    sectionForceFibers[i] = commitedSectionForceFibers[i];
    sectionDefFibers[i] = commitedSectionDefFibers[i];
    sectionFlexibility[i] = commitedSectionFlexibility[i];
  }

  // Reset iteration counter
  itr = 0;

  return err;
}


int MixedBeamColumn3d::revertToStart()
{
  int err;
  int i,j,k; // for loops
  i = 0;

  // revert the sections state to start
  do {
     err = sections[i++]->revertToStart();
  } while (err == 0 && i < numSections);

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
  Vector myZeros(NDM_NATURAL);
  myZeros.Zero();

  // Set initial shape functions
  for ( i = 0; i < numSections; i++ ){
    nldhat[i] = this->getNld_hat(i, myZeros, initialLength, geomLinear);
    nd1[i] = this->getNd1(i, myZeros, initialLength, geomLinear);
    nd2[i] = this->getNd2(i, 0, initialLength);

    for( j = 0; j < NDM_SECTION; j++ ){
      for( k = 0; k < NDM_NATURAL; k++ ){
        nd1T[i](k,j) = nd1[i](j,k);
        nd2T[i](k,j) = nd2[i](j,k);
      }
    }
  }

  // Set initial and committed section flexibility and GJ
  Matrix ks(NDM_SECTION,NDM_SECTION);
  double GJ;
  for ( i = 0; i < numSections; i++ ){
    getSectionTangent(i,2,ks,GJ);
    invertMatrix(NDM_SECTION,ks,sectionFlexibility[i]);
    commitedSectionFlexibility[i] = sectionFlexibility[i];
  }

  // Set initial and committed section forces and deformations
  for ( i = 0; i < numSections; i++ ){
    sectionForceFibers[i].Zero();
    commitedSectionForceFibers[i].Zero();
    sectionDefFibers[i].Zero();
    commitedSectionDefFibers[i].Zero();
  }

  // Compute the following matrices: G, G2, H, H12, H22, Md, Kg
  Matrix G(NDM_NATURAL,NDM_NATURAL);
  Matrix G2(NDM_NATURAL,NDM_NATURAL);
  Matrix H(NDM_NATURAL,NDM_NATURAL);
  Matrix H12(NDM_NATURAL,NDM_NATURAL);
  Matrix H22(NDM_NATURAL,NDM_NATURAL);
  Matrix Md(NDM_NATURAL,NDM_NATURAL);
  Matrix Kg(NDM_NATURAL,NDM_NATURAL);

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
  commitedHinv = Hinv;

  // Compute the GMH matrix ( G + Md - H12 ) and its transpose
  GMH = G + Md - H12;
  //GMH = G; // Omit P-small delta
  commitedGMH = GMH;

  // Compute the transposes of the following matrices: G2, GMH
  Matrix G2T(NDM_NATURAL,NDM_NATURAL);
  Matrix GMHT(NDM_NATURAL,NDM_NATURAL);
  for( i = 0; i < NDM_NATURAL; i++ ){
    for( j = 0; j < NDM_NATURAL; j++ ){
      G2T(i,j) = G2(j,i);
      GMHT(i,j) = GMH(j,i);
    }
  }

  // Compute the stiffness matrix without the torsion term
  Matrix K_temp_noT(NDM_NATURAL,NDM_NATURAL);
  K_temp_noT = ( Kg + G2 + G2T - H22 ) + GMHT * Hinv * GMH;
  //K_temp_noT = ( Kg ) + GMHT * Hinv * GMH; // Omit P-small delta

  // Add in the torsional stiffness term
  kv.Zero();
  for( i = 0; i < NDM_NATURAL; i++ ) {
    for( j = 0; j < NDM_NATURAL; j++ ) {
      kv(i,j) = K_temp_noT(i,j);
    }
  }

  kv(5,5) =  GJ/initialLength; // Torsional Stiffness GJ/L
  kvcommit = kv;

  Matrix kvOpenSees = transformNaturalCoordsT*kv*transformNaturalCoords;
  Ki = new Matrix(crdTransf->getInitialGlobalStiffMatrix(kvOpenSees));

  // Vector V is zero at initial state
  V.Zero();
  committedV.Zero();

  // Internal force is zero at initial state
  internalForceOpenSees.Zero();
  committedInternalForceOpenSees.Zero();
  naturalForce.Zero();
  commitedNaturalForce.Zero();

  // Last natural displacement is zero at initial state
  lastNaturalDisp.Zero();
  commitedLastNaturalDisp.Zero();

  // Reset iteration counter
  itr = 0;

  // Set initialFlag to 1 so update doesn't call again
  initialFlag = 1;

  return err;
}

const Matrix & MixedBeamColumn3d::getInitialStiff(void) {
  // If things haven't be initialized, then do so
  if (initialFlag == 0) {
    this->revertToStart();
  }
  return *Ki;
}

const Matrix & MixedBeamColumn3d::getTangentStiff(void) {
  // If things haven't be initialized, then do so
  if (initialFlag == 0) {
    this->revertToStart();
  }
  crdTransf->update();  // Will remove once we clean up the corotational 3d transformation -- MHS
  Matrix ktOpenSees = transformNaturalCoordsT*kv*transformNaturalCoords;
  return crdTransf->getGlobalStiffMatrix(ktOpenSees,internalForceOpenSees);
}

const Vector & MixedBeamColumn3d::getResistingForce(void) {
  crdTransf->update();  // Will remove once we clean up the corotational 3d transformation -- MHS
  Vector p0Vec(p0, NDM_NATURAL);
  return crdTransf->getGlobalResistingForce(internalForceOpenSees, p0Vec);
}

int MixedBeamColumn3d::update() {

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
  Vector naturalDispWithTorsion = crdTransf->getBasicTrialDisp();
  naturalDispWithTorsion = transformNaturalCoords*naturalDispWithTorsion;
    // convert to the arrangement of natural deformations that the element likes

  Vector naturalDisp(NDM_NATURAL);
  for ( i = 0; i < NDM_NATURAL; i++ ) {
    naturalDisp(i) = naturalDispWithTorsion(i); //all but the torsional component
  }
  double twist = naturalDispWithTorsion(5);

  Vector naturalIncrDeltaDisp(NDM_NATURAL);
  naturalIncrDeltaDisp = naturalDisp - lastNaturalDisp;
  lastNaturalDisp = naturalDisp;

  // Get the numerical integration weights
  double wt[MAX_NUM_SECTIONS]; // weights of sections or gauss points of integration points
  beamIntegr->getSectionWeights(numSections, initialLength, wt);

  // Define Variables
  double GJ;
  double torsionalForce;
  for ( i = 0; i < numSections; i++ ) {
    sectionForceShapeFcn[i] = Vector(NDM_SECTION);
  }

  // Compute shape functions and their transposes
  for ( i = 0; i < numSections; i++ ){
    // Shape Functions
    nldhat[i] = this->getNld_hat(i, naturalDisp, currentLength, geomLinear);
    sectionDefShapeFcn[i] = this->getd_hat(i, naturalDisp, currentLength, geomLinear);
    nd1[i] = this->getNd1(i, naturalDisp, currentLength, geomLinear);
    if (geomLinear) {
      nd2[i].Zero();
    } else {
      nd2[i] = this->getNd2(i, internalForceOpenSees(0), currentLength);
    }

    // Transpose of shape functions
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
    double torsionalStrain = twist/currentLength;
    setSectionDeformation(i,sectionDefFibers[i],torsionalStrain);

    // Get section force vector
    double tempTorsionalForce;
    getSectionStress(i,sectionForceFibers[i],tempTorsionalForce);
    if (i == 0) {
      torsionalForce = tempTorsionalForce;
    }

    // Get section tangent matrix
    Matrix ks(NDM_SECTION,NDM_SECTION);
    getSectionTangent(i,1,ks,GJ);

    // Compute section flexibility matrix
    invertMatrix(NDM_SECTION,ks,sectionFlexibility[i]);
  }

  // Compute the following matrices: V, V2, G, G2, H, H12, H22, Md, Kg
  Vector V2(NDM_NATURAL);
  Matrix G(NDM_NATURAL,NDM_NATURAL);
  Matrix G2(NDM_NATURAL,NDM_NATURAL);
  Matrix H(NDM_NATURAL,NDM_NATURAL);
  Matrix H12(NDM_NATURAL,NDM_NATURAL);
  Matrix H22(NDM_NATURAL,NDM_NATURAL);
  Matrix Md(NDM_NATURAL,NDM_NATURAL);
  Matrix Kg(NDM_NATURAL,NDM_NATURAL);

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
  Matrix GT(NDM_NATURAL,NDM_NATURAL);
  Matrix G2T(NDM_NATURAL,NDM_NATURAL);
  Matrix GMHT(NDM_NATURAL,NDM_NATURAL);
  for( i = 0; i < NDM_NATURAL; i++ ){
    for( j = 0; j < NDM_NATURAL; j++ ){
      GT(i,j) = G(j,i);
      G2T(i,j) = G2(j,i);
      GMHT(i,j) = GMH(j,i);
    }
  }

  // Compute new internal force
  Vector internalForce(NDM_NATURAL);
  internalForce.Zero();

  if (geomLinear) {
    internalForce = GT * naturalForce + V2 + GMHT * Hinv * V;
  } else {
    internalForce = GT * naturalForce + V2 + GMHT * Hinv * V;
  }

  // Compute internal force for OpenSees ( i.e., add torsion and rearrange )
  for ( i = 0; i < NDM_NATURAL; i++ ) {
    internalForceOpenSees(i) = internalForce(i);
  }
  internalForceOpenSees(5) = torsionalForce; // Add in torsional force
  internalForceOpenSees = transformNaturalCoordsT*internalForceOpenSees;

  // Compute the stiffness matrix without the torsion term
  Matrix K_temp(NDM_NATURAL,NDM_NATURAL);
  if (geomLinear) {
    K_temp = ( Kg + G2 + G2T - H22 ) + GMHT * Hinv * GMH;
  } else {
    K_temp = ( Kg + G2 + G2T - H22 ) + GMHT * Hinv * GMH;
  }

  // Add in the torsional stiffness term
  kv.Zero();
  for( i = 0; i < NDM_NATURAL; i++ ) {
    for( j = 0; j< NDM_NATURAL; j++ ) {
      kv(i,j) = K_temp(i,j);
    }
  }

  kv(5,5) =  GJ/currentLength; // Torsional Stiffness GJ/L

  return 0;
}

const Matrix & MixedBeamColumn3d::getMass(void) {
  theMatrix.Zero();

  if (rho != 0.0) {
    theMatrix(0,0) = theMatrix(1,1) = theMatrix(2,2) =
    theMatrix(6,6) = theMatrix(7,7) = theMatrix(8,8) = 0.5*initialLength*rho;
  }

  return theMatrix;
}

const Matrix & MixedBeamColumn3d::getDamp(void) {
  theMatrix.Zero();

  // Add the damping forces
  if ( doRayleigh == 1 && (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0) ) {
    theMatrix = this->Element::getDamp();
  }

  return theMatrix;
}

void MixedBeamColumn3d::zeroLoad(void) {
  if (sp != 0)
    sp->Zero();

  p0[0] = 0.0;
  p0[1] = 0.0;
  p0[2] = 0.0;
  p0[3] = 0.0;
  p0[4] = 0.0;
}

int MixedBeamColumn3d::addLoad(ElementalLoad *theLoad, double loadFactor) {

  int type;
  const Vector &data = theLoad->getData(type, loadFactor);

  if (sp == 0) {
    sp = new Matrix(NDM_SECTION,numSections);
    if (sp == 0) {
      opserr << "MixedBeamColumn3d::addLoad -- out of memory\n";
      exit(-1);
    }
  }

  double L = crdTransf->getInitialLength();

  double xi[MAX_NUM_SECTIONS];
  beamIntegr->getSectionLocations(numSections, L, xi);

  if (type == LOAD_TAG_Beam3dUniformLoad) {
    double wy = data(0)*loadFactor;  // Transverse
    double wz = data(1)*loadFactor;  // Transverse
    double wx = data(2)*loadFactor;  // Axial

    Matrix &s_p = *sp;

    // Accumulate applied section forces due to element loads
    for (int i = 0; i < numSections; i++) {
      double x = xi[i]*L;
      // Axial
      s_p(0,i) += wx*(L-x);
      // Moment
      s_p(1,i) += wy*0.5*x*(x-L);
      // Moment
      s_p(2,i) += wz*0.5*x*(L-x);
    }

    // Accumulate reactions in basic system
    p0[0] -= wx*L;
    double V;
    V = 0.5*wy*L;
    p0[1] -= V;
    p0[2] -= V;
    V = 0.5*wz*L;
    p0[3] -= V;
    p0[4] -= V;


  } else if (type == LOAD_TAG_Beam3dPointLoad) {
    double Py = data(0)*loadFactor;
    double Pz = data(1)*loadFactor;
    double N  = data(2)*loadFactor;
    double aOverL = data(3);

    if (aOverL < 0.0 || aOverL > 1.0)
      return 0;

    double a = aOverL*L;

    double Vy2 = Py*aOverL;
    double Vy1 = Py-Vy2;

    double Vz2 = Pz*aOverL;
    double Vz1 = Pz-Vz2;

    Matrix &s_p = *sp;

    // Accumulate applied section forces due to element loads
    for (int i = 0; i < numSections; i++) {
      double x = xi[i]*L;
      if (x <= a) {
        s_p(0,i) += N;
        s_p(1,i) -= x*Vy1;
        s_p(2,i) += x*Vz1;
      }
      else {
        s_p(1,i) -= (L-x)*Vy2;
        s_p(2,i) += (L-x)*Vz2;
      }
    }

    // Accumulate reactions in basic system
    p0[0] -= N;
    p0[1] -= Vy1;
    p0[2] -= Vy2;
    p0[3] -= Vz1;
    p0[4] -= Vz2;


  } else {
    opserr << "MixedBeamColumn3d::addLoad() -- load type unknown for element with tag: " <<
        this->getTag() << endln;

    return -1;
  }

  return 0;
}

const Vector & MixedBeamColumn3d::getResistingForceIncInertia() {

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
    theVector(2) += m*accel1(2);
    theVector(6) += m*accel2(0);
    theVector(7) += m*accel2(1);
    theVector(8) += m*accel2(2);
  }

  // Add the damping forces
  if ( doRayleigh == 1 && (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0) ) {
    theVector += this->getRayleighDampingForces();
  }

  return theVector;
}


void MixedBeamColumn3d::Print(OPS_Stream &s, int flag) {

  if (flag == 1) {
    s << "\nElement: " << this->getTag() << " Type: MixedBeamColumn3d ";
    s << "\tConnected Nodes: " << connectedExternalNodes ;
    s << "\tNumber of Sections: " << numSections;
    s << "\tMass density: " << rho;
    for (int i = 0; i < numSections; i++)
      s << "\nSection "<<i<<" :" << *sections[i];

  } else if (flag == 33) {
    s << "\nElement: " << this->getTag() << " Type: MixedBeamColumn3d ";
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
    s << "\"type\": \"mixedBeamColumn2d\", ";
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
    s << "\nElement: " << this->getTag() << " Type: MixedBeamColumn3d ";
    s << "\tConnected Nodes: " << connectedExternalNodes ;
    s << "\tNumber of Sections: " << numSections;
    s << "\tMass density: " << rho << endln;
  }

}


OPS_Stream &operator<<(OPS_Stream &s, MixedBeamColumn3d &E) {
  E.Print(s);
  return s;
}


Response* MixedBeamColumn3d::setResponse(const char **argv, int argc,
                                         OPS_Stream &output) {

  Response *theResponse = 0;

  output.tag("ElementOutput");
  output.attr("eleType","MixedBeamColumn3d");
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
  } else if (strcmp(argv[0],"localForce") == 0 ||
             strcmp(argv[0],"localForces") == 0) {

    output.tag("ResponseType","N_ 1");
    output.tag("ResponseType","Vy_1");
    output.tag("ResponseType","Vz_1");
    output.tag("ResponseType","T_1");
    output.tag("ResponseType","My_1");
    output.tag("ResponseType","Tz_1");
    output.tag("ResponseType","N_2");
    output.tag("ResponseType","Py_2");
    output.tag("ResponseType","Pz_2");
    output.tag("ResponseType","T_2");
    output.tag("ResponseType","My_2");
    output.tag("ResponseType","Mz_2");

    theResponse = new ElementResponse(this, 2, theVector);

  // basic or natural forces
  } else if (strcmp(argv[0],"basicForce") == 0 ||
             strcmp(argv[0],"basicForces") == 0) {

    output.tag("ResponseType","N");
    output.tag("ResponseType","Mz_1");
    output.tag("ResponseType","Mz_2");
    output.tag("ResponseType","My_1");
    output.tag("ResponseType","My_2");
    output.tag("ResponseType","T");

    theResponse = new ElementResponse(this, 3, Vector(6));
  } else if (strcmp(argv[0],"sectionDeformation_Force") == 0) {

    int i;
    char *q  = new char[15];
    for ( i = 0; i < numSections; i++ ){
      sprintf(q,"axialStrain_%i",i+1);
      output.tag("ResponseType",q);
      sprintf(q,"curvatureZ_%i",i+1);
      output.tag("ResponseType",q);
      sprintf(q,"curvatureY_%i",i+1);
      output.tag("ResponseType",q);
    }
    delete [] q;

    theResponse =  new ElementResponse(this, 4, Vector(3*numSections));

  } else if (strcmp(argv[0],"plasticSectionDeformation_Force") == 0) {

    int i;
    char *q  = new char[25];
    for ( i = 0; i < numSections; i++ ){
      sprintf(q,"plasticAxialStrain_%i",i+1);
      output.tag("ResponseType",q);
      sprintf(q,"plasticCurvatureZ_%i",i+1);
      output.tag("ResponseType",q);
      sprintf(q,"plasticCurvatureY_%i",i+1);
      output.tag("ResponseType",q);
    }
    delete [] q;

    theResponse =  new ElementResponse(this, 5, Vector(3*numSections));

  } else if (strcmp(argv[0],"integrationPoints") == 0) {
    theResponse =  new ElementResponse(this, 100, Vector(numSections));

  } else if (strcmp(argv[0],"integrationWeights") == 0) {
    theResponse =  new ElementResponse(this, 101, Vector(numSections));

  } else if (strcmp(argv[0],"connectedNodes") == 0) {
    theResponse =  new ElementResponse(this, 102, Vector(2));

  } else if (strcmp(argv[0],"numSections") == 0 ||
             strcmp(argv[0],"numberOfSections") == 0 ) {
    theResponse =  new ElementResponse(this, 103, Vector(1));

  }

    else if (strcmp(argv[0],"xaxis") == 0 || strcmp(argv[0],"xlocal") == 0)
      theResponse = new ElementResponse(this, 201, Vector(3));

    else if (strcmp(argv[0],"yaxis") == 0 || strcmp(argv[0],"ylocal") == 0)
      theResponse = new ElementResponse(this, 202, Vector(3));

    else if (strcmp(argv[0],"zaxis") == 0 || strcmp(argv[0],"zlocal") == 0)
      theResponse = new ElementResponse(this, 203, Vector(3));  

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

        theResponse =  sections[sectionNum-1]->setResponse(&argv[2], argc-2, output);

        output.endTag();
      }
    }
  }

  output.endTag();
  return theResponse;
}


int MixedBeamColumn3d::getResponse(int responseID, Information &eleInfo) {
  if (responseID == 1) { // global forces
    return eleInfo.setVector(this->getResistingForce());

  } else if (responseID == 2) { // local forces
    // Axial
    double N = internalForceOpenSees(0);
    theVector(6) =  N;
    theVector(0) = -N+p0[0];

    // Torsion
    double T = internalForceOpenSees(5);
    theVector(9) =  T;
    theVector(3) = -T;

    // Moments about z and shears along y
    double M1 = internalForceOpenSees(1);
    double M2 = internalForceOpenSees(2);
    theVector(5)  = M1;
    theVector(11) = M2;
    double L = crdTransf->getInitialLength();
    double V = (M1+M2)/L;
    theVector(1) =  V+p0[1];
    theVector(7) = -V+p0[2];

    // Moments about y and shears along z
    M1 = internalForceOpenSees(3);
    M2 = internalForceOpenSees(4);
    theVector(4)  = M1;
    theVector(10) = M2;
    V = -(M1+M2)/L;
    theVector(2) = -V+p0[3];
    theVector(8) =  V+p0[4];

    return eleInfo.setVector(theVector);

  } else if (responseID == 3) { // basic forces
    return eleInfo.setVector(internalForceOpenSees);

  } else if (responseID == 4) { // section deformation (from forces)

    int i;
    Vector tempVector(3*numSections);
    tempVector.Zero();
    for ( i = 0; i < numSections; i++ ){
      tempVector(3*i)   = sectionDefFibers[i](0);
      tempVector(3*i+1) = sectionDefFibers[i](1);
      tempVector(3*i+2) = sectionDefFibers[i](2);
    }

    return eleInfo.setVector(tempVector);

  } else if (responseID == 5) { // plastic section deformation (from forces)

    int i;
    Vector tempVector(3*numSections);
    Vector sectionForce(NDM_SECTION);
    Vector plasticSectionDef(NDM_SECTION);
    Matrix ks(3,3);
    Matrix fs(3,3);
    tempVector.Zero();
    double scratch = 0.0;
    for ( i = 0; i < numSections; i++ ){

      getSectionStress(i,sectionForce,scratch);
      getSectionTangent(i,2,ks,scratch);
      invertMatrix(NDM_SECTION,ks,fs);

      plasticSectionDef = sectionDefFibers[i] - fs*sectionForce;

      tempVector(3*i)   = plasticSectionDef(0);
      tempVector(3*i+1) = plasticSectionDef(1);
      tempVector(3*i+2) = plasticSectionDef(2);
    }

    return eleInfo.setVector(tempVector);

  } else if (responseID == 100) { // integration points

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

  } else if (responseID == 102) { // connected nodes
    Vector tempVector(2);
    tempVector(0) = connectedExternalNodes(0);
    tempVector(1) = connectedExternalNodes(1);
    return eleInfo.setVector(tempVector);

  } else if (responseID == 103) { // number of sections
    Vector tempVector(1);
    tempVector(0) = numSections;
    return eleInfo.setVector(tempVector);

  }

  else if (responseID >= 201 && responseID <= 203) {
    static Vector xlocal(3);
    static Vector ylocal(3);
    static Vector zlocal(3);

    crdTransf->getLocalAxes(xlocal,ylocal,zlocal);
    
    if (responseID == 201)
      return eleInfo.setVector(xlocal);
    if (responseID == 202)
      return eleInfo.setVector(ylocal);
    if (responseID == 203)
      return eleInfo.setVector(zlocal);    
  }

  else {
    return -1;
  }
}

Vector MixedBeamColumn3d::getd_hat(int sec, const Vector &v, double L, bool geomLinear) {
  double xi[MAX_NUM_SECTIONS];
  beamIntegr->getSectionLocations(numSections, L, xi);

  double x, C, E, F;
  Vector D_hat(NDM_SECTION);
  D_hat.Zero();

  x = L*xi[sec];
  C =  1/L;
  E = -4/L + 6*x/(L*L);
  F = -2/L + 6*x/(L*L);

  if (geomLinear) {

    D_hat(0) = C*v(0);
    D_hat(1) = E*v(1) + F*v(3);
    D_hat(2) = E*v(2) + F*v(4);

  } else {

    double A,B;
    A = 1 - 4*(x/L) + 3*pow(x/L,2);
    B =   - 2*(x/L) + 3*pow(x/L,2);

    D_hat(0) =  C * v(0) + 0.5 * ( C*C*v(0) ) * v(0) +
                0.5 * ( A*A*v(1) + A*B*v(3) ) * v(1) +
                0.5 * ( A*A*v(2) + A*B*v(4) ) * v(2) +
                0.5 * ( A*B*v(1) + B*B*v(3) ) * v(3) +
                0.5 * ( A*B*v(2) + B*B*v(4) ) * v(4);
    D_hat(1) =  E*v(1) + F*v(3);
    D_hat(2) =  E*v(2) + F*v(4);

  }

  return D_hat;
}

Matrix MixedBeamColumn3d::getKg(int sec, double P, double L) {
  double xi[MAX_NUM_SECTIONS];
  beamIntegr->getSectionLocations(numSections, L, xi);

  double temp_x, temp_A, temp_B;

  temp_x = L * xi[sec];

  Matrix kg(NDM_NATURAL,NDM_NATURAL);
  kg.Zero();

  temp_A = 1 - 4 * temp_x / L + 3 * ( temp_x * temp_x ) / ( L * L );
  temp_B = - 2 * temp_x / L + 3 * ( temp_x * temp_x ) / ( L * L );

  kg(0,0) = P / ( L * L );
  kg(1,1) = P * temp_A * temp_A;
  kg(1,3) = P * temp_A * temp_B;
  kg(2,2) = P * temp_A * temp_A;
  kg(2,4) = P * temp_A * temp_B;
  kg(3,1) = P * temp_A * temp_B;
  kg(3,3) = P * temp_B * temp_B;
  kg(4,2) = P * temp_A * temp_B;
  kg(4,4) = P * temp_B * temp_B;

  return kg;
}

Matrix MixedBeamColumn3d::getMd(int sec, Vector dShapeFcn, Vector dFibers, double L) {
  double xi[MAX_NUM_SECTIONS];
  beamIntegr->getSectionLocations(numSections, L, xi);

  double x, A, B;

  Matrix md(NDM_NATURAL,NDM_NATURAL);
  md.Zero();

  x = L*xi[sec];
  A =  ( x/L - 2*pow(x/L,2) + pow(x/L,3) )*L;
  B =          (-pow(x/L,2) + pow(x/L,3) )*L;

  md(0,1) = A * ( dShapeFcn(1) - dFibers(1) );
  md(0,2) = A * ( dShapeFcn(2) - dFibers(2) );
  md(0,3) = B * ( dShapeFcn(1) - dFibers(1) );
  md(0,4) = B * ( dShapeFcn(2) - dFibers(2) );

  return md;
}

Matrix MixedBeamColumn3d::getNld_hat(int sec, const Vector &v, double L, bool geomLinear) {
  double xi[MAX_NUM_SECTIONS];
  beamIntegr->getSectionLocations(numSections, L, xi);

  double x, C, E, F;
  Matrix Nld_hat(NDM_SECTION,NDM_NATURAL);
  Nld_hat.Zero();

  x = L*xi[sec];

  C =  1/L;
  E = -4/L + 6*x/(L*L);
  F = -2/L + 6*x/(L*L);

  if (geomLinear) {

    Nld_hat(0,0) = C;
    Nld_hat(1,1) = E;
    Nld_hat(1,3) = F;
    Nld_hat(2,2) = E;
    Nld_hat(2,4) = F;

  } else {

    double A,B;
    A = 1 - 4 * ( x / L ) + 3 * pow ( ( x / L ) , 2 );
    B =   - 2 * ( x / L ) + 3 * pow ( ( x / L ) , 2 );

    Nld_hat(0,0) = C + C*C*v(0);
    Nld_hat(0,1) = A*A*v(1) + A*B*v(3);
    Nld_hat(0,2) = A*A*v(2) + A*B*v(4);
    Nld_hat(0,3) = A*B*v(1) + B*B*v(3);
    Nld_hat(0,4) = A*B*v(2) + B*B*v(4);
    Nld_hat(1,1) = E;
    Nld_hat(1,3) = F;
    Nld_hat(2,2) = E;
    Nld_hat(2,4) = F;
  }

  return Nld_hat;
}

Matrix MixedBeamColumn3d::getNd2(int sec, double P, double L) {
  double xi[MAX_NUM_SECTIONS];
  beamIntegr->getSectionLocations(numSections, L, xi);

  double temp_x, temp_A, temp_B;

  temp_x = L * xi[sec];

  Matrix Nd2(NDM_SECTION,NDM_NATURAL);
  Nd2.Zero();

  temp_A = L * ( temp_x / L - 2 * pow( temp_x / L, 2 ) + pow( temp_x / L, 3 ) );
  temp_B = L * ( -pow( temp_x / L, 2 ) + pow( temp_x / L, 3 ) );

  Nd2(1,1) = P * temp_A;
  Nd2(1,3) = P * temp_B;
  Nd2(2,2) = P * temp_A;
  Nd2(2,4) = P * temp_B;

  return Nd2;
}

Matrix MixedBeamColumn3d::getNd1(int sec, const Vector &v, double L, bool geomLinear) {
  double xi[MAX_NUM_SECTIONS];
  beamIntegr->getSectionLocations(numSections, L, xi);

  double x = L*xi[sec];

  Matrix Nd1(NDM_SECTION,NDM_NATURAL);
  Nd1.Zero();

  if (geomLinear) {

    Nd1(0,0)   = 1.0;
    Nd1(1,1)   = -x/L + 1.0;
    Nd1(1,3)   =  x/L;
    Nd1(2,2)   = -x/L + 1.0;
    Nd1(2,4)   =  x/L;

  } else {

    double A,B;

    A = L * ( x/L - 2*pow(x/L,2) + pow(x/L,3) ) * v[1]
             + L * ( -pow(x/L,2) + pow(x/L,3) ) * v[3];

    B = L * ( x/L - 2*pow(x/L,2) + pow(x/L,3) ) * v[2]
             + L * ( -pow(x/L,2) + pow(x/L,3) ) * v[4];

    Nd1(0,0)   = 1.0;
    Nd1(1,0)   = A;
    Nd1(1,1)   = -x/L + 1.0;
    Nd1(1,3)   =  x/L;
    Nd1(2,0)   = B;
    Nd1(2,2)   = -x/L + 1.0;
    Nd1(2,4)   =  x/L;
  }

  return Nd1;
}

void MixedBeamColumn3d::getSectionTangent(int sec,int type,Matrix &kSection,
                                          double &GJ) {
  int order = sections[sec]->getOrder();
  const ID &code = sections[sec]->getType();

  // Initialize formulation friendly variables
  kSection.Zero();
  GJ = 0.0;

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
            case SECTION_RESPONSE_MY:
              kSection(0,2) = sectionTangent(i,j);
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
            case SECTION_RESPONSE_MY:
              kSection(1,2) = sectionTangent(i,j);
              break;
            default:
              break;
          }
          break;
        case SECTION_RESPONSE_MY:
          switch(code(j)) {
            case SECTION_RESPONSE_P:
              kSection(2,0) = sectionTangent(i,j);
              break;
            case SECTION_RESPONSE_MZ:
              kSection(2,1) = sectionTangent(i,j);
              break;
            case SECTION_RESPONSE_MY:
              kSection(2,2) = sectionTangent(i,j);
              break;
            default:
              break;
          }
          break;
        case SECTION_RESPONSE_T:
          GJ = sectionTangent(i,i);
          break;
        default:
          break;
      }
    }
  }
}

void MixedBeamColumn3d::getSectionStress(int sec,Vector &fSection,
                                         double &torsion) {
  int order = sections[sec]->getOrder();
  const ID &code = sections[sec]->getType();

  // Get the stress resultant from section
  Vector stressResultant = sections[sec]->getStressResultant();

  // Initialize formulation friendly variables
  fSection.Zero();
  torsion = 0.0;

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
      case SECTION_RESPONSE_MY:
        fSection(2) = stressResultant(j);
        break;
      case SECTION_RESPONSE_T:
        torsion = stressResultant(j);
        break;
      default:
        break;
    }
  }
}

void MixedBeamColumn3d::setSectionDeformation(int sec,Vector &defSection,
                                              double &twist) {
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
      case SECTION_RESPONSE_MY:
        sectionDeformation(j) = defSection(2);
        break;
      case SECTION_RESPONSE_T:
        sectionDeformation(j) = twist;
        break;
      default:
        break;
    }
  }

  // Set the section deformations
  int res = sections[sec]->setTrialSectionDeformation(sectionDeformation);
}


int MixedBeamColumn3d::sendSelf(int commitTag, Channel &theChannel) {
  // @todo write MixedBeamColumn3d::sendSelf
  opserr << "Error: MixedBeamColumn3d::sendSelf -- not yet implemented for MixedBeamColumn3d element";
  return -1;
}

int MixedBeamColumn3d::recvSelf(int commitTag, Channel &theChannel,
                                FEM_ObjectBroker &theBroker) {
  // @todo write MixedBeamColumn3d::recvSelf
  opserr << "Error: MixedBeamColumn3d::sendSelf -- not yet implemented for MixedBeamColumn3d element";
  return -1;
}

int MixedBeamColumn3d::displaySelf(Renderer& theViewer, int displayMode, float fact, const char** modes, int numMode)
{
    static Vector v1(3);
    static Vector v2(3);

    theNodes[0]->getDisplayCrds(v1, fact, displayMode);
    theNodes[1]->getDisplayCrds(v2, fact, displayMode);

    return theViewer.drawLine(v1, v2, 1.0, 1.0, this->getTag());
}
