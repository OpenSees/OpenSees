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
                                                                        
// $Revision: 1.35 $
// $Date: 2009-10-13 21:14:21 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/fourNodeQuad/FourNodeQuad.cpp,v $

// Written: MHS
// Created: Feb 2000
// Revised: Dec 2000 for efficiency
//
// Description: This file contains the class definition for FourNodeQuad.

#include <FourNodeQuad.h>
#include <Node.h>
#include <NDMaterial.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <Renderer.h>
#include <Domain.h>
#include <string.h>
#include <Information.h>
#include <Parameter.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <ElementResponse.h>
#include <ElementalLoad.h>
#include <elementAPI.h>

void* OPS_FourNodeQuad()
{
    int dampingTag = 0;
    Damping* theDamping = 0;

    int ndm = OPS_GetNDM();
    int ndf = OPS_GetNDF();

    if (ndm != 2 || ndf != 2) {
	opserr << "WARNING -- model dimensions and/or nodal DOF not compatible with quad element\n";
	return 0;
    }
    
    if (OPS_GetNumRemainingInputArgs() < 8) {
	opserr << "WARNING insufficient arguments\n";
	opserr << "Want: element FourNodeQuad eleTag? iNode? jNode? kNode? lNode? thk? type? matTag? <pressure? rho? b1? b2?>\n";
	return 0;
    }

    // FourNodeQuadId, iNode, jNode, kNode, lNode
    int idata[5];
    int num = 5;
    if (OPS_GetIntInput(&num,idata) < 0) {
	opserr<<"WARNING: invalid integer inputs\n";
	return 0;
    }

    double thk = 1.0;
    num = 1;
    if (OPS_GetDoubleInput(&num,&thk) < 0) {
	opserr<<"WARNING: invalid double inputs\n";
	return 0;
    }

    const char* type = OPS_GetString();

    int matTag;
    num = 1;
    if (OPS_GetIntInput(&num,&matTag) < 0) {
	opserr<<"WARNING: invalid matTag\n";
	return 0;
    }

    NDMaterial* mat = OPS_getNDMaterial(matTag);
    if (mat == 0) {
	opserr << "WARNING material not found\n";
	opserr << "Material: " << matTag;
	opserr << "\nFourNodeQuad element: " << idata[0] << endln;
	return 0;
    }

    // p, rho, b1, b2
    double data[4] = {0,0,0,0};
    num = OPS_GetNumRemainingInputArgs();
    if (num > 4) {
	num = 4;
    }
    if (num > 0) {
	if (OPS_GetDoubleInput(&num,data) < 0) {
	    opserr<<"WARNING: invalid integer data\n";
	    return 0;
	}	
    }

    //option,written by Tang.S
    int numData = 1;
    while (OPS_GetNumRemainingInputArgs() > 0) {
        std::string theType = OPS_GetString();
        if (theType == "-damp") {

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
    //Add theDamping at the end,Tang.S
    return new FourNodeQuad(idata[0], idata[1], idata[2], idata[3], idata[4],
        *mat, type, thk, data[0], data[1], data[2], data[3], theDamping);
}

void *OPS_FourNodeQuad(const ID &info) {
    if (info.Size() == 0) {
        opserr << "WARNING: info is empty -- FourNodeQuad\n";
        return 0;
    }

    int dampingTag = 0;
    Damping* theDamping = 0;

    int ndm = OPS_GetNDM();
    int ndf = OPS_GetNDF();

    // mesh data
    static std::map<int, Vector> meshdata;

    // save data
    if (info(0) == 1) {
        // check input
        if (info.Size() < 2) {
            opserr << "WARNING: need info -- inmesh, meshtag\n";
            return 0;
        }
        if (OPS_GetNumRemainingInputArgs() < 3) {
            opserr << "WARNING: insuficient arguments -- thk? type? "
                      "matTag? <pressure? rho? "
                      "b1? b2?>\n";
            return 0;
        }
        if (ndm != 2 || ndf != 2) {
            opserr
                << "WARNING -- model dimensions and/or nodal DOF not "
                   "compatible with quad element\n";
            return 0;
        }

        // save data
        Vector &mdata = meshdata[info(1)];
        mdata.resize(7);
        mdata.Zero();

        // get thk
        double thk = 1.0;
        int num = 1;
        if (OPS_GetDoubleInput(&num, &thk) < 0) {
            opserr << "WARNING: failed to get thk -- FourNodeQuad\n";
            return 0;
        }
        mdata(0) = thk;

        // get type
        const char *type = OPS_GetString();
        if (strcmp(type, "PlaneStrain") == 0 ||
            strcmp(type, "PlaneStrain2D") == 0) {
            mdata(1) = 1;
        } else if (strcmp(type, "PlaneStress") == 0 ||
                   strcmp(type, "PlaneStress2D") == 0) {
            mdata(1) = 2;
        }

        // get matTag
        int matTag;
        num = 1;
        if (OPS_GetIntInput(&num, &matTag) < 0) {
            opserr << "WARNING: invalid matTag\n";
            return 0;
        }
        mdata(2) = matTag;

        // optional
        // p, rho, b1, b2
        double data[4] = {0, 0, 0, 0};
        num = OPS_GetNumRemainingInputArgs();
        if (num > 4) {
            num = 4;
        }
        if (num > 0) {
            if (OPS_GetDoubleInput(&num, data) < 0) {
                opserr << "WARNING: invalid integer data\n";
                return 0;
            }
        }
        for (int i = 0; i < 4; ++i) {
            mdata(3 + i) = data[i];
        }

        return &meshdata;
    }

    // load data
    if (info(0) == 2) {
        if (info.Size() < 7) {
            opserr << "WARNING: need info -- inmesh, meshtag, "
                      "eleTag, nd1, nd2, nd3, nd4\n";
            return 0;
        }

        int eleTag = info(2);

        // get data
        Vector &mdata = meshdata[info(1)];
        if (mdata.Size() < 7) {
            return 0;
        }

        double thk = mdata(0);
        const char *type = "PlaneStrain";
        if (mdata(1) == 1) {
            type = "PlaneStrain";
        } else {
            type = "PlaneStress";
        }
        int matTag = (int)mdata(2);
        double pressure = mdata(3);
        double rho = mdata(4);
        double b1 = mdata(5);
        double b2 = mdata(6);

        NDMaterial *mat = OPS_getNDMaterial(matTag);
        if (mat == 0) {
            opserr << "WARNING material not found\n";
            opserr << "Material: " << matTag;
            opserr << "\nFourNodeQuad element: " << eleTag << endln;
            return 0;
        }

        //option,written by Tang.S
        int numData = 1;
        while (OPS_GetNumRemainingInputArgs() > 0) {
            std::string theType = OPS_GetString();
            if (theType == "-damp") {

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
        //Add theDamping at the end,Tang.S
        return new FourNodeQuad(eleTag, info(3), info(4), info(5),
                                info(6), *mat, type, thk, pressure,
                                rho, b1, b2, theDamping);
    }

    return 0;
}

double FourNodeQuad::matrixData[64];
Matrix FourNodeQuad::K(matrixData, 8, 8);
Vector FourNodeQuad::P(8);
double FourNodeQuad::shp[3][4];
double FourNodeQuad::pts[4][2];
double FourNodeQuad::wts[4];

FourNodeQuad::FourNodeQuad(int tag, int nd1, int nd2, int nd3, int nd4,
			   NDMaterial &m, const char *type, double t,
			   double p, double r, double b1, double b2,
         Damping *damping)
:Element (tag, ELE_TAG_FourNodeQuad), 
  theMaterial(0), connectedExternalNodes(4), 
 Q(8), pressureLoad(8), thickness(t), applyLoad(0), pressure(p), rho(r), Ki(0)
{
	pts[0][0] = -0.5773502691896258;
	pts[0][1] = -0.5773502691896258;
	pts[1][0] =  0.5773502691896258;
	pts[1][1] = -0.5773502691896258;
	pts[2][0] =  0.5773502691896258;
	pts[2][1] =  0.5773502691896258;
	pts[3][0] = -0.5773502691896258;
	pts[3][1] =  0.5773502691896258;

	wts[0] = 1.0;
	wts[1] = 1.0;
	wts[2] = 1.0;
	wts[3] = 1.0;

	if (strcmp(type,"PlaneStrain") != 0 && strcmp(type,"PlaneStress") != 0
	    && strcmp(type,"PlaneStrain2D") != 0 && strcmp(type,"PlaneStress2D") != 0) {
	  opserr << "FourNodeQuad::FourNodeQuad -- improper material type: " << type << "for FourNodeQuad\n";
	  exit(-1);
	}

	// Body forces
	b[0] = b1;
	b[1] = b2;

    // Allocate arrays of pointers to NDMaterials
    theMaterial = new NDMaterial *[4];
    
    if (theMaterial == 0) {
      opserr << "FourNodeQuad::FourNodeQuad - failed allocate material model pointer\n";
      exit(-1);
    }

	int i;
    for (i = 0; i < 4; i++) {
      
      // Get copies of the material model for each integration point
      theMaterial[i] = m.getCopy(type);
			
      // Check allocation
      if (theMaterial[i] == 0) {
	opserr << "FourNodeQuad::FourNodeQuad -- failed to get a copy of material model\n";
	exit(-1);
      }
    }

    if (damping)
    {
      for (i = 0; i < 4; i++)
      {
        theDamping[i] =(*damping).getCopy();
      
        if (!theDamping[i]) {
          opserr << "FourNodeQuad::FourNodeQuad -- failed to get copy of damping\n";
          exit(-1);
        }
      }
    }
    else
    {
      for (i = 0; i < 4; i++) theDamping[i] = 0;
    }

    // Set connected external node IDs
    connectedExternalNodes(0) = nd1;
    connectedExternalNodes(1) = nd2;
    connectedExternalNodes(2) = nd3;
    connectedExternalNodes(3) = nd4;
    
    for (i=0; i<4; i++)
      theNodes[i] = 0;
}

FourNodeQuad::FourNodeQuad()
:Element (0,ELE_TAG_FourNodeQuad),
  theMaterial(0), connectedExternalNodes(4), 
 Q(8), pressureLoad(8), thickness(0.0), applyLoad(0), pressure(0.0), Ki(0)
{
  pts[0][0] = -0.577350269189626;
  pts[0][1] = -0.577350269189626;
  pts[1][0] =  0.577350269189626;
  pts[1][1] = -0.577350269189626;
  pts[2][0] =  0.577350269189626;
  pts[2][1] =  0.577350269189626;
  pts[3][0] = -0.577350269189626;
  pts[3][1] =  0.577350269189626;
  
  wts[0] = 1.0;
  wts[1] = 1.0;
  wts[2] = 1.0;
  wts[3] = 1.0;

    for (int i=0; i<4; i++)
      theNodes[i] = 0;

  for (int i = 0 ;  i < 4; i++ ) 
    theDamping[i] = 0;
}

FourNodeQuad::~FourNodeQuad()
{    
  for (int i = 0; i < 4; i++) {
    if (theMaterial[i])
      delete theMaterial[i];
  }

  // Delete the array of pointers to NDMaterial pointer arrays
  if (theMaterial)
    delete [] theMaterial;

  for (int i = 0; i < 4; i++)
  {
    if (theDamping[i])
    {
      delete theDamping[i];
      theDamping[i] = 0;
    }
  }

  if (Ki != 0)
    delete Ki;
}

int
FourNodeQuad::getNumExternalNodes() const
{
    return 4;
}

const ID&
FourNodeQuad::getExternalNodes()
{
    return connectedExternalNodes;
}


Node **
FourNodeQuad::getNodePtrs(void) 
{
  return theNodes;
}

int
FourNodeQuad::getNumDOF()
{
    return 8;
}

void
FourNodeQuad::setDomain(Domain *theDomain)
{
	// Check Domain is not null - invoked when object removed from a domain
    if (theDomain == 0) {
	theNodes[0] = 0;
	theNodes[1] = 0;
	theNodes[2] = 0;
	theNodes[3] = 0;
	return;
    }

    int Nd1 = connectedExternalNodes(0);
    int Nd2 = connectedExternalNodes(1);
    int Nd3 = connectedExternalNodes(2);
    int Nd4 = connectedExternalNodes(3);

    theNodes[0] = theDomain->getNode(Nd1);
    theNodes[1] = theDomain->getNode(Nd2);
    theNodes[2] = theDomain->getNode(Nd3);
    theNodes[3] = theDomain->getNode(Nd4);

    if (theNodes[0] == 0 || theNodes[1] == 0 || theNodes[2] == 0 || theNodes[3] == 0) {
	//opserr << "FATAL ERROR FourNodeQuad (tag: %d), node not found in domain",
	//	this->getTag());
	
	return;
    }

    int dofNd1 = theNodes[0]->getNumberDOF();
    int dofNd2 = theNodes[1]->getNumberDOF();
    int dofNd3 = theNodes[2]->getNumberDOF();
    int dofNd4 = theNodes[3]->getNumberDOF();
    
    if (dofNd1 != 2 || dofNd2 != 2 || dofNd3 != 2 || dofNd4 != 2) {
	//opserr << "FATAL ERROR FourNodeQuad (tag: %d), has differing number of DOFs at its nodes",
	//	this->getTag());
	
	return;
    }
    this->DomainComponent::setDomain(theDomain);

    // Compute consistent nodal loads due to pressure
    this->setPressureLoadAtNodes();
    
    for (int i = 0; i < 4; i++)
    {
      if (theDamping[i] && theDamping[i]->setDomain(theDomain, 3)) {
        opserr << "FourNodeQuad::setDomain -- Error initializing damping\n";
        return;
      }
    }

}

int
FourNodeQuad::setDamping(Domain *theDomain, Damping *damping)
{
  if (theDomain && damping)
  {
    for (int i = 0; i < 4; i++)
    {
      if (theDamping[i]) delete theDamping[i];

      theDamping[i] =(*damping).getCopy();
      
      if (!theDamping[i]) {
        opserr << "FourNodeQuad::setDamping -- failed to get copy of damping\n";
        return -1;
      }
      if (theDamping[i]->setDomain(theDomain, 3)) {
        opserr << "FourNodeQuad::setDamping -- Error initializing damping\n";
        return -2;
      }
    }
  }
  return 0;
}

int
FourNodeQuad::commitState()
{
    int retVal = 0;

    // call element commitState to do any base class stuff
    if ((retVal = this->Element::commitState()) != 0) {
      opserr << "FourNodeQuad::commitState () - failed in base class";
    }    

    // Loop over the integration points and commit the material states
    for (int i = 0; i < 4; i++)
      retVal += theMaterial[i]->commitState();

    for (int i = 0; i < 4; i++ )
      if (theDamping[i]) retVal += theDamping[i]->commitState();

    return retVal;
}

int
FourNodeQuad::revertToLastCommit()
{
    int retVal = 0;

    // Loop over the integration points and revert to last committed state
    for (int i = 0; i < 4; i++)
		retVal += theMaterial[i]->revertToLastCommit();

    for (int i = 0; i < 4; i++ )
      if (theDamping[i]) retVal += theDamping[i]->revertToLastCommit();

    return retVal;
}

int
FourNodeQuad::revertToStart()
{
    int retVal = 0;

    // Loop over the integration points and revert states to start
    for (int i = 0; i < 4; i++)
		retVal += theMaterial[i]->revertToStart();

    for (int i = 0; i < 4; i++ )
      if (theDamping[i]) retVal += theDamping[i]->revertToStart();

    return retVal;
}


int
FourNodeQuad::update()
{
	const Vector &disp1 = theNodes[0]->getTrialDisp();
	const Vector &disp2 = theNodes[1]->getTrialDisp();
	const Vector &disp3 = theNodes[2]->getTrialDisp();
	const Vector &disp4 = theNodes[3]->getTrialDisp();
	
	static double u[2][4];

	u[0][0] = disp1(0);
	u[1][0] = disp1(1);
	u[0][1] = disp2(0);
	u[1][1] = disp2(1);
	u[0][2] = disp3(0);
	u[1][2] = disp3(1);
	u[0][3] = disp4(0);
	u[1][3] = disp4(1);

	static Vector eps(3);

	int ret = 0;

	// Loop over the integration points
	for (int i = 0; i < 4; i++) {

		// Determine Jacobian for this integration point
		this->shapeFunction(pts[i][0], pts[i][1]);

		// Interpolate strains
		//eps = B*u;
		//eps.addMatrixVector(0.0, B, u, 1.0);
		eps.Zero();
		for (int beta = 0; beta < 4; beta++) {
			eps(0) += shp[0][beta]*u[0][beta];
			eps(1) += shp[1][beta]*u[1][beta];
			eps(2) += shp[0][beta]*u[1][beta] + shp[1][beta]*u[0][beta];
		}

		// Set the material strain
		ret += theMaterial[i]->setTrialStrain(eps);
	}

	return ret;
}


const Matrix&
FourNodeQuad::getTangentStiff()
{
  static Matrix D(3,3);

	K.Zero();

	double dvol;
	double DB[3][2];

	// Loop over the integration points
	for (int i = 0; i < 4; i++) {

	  // Determine Jacobian for this integration point
	  dvol = this->shapeFunction(pts[i][0], pts[i][1]);
	  dvol *= (thickness*wts[i]);
	  
	  // Get the material tangent
	  D = theMaterial[i]->getTangent();
    if(theDamping[i]) D *= theDamping[i]->getStiffnessMultiplier();
	  
	  // Perform numerical integration
	  //K = K + (B^ D * B) * intWt(i)*intWt(j) * detJ;
	  //K.addMatrixTripleProduct(1.0, B, D, intWt(i)*intWt(j)*detJ);
	  
	  double D00 = D(0,0); double D01 = D(0,1); double D02 = D(0,2);
	  double D10 = D(1,0); double D11 = D(1,1); double D12 = D(1,2);
	  double D20 = D(2,0); double D21 = D(2,1); double D22 = D(2,2);

	  //	  for (int beta = 0, ib = 0, colIb =0, colIbP1 = 8; 
	  //   beta < 4; 
	  //   beta++, ib += 2, colIb += 16, colIbP1 += 16) {
	    
	  for (int alpha = 0, ia = 0; alpha < 4; alpha++, ia += 2) {
	    for (int beta = 0, ib = 0; beta < 4; beta++, ib += 2) {
	      
	      DB[0][0] = dvol * (D00 * shp[0][beta] + D02 * shp[1][beta]);
	      DB[1][0] = dvol * (D10 * shp[0][beta] + D12 * shp[1][beta]);
	      DB[2][0] = dvol * (D20 * shp[0][beta] + D22 * shp[1][beta]);
	      DB[0][1] = dvol * (D01 * shp[1][beta] + D02 * shp[0][beta]);
	      DB[1][1] = dvol * (D11 * shp[1][beta] + D12 * shp[0][beta]);
	      DB[2][1] = dvol * (D21 * shp[1][beta] + D22 * shp[0][beta]);
	      

	      K(ia,ib) += shp[0][alpha]*DB[0][0] + shp[1][alpha]*DB[2][0];
	      K(ia,ib+1) += shp[0][alpha]*DB[0][1] + shp[1][alpha]*DB[2][1];
	      K(ia+1,ib) += shp[1][alpha]*DB[1][0] + shp[0][alpha]*DB[2][0];
	      K(ia+1,ib+1) += shp[1][alpha]*DB[1][1] + shp[0][alpha]*DB[2][1];
	      //	      matrixData[colIb   +   ia] += shp[0][alpha]*DB[0][0] + shp[1][alpha]*DB[2][0];
	      //matrixData[colIbP1 +   ia] += shp[0][alpha]*DB[0][1] + shp[1][alpha]*DB[2][1];
	      //matrixData[colIb   + ia+1] += shp[1][alpha]*DB[1][0] + shp[0][alpha]*DB[2][0];
	      //matrixData[colIbP1 + ia+1] += shp[1][alpha]*DB[1][1] + shp[0][alpha]*DB[2][1];
	      
	    }
	  }
	}
	
	return K;
}


const Matrix&
FourNodeQuad::getInitialStiff()
{
  static Matrix D(3,3);
  if (Ki != 0)
    return *Ki;

  K.Zero();
  
  double dvol;
  double DB[3][2];
  
  // Loop over the integration points
  for (int i = 0; i < 4; i++) {
    
    // Determine Jacobian for this integration point
    dvol = this->shapeFunction(pts[i][0], pts[i][1]);
    dvol *= (thickness*wts[i]);
    
    // Get the material tangent
    D = theMaterial[i]->getInitialTangent();
    if(theDamping[i]) D *= theDamping[i]->getStiffnessMultiplier();

    double D00 = D(0,0); double D01 = D(0,1); double D02 = D(0,2);
    double D10 = D(1,0); double D11 = D(1,1); double D12 = D(1,2);
    double D20 = D(2,0); double D21 = D(2,1); double D22 = D(2,2);
    
    // Perform numerical integration
    //K = K + (B^ D * B) * intWt(i)*intWt(j) * detJ;
    //K.addMatrixTripleProduct(1.0, B, D, intWt(i)*intWt(j)*detJ);
    for (int beta = 0, ib = 0, colIb =0, colIbP1 = 8; 
	 beta < 4; 
	 beta++, ib += 2, colIb += 16, colIbP1 += 16) {
      
      for (int alpha = 0, ia = 0; alpha < 4; alpha++, ia += 2) {
	
	DB[0][0] = dvol * (D00 * shp[0][beta] + D02 * shp[1][beta]);
	DB[1][0] = dvol * (D10 * shp[0][beta] + D12 * shp[1][beta]);
	DB[2][0] = dvol * (D20 * shp[0][beta] + D22 * shp[1][beta]);
	DB[0][1] = dvol * (D01 * shp[1][beta] + D02 * shp[0][beta]);
	DB[1][1] = dvol * (D11 * shp[1][beta] + D12 * shp[0][beta]);
	DB[2][1] = dvol * (D21 * shp[1][beta] + D22 * shp[0][beta]);
	
	matrixData[colIb   +   ia] += shp[0][alpha]*DB[0][0] + shp[1][alpha]*DB[2][0];
	matrixData[colIbP1 +   ia] += shp[0][alpha]*DB[0][1] + shp[1][alpha]*DB[2][1];
	matrixData[colIb   + ia+1] += shp[1][alpha]*DB[1][0] + shp[0][alpha]*DB[2][0];
	matrixData[colIbP1 + ia+1] += shp[1][alpha]*DB[1][1] + shp[0][alpha]*DB[2][1];
      }
    }
  }

  Ki = new Matrix(K);
  return K;
}

const Matrix&
FourNodeQuad::getMass()
{
	K.Zero();

	int i;
	static double rhoi[4];
	double sum = 0.0;
	for (i = 0; i < 4; i++) {
	  if (rho == 0)
	    rhoi[i] = theMaterial[i]->getRho();
	  else
	    rhoi[i] = rho;	    
	  sum += rhoi[i];
	}

	if (sum == 0.0)
	  return K;

	double rhodvol, Nrho;

	// Compute a lumped mass matrix
	for (i = 0; i < 4; i++) {

		// Determine Jacobian for this integration point
		rhodvol = this->shapeFunction(pts[i][0], pts[i][1]);

		// Element plus material density ... MAY WANT TO REMOVE ELEMENT DENSITY
		rhodvol *= (rhoi[i]*thickness*wts[i]);

		for (int alpha = 0, ia = 0; alpha < 4; alpha++, ia++) {
			Nrho = shp[2][alpha]*rhodvol;
			K(ia,ia) += Nrho;
			ia++;
			K(ia,ia) += Nrho;
		}
	}

	return K;
}

void
FourNodeQuad::zeroLoad(void)
{
	Q.Zero();

	applyLoad = 0;

	appliedB[0] = 0.0;
	appliedB[1] = 0.0;

  	return;
}

int 
FourNodeQuad::addLoad(ElementalLoad *theLoad, double loadFactor)
{
	// Added option for applying body forces in load pattern: C.McGann, U.Washington
	int type;
	const Vector &data = theLoad->getData(type, loadFactor);

	if (type == LOAD_TAG_SelfWeight) {
		applyLoad = 1;
		appliedB[0] += loadFactor*data(0)*b[0];
		appliedB[1] += loadFactor*data(1)*b[1];
		return 0;
	} else {
		opserr << "FourNodeQuad::addLoad - load type unknown for ele with tag: " << this->getTag() << endln;
		return -1;
	} 

	return -1;
}

int 
FourNodeQuad::addInertiaLoadToUnbalance(const Vector &accel)
{
  int i;
  static double rhoi[4];
  double sum = 0.0;
  for (i = 0; i < 4; i++) {
    if (rho == 0)
      rhoi[i] = theMaterial[i]->getRho();
    else
      rhoi[i] = rho;
    sum += rhoi[i];
  }
  
  if (sum == 0.0)
    return 0;
  
  // Get R * accel from the nodes
  const Vector &Raccel1 = theNodes[0]->getRV(accel);
  const Vector &Raccel2 = theNodes[1]->getRV(accel);
  const Vector &Raccel3 = theNodes[2]->getRV(accel);
  const Vector &Raccel4 = theNodes[3]->getRV(accel);
  
  if (2 != Raccel1.Size() || 2 != Raccel2.Size() || 2 != Raccel3.Size() ||
      2 != Raccel4.Size()) {
    opserr << "FourNodeQuad::addInertiaLoadToUnbalance matrix and vector sizes are incompatible\n";
    return -1;
  }
  
  static double ra[8];
  
  ra[0] = Raccel1(0);
  ra[1] = Raccel1(1);
  ra[2] = Raccel2(0);
  ra[3] = Raccel2(1);
  ra[4] = Raccel3(0);
  ra[5] = Raccel3(1);
  ra[6] = Raccel4(0);
  ra[7] = Raccel4(1);
  
  // Compute mass matrix
  this->getMass();
  
  // Want to add ( - fact * M R * accel ) to unbalance
  // Take advantage of lumped mass matrix
  for (i = 0; i < 8; i++)
    Q(i) += -K(i,i)*ra[i];
  
  return 0;
}

const Vector&
FourNodeQuad::getResistingForce()
{
  static Vector sigma(3);
	P.Zero();

	double dvol;

	// Loop over the integration points
	for (int i = 0; i < 4; i++) {

		// Determine Jacobian for this integration point
		dvol = this->shapeFunction(pts[i][0], pts[i][1]);
		dvol *= (thickness*wts[i]);

		// Get material stress response
		sigma = theMaterial[i]->getStress();

    if (theDamping[i])
    {
      theDamping[i]->update(sigma);
      sigma += theDamping[i]->getDampingForce();
    }

		// Perform numerical integration on internal force
		//P = P + (B^ sigma) * intWt(i)*intWt(j) * detJ;
		//P.addMatrixTransposeVector(1.0, B, sigma, intWt(i)*intWt(j)*detJ);
		for (int alpha = 0, ia = 0; alpha < 4; alpha++, ia += 2) {
			
			P(ia) += dvol*(shp[0][alpha]*sigma(0) + shp[1][alpha]*sigma(2));
			
			P(ia+1) += dvol*(shp[1][alpha]*sigma(1) + shp[0][alpha]*sigma(2));

			// Subtract equiv. body forces from the nodes
			//P = P - (N^ b) * intWt(i)*intWt(j) * detJ;
			//P.addMatrixTransposeVector(1.0, N, b, -intWt(i)*intWt(j)*detJ);
			if (applyLoad == 0) {
				P(ia) -= dvol*(shp[2][alpha]*b[0]);
				P(ia+1) -= dvol*(shp[2][alpha]*b[1]);
			} else {
				P(ia) -= dvol*(shp[2][alpha]*appliedB[0]);
				P(ia+1) -= dvol*(shp[2][alpha]*appliedB[1]);
			}
		}
	}

	// Subtract pressure loading from resisting force
	if (pressure != 0.0) {
		//P = P - pressureLoad;
		P.addVector(1.0, pressureLoad, -1.0);
	}
	
	// Subtract other external nodal loads ... P_res = P_int - P_ext
	//P = P - Q;
	P.addVector(1.0, Q, -1.0);

	return P;
}

const Vector&
FourNodeQuad::getResistingForceIncInertia()
{
	int i;
	static double rhoi[4];
	double sum = 0.0;
	for (i = 0; i < 4; i++) {
	  if (rho == 0)
	    rhoi[i] = theMaterial[i]->getRho();
	  else
	    rhoi[i] = rho;
	  sum += rhoi[i];
	}

	// if no mass terms .. just add damping terms
	if (sum == 0.0) {
	  this->getResistingForce();

	  // add the damping forces if rayleigh damping
	  if (betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
	    P += this->getRayleighDampingForces();

	  return P;
	}

	const Vector &accel1 = theNodes[0]->getTrialAccel();
	const Vector &accel2 = theNodes[1]->getTrialAccel();
	const Vector &accel3 = theNodes[2]->getTrialAccel();
	const Vector &accel4 = theNodes[3]->getTrialAccel();
	
	static double a[8];

	a[0] = accel1(0);
	a[1] = accel1(1);
	a[2] = accel2(0);
	a[3] = accel2(1);
	a[4] = accel3(0);
	a[5] = accel3(1);
	a[6] = accel4(0);
	a[7] = accel4(1);

	// Compute the current resisting force
	this->getResistingForce();

	// Compute the mass matrix
	this->getMass();

	// Take advantage of lumped mass matrix
	for (i = 0; i < 8; i++)
		P(i) += K(i,i)*a[i];

	// add the damping forces if rayleigh damping
	if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
	  P += this->getRayleighDampingForces();

	return P;
}

int
FourNodeQuad::sendSelf(int commitTag, Channel &theChannel)
{
  int res = 0;
  
  // note: we don't check for dataTag == 0 for Element
  // objects as that is taken care of in a commit by the Domain
  // object - don't want to have to do the check if sending data
  int dataTag = this->getDbTag();
  
  // Quad packs its data into a Vector and sends this to theChannel
  // along with its dbTag and the commitTag passed in the arguments
  static Vector data(11);
  data(0) = this->getTag();
  data(1) = thickness;
  data(2) = b[0];
  data(3) = b[1];
  data(4) = pressure;

  data(5) = alphaM;
  data(6) = betaK;
  data(7) = betaK0;
  data(8) = betaKc;
  
  data(9) = 0;
  data(10) = 0;
  if (theDamping[0]) {
    data(9) = theDamping[0]->getClassTag();
    int dbTag = theDamping[0]->getDbTag();
    if (dbTag == 0) {
      dbTag = theChannel.getDbTag();
      if (dbTag != 0)
        for (int i = 0 ;  i < 4; i++)
	        theDamping[i]->setDbTag(dbTag);
	  }
    data(10) = dbTag;
  }

  res += theChannel.sendVector(dataTag, commitTag, data);
  if (res < 0) {
    opserr << "WARNING FourNodeQuad::sendSelf() - " << this->getTag() << " failed to send Vector\n";
    return res;
  }	      
  

  // Now quad sends the ids of its materials
  int matDbTag;
  
  static ID idData(12);
  
  int i;
  for (i = 0; i < 4; i++) {
    idData(i) = theMaterial[i]->getClassTag();
    matDbTag = theMaterial[i]->getDbTag();
    // NOTE: we do have to ensure that the material has a database
    // tag if we are sending to a database channel.
    if (matDbTag == 0) {
      matDbTag = theChannel.getDbTag();
			if (matDbTag != 0)
			  theMaterial[i]->setDbTag(matDbTag);
    }
    idData(i+4) = matDbTag;
  }
  
  idData(8) = connectedExternalNodes(0);
  idData(9) = connectedExternalNodes(1);
  idData(10) = connectedExternalNodes(2);
  idData(11) = connectedExternalNodes(3);

  res += theChannel.sendID(dataTag, commitTag, idData);
  if (res < 0) {
    opserr << "WARNING FourNodeQuad::sendSelf() - " << this->getTag() << " failed to send ID\n";
    return res;
  }

  // Finally, quad asks its material objects to send themselves
  for (i = 0; i < 4; i++) {
    res += theMaterial[i]->sendSelf(commitTag, theChannel);
    if (res < 0) {
      opserr << "WARNING FourNodeQuad::sendSelf() - " << this->getTag() << " failed to send its Material\n";
      return res;
    }
  }
  
  // Ask the Damping to send itself
  if (theDamping[0]) {
    for (int i = 0 ;  i < 4; i++) {
      res += theDamping[i]->sendSelf(commitTag, theChannel);
      if (res < 0) {
        opserr << "FourNodeQuad::sendSelf -- could not send Damping\n";
        return res;
      }
    }
  }

  return res;
}

int
FourNodeQuad::recvSelf(int commitTag, Channel &theChannel,
                       FEM_ObjectBroker &theBroker)
{
  int res = 0;
  
  int dataTag = this->getDbTag();

  // Quad creates a Vector, receives the Vector and then sets the 
  // internal data with the data in the Vector
  static Vector data(11);
  res += theChannel.recvVector(dataTag, commitTag, data);
  if (res < 0) {
    opserr << "WARNING FourNodeQuad::recvSelf() - failed to receive Vector\n";
    return res;
  }
  
  this->setTag((int)data(0));
  thickness = data(1);
  b[0] = data(2);
  b[1] = data(3);
  pressure = data(4);

  alphaM = data(5);
  betaK = data(6);
  betaK0 = data(7);
  betaKc = data(8);

  static ID idData(12);
  // Quad now receives the tags of its four external nodes
  res += theChannel.recvID(dataTag, commitTag, idData);
  if (res < 0) {
    opserr << "WARNING FourNodeQuad::recvSelf() - " << this->getTag() << " failed to receive ID\n";
    return res;
  }

  connectedExternalNodes(0) = idData(8);
  connectedExternalNodes(1) = idData(9);
  connectedExternalNodes(2) = idData(10);
  connectedExternalNodes(3) = idData(11);
  
  if (theMaterial == 0) {
    // Allocate new materials
    theMaterial = new NDMaterial *[4];
    if (theMaterial == 0) {
      opserr << "FourNodeQuad::recvSelf() - Could not allocate NDMaterial* array\n";
      return -1;
    }
    for (int i = 0; i < 4; i++) {
      int matClassTag = idData(i);
      int matDbTag = idData(i+4);
      // Allocate new material with the sent class tag
      theMaterial[i] = theBroker.getNewNDMaterial(matClassTag);
      if (theMaterial[i] == 0) {
	    opserr << "FourNodeQuad::recvSelf() - Broker could not create NDMaterial of class type " << matClassTag << endln;
	    return -1;
      }
      // Now receive materials into the newly allocated space
      theMaterial[i]->setDbTag(matDbTag);
      res += theMaterial[i]->recvSelf(commitTag, theChannel, theBroker);
      if (res < 0) {
        opserr << "FourNodeQuad::recvSelf() - material " << i << "failed to recv itself\n";
	    return res;
      }
    }
  }

  // materials exist , ensure materials of correct type and recvSelf on them
  else {
    for (int i = 0; i < 4; i++) {
      int matClassTag = idData(i);
      int matDbTag = idData(i+4);
      // Check that material is of the right type; if not,
      // delete it and create a new one of the right type
      if (theMaterial[i]->getClassTag() != matClassTag) {
	    delete theMaterial[i];
	    theMaterial[i] = theBroker.getNewNDMaterial(matClassTag);
	    if (theMaterial[i] == 0) {
          opserr << "FourNodeQuad::recvSelf() - material " << i << "failed to create\n";
	      return -1;
	    }
      }
      // Receive the material
      theMaterial[i]->setDbTag(matDbTag);
      res += theMaterial[i]->recvSelf(commitTag, theChannel, theBroker);
      if (res < 0) {
        opserr << "FourNodeQuad::recvSelf() - material " << i << "failed to recv itself\n";
	    return res;
      }
    }
  }
  
  int dmpTag = (int)data(9);
  if (dmpTag) {
    for (int i = 0 ;  i < 4; i++) {
      // Check if the Damping is null; if so, get a new one
      if (theDamping[i] == 0) {
        theDamping[i] = theBroker.getNewDamping(dmpTag);
        if (theDamping[i] == 0) {
          opserr << "FourNodeQuad::recvSelf -- could not get a Damping\n";
          exit(-1);
        }
      }
  
      // Check that the Damping is of the right type; if not, delete
      // the current one and get a new one of the right type
      if (theDamping[i]->getClassTag() != dmpTag) {
        delete theDamping[i];
        theDamping[i] = theBroker.getNewDamping(dmpTag);
        if (theDamping[i] == 0) {
          opserr << "FourNodeQuad::recvSelf -- could not get a Damping\n";
          exit(-1);
        }
      }
  
      // Now, receive the Damping
      theDamping[i]->setDbTag((int)data(10));
      res += theDamping[i]->recvSelf(commitTag, theChannel, theBroker);
      if (res < 0) {
        opserr << "FourNodeQuad::recvSelf -- could not receive Damping\n";
        return res;
      }
    }
  }
  else {
    for (int i = 0; i < 4; i++)
    {
      if (theDamping[i])
      {
        delete theDamping[i];
        theDamping[i] = 0;
      }
    }
  }
    
  return res;
}

void
FourNodeQuad::Print(OPS_Stream &s, int flag)
{
  if (flag == 2) {

    s << "#FourNodeQuad\n";
    
    int i;
    const int numNodes = 4;
    const int nstress = 3 ;
    
    for (i=0; i<numNodes; i++) {
      const Vector &nodeCrd = theNodes[i]->getCrds();
      // const Vector &nodeDisp = theNodes[i]->getDisp();
      s << "#NODE " << nodeCrd(0) << " " << nodeCrd(1) << " " << endln;
     }
    
    // spit out the section location & invoke print on the scetion
    const int numMaterials = 4;

    static Vector avgStress(nstress);
    static Vector avgStrain(nstress);
    avgStress.Zero();
    avgStrain.Zero();
    for (i=0; i<numMaterials; i++) {
      avgStress += theMaterial[i]->getStress();
      avgStrain += theMaterial[i]->getStrain();
    }
    avgStress /= numMaterials;
    avgStrain /= numMaterials;

    s << "#AVERAGE_STRESS ";
    for (i=0; i<nstress; i++)
      s << avgStress(i) << " " ;
    s << endln;

    s << "#AVERAGE_STRAIN ";
    for (i=0; i<nstress; i++)
      s << avgStrain(i) << " " ;
    s << endln;
  }
  
  if (flag == OPS_PRINT_CURRENTSTATE) {
    s << "\nFourNodeQuad, element id:  " << this->getTag() << endln;
	s << "\tConnected external nodes:  " << connectedExternalNodes;
	s << "\tthickness:  " << thickness << endln;
	s << "\tsurface pressure:  " << pressure << endln;
    s << "\tmass density:  " << rho << endln;
    s << "\tbody forces:  " << b[0] << " " << b[1] << endln;
	theMaterial[0]->Print(s,flag);
	s << "\tStress (xx yy xy)" << endln;
	for (int i = 0; i < 4; i++)
		s << "\t\tGauss point " << i+1 << ": " << theMaterial[i]->getStress();
  }
  
  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
      s << "\t\t\t{";
      s << "\"name\": " << this->getTag() << ", ";
      s << "\"type\": \"FourNodeQuad\", ";
      s << "\"nodes\": [" << connectedExternalNodes(0) << ", ";
      s << connectedExternalNodes(1) << ", ";
      s << connectedExternalNodes(2) << ", ";
      s << connectedExternalNodes(3) << "], ";
      s << "\"thickness\": " << thickness << ", ";
      s << "\"surfacePressure\": " << pressure << ", ";
      s << "\"masspervolume\": " << rho << ", ";
      s << "\"bodyForces\": [" << b[0] << ", " << b[1] << "], ";
      s << "\"material\": \"" << theMaterial[0]->getTag() << "\"}";
  }
}

int
FourNodeQuad::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode)
{
	// get the end point display coords
	static Vector v1(3);
	static Vector v2(3);
	static Vector v3(3);
	static Vector v4(3);
	theNodes[0]->getDisplayCrds(v1, fact, displayMode);
	theNodes[1]->getDisplayCrds(v2, fact, displayMode);
	theNodes[2]->getDisplayCrds(v3, fact, displayMode);
	theNodes[3]->getDisplayCrds(v4, fact, displayMode);

	// place values in coords matrix
	static Matrix coords(4, 3);
	for (int i = 0; i < 3; i++) {
		coords(0, i) = v1(i);
		coords(1, i) = v2(i);
		coords(2, i) = v3(i);
		coords(3, i) = v4(i);
	}

    // set the quantity to be displayed at the nodes;
    // if displayMode is 1 through 3 we will plot material stresses otherwise 0.0
	static Vector values(4);
	if (displayMode < 4 && displayMode > 0) {
		for (int i = 0; i < 4; i++) {
			const Vector& stress = theMaterial[i]->getStress();
			values(i) = stress(displayMode - 1);
		}
	}
	else {
		for (int i = 0; i < 4; i++)
			values(i) = 0.0;
	}

	// draw the polygon
	return theViewer.drawPolygon(coords, values, this->getTag());
}

Response*
FourNodeQuad::setResponse(const char **argv, int argc, 
			  OPS_Stream &output)
{
  Response *theResponse =0;

  output.tag("ElementOutput");
  output.attr("eleType","FourNodeQuad");
  output.attr("eleTag",this->getTag());
  output.attr("node1",connectedExternalNodes[0]);
  output.attr("node2",connectedExternalNodes[1]);
  output.attr("node3",connectedExternalNodes[2]);
  output.attr("node4",connectedExternalNodes[3]);

  char dataOut[10];
  if (strcmp(argv[0],"force") == 0 || strcmp(argv[0],"forces") == 0) {

    for (int i=1; i<=4; i++) {
      sprintf(dataOut,"P1_%d",i);
      output.tag("ResponseType",dataOut);
      sprintf(dataOut,"P2_%d",i);
      output.tag("ResponseType",dataOut);
    }
    
    theResponse =  new ElementResponse(this, 1, P);
  }   

  else if (strcmp(argv[0],"material") == 0 || strcmp(argv[0],"integrPoint") == 0) {
    
    int pointNum = atoi(argv[1]);
    if (pointNum > 0 && pointNum <= 4) {
      
      output.tag("GaussPoint");
      output.attr("number",pointNum);
      output.attr("eta",pts[pointNum-1][0]);
      output.attr("neta",pts[pointNum-1][1]);

      theResponse =  theMaterial[pointNum-1]->setResponse(&argv[2], argc-2, output);
      
      output.endTag();

    }
  }
  else if ((strcmp(argv[0],"stresses") ==0) || (strcmp(argv[0],"stress") ==0)) {
    for (int i=0; i<4; i++) {
      output.tag("GaussPoint");
      output.attr("number",i+1);
      output.attr("eta",pts[i][0]);
      output.attr("neta",pts[i][1]);

      output.tag("NdMaterialOutput");
      output.attr("classType", theMaterial[i]->getClassTag());
      output.attr("tag", theMaterial[i]->getTag());
      
      output.tag("ResponseType","sigma11");
      output.tag("ResponseType","sigma22");
      output.tag("ResponseType","sigma12");
      
      output.endTag(); // GaussPoint
      output.endTag(); // NdMaterialOutput
      }
    theResponse =  new ElementResponse(this, 3, Vector(12));
  }

  else if ((strcmp(argv[0],"stressesAtNodes") ==0) || (strcmp(argv[0],"stressAtNodes") ==0)) {
    for (int i=0; i<4; i++) { // nnodes
      output.tag("NodalPoint");
      output.attr("number",i+1);
      // output.attr("eta",pts[i][0]);
      // output.attr("neta",pts[i][1]);

      // output.tag("NdMaterialOutput");
      // output.attr("classType", theMaterial[i]->getClassTag());
      // output.attr("tag", theMaterial[i]->getTag());

      output.tag("ResponseType","sigma11");
      output.tag("ResponseType","sigma22");
      output.tag("ResponseType","sigma12");

      output.endTag(); // GaussPoint
      // output.endTag(); // NdMaterialOutput
      }
    theResponse =  new ElementResponse(this, 11, Vector(12)); // 3 * nnodes
  }

  else if ((strcmp(argv[0],"strain") ==0) || (strcmp(argv[0],"strains") ==0)) {
    for (int i=0; i<4; i++) {
      output.tag("GaussPoint");
      output.attr("number",i+1);
      output.attr("eta",pts[i][0]);
      output.attr("neta",pts[i][1]);

      output.tag("NdMaterialOutput");
      output.attr("classType", theMaterial[i]->getClassTag());
      output.attr("tag", theMaterial[i]->getTag());
      
      output.tag("ResponseType","eta11");
      output.tag("ResponseType","eta22");
      output.tag("ResponseType","eta12");
      
      output.endTag(); // GaussPoint
      output.endTag(); // NdMaterialOutput
      }
    theResponse =  new ElementResponse(this, 4, Vector(12));
  }

  else if (theDamping[0] && strcmp(argv[0],"dampingStresses") ==0) {
    for (int i=0; i<4; i++) {
      output.tag("GaussPoint");
      output.attr("number",i+1);
      output.attr("eta",pts[i][0]);
      output.attr("neta",pts[i][1]);

      output.tag("NdMaterialOutput");
      output.attr("classType", theMaterial[i]->getClassTag());
      output.attr("tag", theMaterial[i]->getTag());
      
      output.tag("ResponseType","sigma11");
      output.tag("ResponseType","sigma22");
      output.tag("ResponseType","sigma12");
      
      output.endTag(); // GaussPoint
      output.endTag(); // NdMaterialOutput
      }
    theResponse =  new ElementResponse(this, 5, Vector(12));
  }

  output.endTag(); // ElementOutput

  return theResponse;
}

int 
FourNodeQuad::getResponse(int responseID, Information &eleInfo)
{
  if (responseID == 1) {

    return eleInfo.setVector(this->getResistingForce());

  } else if (responseID == 3) {

    // Loop over the integration points
    static Vector stresses(12);
    int cnt = 0;
    for (int i = 0; i < 4; i++) {

      // Get material stress response
      const Vector &sigma = theMaterial[i]->getStress();
      stresses(cnt) = sigma(0);
      stresses(cnt+1) = sigma(1);
      stresses(cnt+2) = sigma(2);
      cnt += 3;
    }
    
    return eleInfo.setVector(stresses);
      
  } else if (responseID == 11) {

    // extrapolate stress from Gauss points to element nodes
    static Vector stressGP(12); // 3*nip
    static Vector stressAtNodes(12); // 3*nnodes
	stressAtNodes.Zero();
    int cnt = 0;
	// first get stress components (xx, yy, xy) at Gauss points
    for (int i = 0; i < 4; i++) { // nip
      // Get material stress response
      const Vector &sigma = theMaterial[i]->getStress();
      stressGP(cnt) = sigma(0);
      stressGP(cnt+1) = sigma(1);
      stressGP(cnt+2) = sigma(2);
      cnt += 3;
    }

	// [nnodes][nip]
	const double We[4][4] = {{1.8660254037844386, -0.5, 0.1339745962155614, -0.5},
							 {-0.5, 1.8660254037844386, -0.5, 0.1339745962155614},
							 {0.1339745962155614, -0.5, 1.8660254037844386, -0.5},
							 {-0.5, 0.1339745962155614, -0.5, 1.8660254037844386}};

	int p, l;
	for (int i = 0; i < 4; i++) { // nnodes
	  for (int k = 0; k < 3; k++) { // number of stress components
		p = 3*i + k;
		for (int j = 0; j < 4; j++) { // nip
		  l = 3*j + k;
		  stressAtNodes(p) += We[i][j] * stressGP(l);
		  // opserr << "stressAtNodes(" << p << ") = We[" << i << "][" << j << "] * stressGP(" << l << ") = " << We[i][j] << " * " << stressGP(l) << " = " << stressAtNodes(p) <<  "\n";
		}
	  }
	}

    return eleInfo.setVector(stressAtNodes);

  } else if (responseID == 4) {

    // Loop over the integration points
    static Vector stresses(12);
    int cnt = 0;
    for (int i = 0; i < 4; i++) {

      // Get material stress response
      const Vector &sigma = theMaterial[i]->getStrain();
      stresses(cnt) = sigma(0);
      stresses(cnt+1) = sigma(1);
      stresses(cnt+2) = sigma(2);
      cnt += 3;
    }

    return eleInfo.setVector(stresses);
	
  } else if (responseID == 5) {

    // Loop over the integration points
    static Vector stresses(12);
    int cnt = 0;
    for (int i = 0; i < 4; i++) {

      // Get material stress response
      const Vector &sigma = theDamping[i]->getDampingForce();
      stresses(cnt) = sigma(0);
      stresses(cnt+1) = sigma(1);
      stresses(cnt+2) = sigma(2);
      cnt += 3;
    }
    
    return eleInfo.setVector(stresses);
      
  } else

    return -1;
}

int
FourNodeQuad::setParameter(const char **argv, int argc, Parameter &param)
{
  if (argc < 1)
    return -1;

  int res = -1;

  // quad pressure loading
  if (strcmp(argv[0],"pressure") == 0) {
    return param.addObject(2, this);
  }
  // a material parameter
  else if ((strstr(argv[0],"material") != 0) && (strcmp(argv[0],"materialState") != 0)) {

    if (argc < 3)
      return -1;

    int pointNum = atoi(argv[1]);
    if (pointNum > 0 && pointNum <= 4)
      return theMaterial[pointNum-1]->setParameter(&argv[2], argc-2, param);
    else 
      return -1;
  }

  // otherwise it could be just a forall material parameter
  else {

    int matRes = res;
    for (int i=0; i<4; i++) {

      matRes =  theMaterial[i]->setParameter(argv, argc, param);

      if (matRes != -1)
	res = matRes;
    }
  }
  
  return res;
}
    
int
FourNodeQuad::updateParameter(int parameterID, Information &info)
{
	int res = -1;
		int matRes = res;
  switch (parameterID) {
    case -1:
      return -1;

	case 1:
		
		for (int i = 0; i<4; i++) {
		matRes = theMaterial[i]->updateParameter(parameterID, info);
		}
		if (matRes != -1) {
			res = matRes;
		}
		return res;
      
	case 2:
		pressure = info.theDouble;
		this->setPressureLoadAtNodes();	// update consistent nodal loads
		return 0;

	default: 
	  /*	  
	  if (parameterID >= 100) { // material parameter
	    int pointNum = parameterID/100;
	    if (pointNum > 0 && pointNum <= 4)
	      return theMaterial[pointNum-1]->updateParameter(parameterID-100*pointNum, info);
	    else
	      return -1;
	  } else // unknown
	  */
	    return -1;
  }
}

double FourNodeQuad::shapeFunction(double xi, double eta)
{
	const Vector &nd1Crds = theNodes[0]->getCrds();
	const Vector &nd2Crds = theNodes[1]->getCrds();
	const Vector &nd3Crds = theNodes[2]->getCrds();
	const Vector &nd4Crds = theNodes[3]->getCrds();

	double oneMinuseta = 1.0-eta;
	double onePluseta = 1.0+eta;
	double oneMinusxi = 1.0-xi;
	double onePlusxi = 1.0+xi;

	shp[2][0] = 0.25*oneMinusxi*oneMinuseta;	// N_1
	shp[2][1] = 0.25*onePlusxi*oneMinuseta;		// N_2
	shp[2][2] = 0.25*onePlusxi*onePluseta;		// N_3
	shp[2][3] = 0.25*oneMinusxi*onePluseta;		// N_4

	double J[2][2];

	J[0][0] = 0.25 * (-nd1Crds(0)*oneMinuseta + nd2Crds(0)*oneMinuseta +
				nd3Crds(0)*(onePluseta) - nd4Crds(0)*(onePluseta));

	J[0][1] = 0.25 * (-nd1Crds(0)*oneMinusxi - nd2Crds(0)*onePlusxi +
				nd3Crds(0)*onePlusxi + nd4Crds(0)*oneMinusxi);

	J[1][0] = 0.25 * (-nd1Crds(1)*oneMinuseta + nd2Crds(1)*oneMinuseta +
				nd3Crds(1)*onePluseta - nd4Crds(1)*onePluseta);

	J[1][1] = 0.25 * (-nd1Crds(1)*oneMinusxi - nd2Crds(1)*onePlusxi +
				nd3Crds(1)*onePlusxi + nd4Crds(1)*oneMinusxi);

	double detJ = J[0][0]*J[1][1] - J[0][1]*J[1][0];
	double oneOverdetJ = 1.0/detJ;
	double L[2][2];

	// L = inv(J)
	L[0][0] =  J[1][1]*oneOverdetJ;
	L[1][0] = -J[0][1]*oneOverdetJ;
	L[0][1] = -J[1][0]*oneOverdetJ;
	L[1][1] =  J[0][0]*oneOverdetJ;

    double L00 = 0.25*L[0][0];
    double L10 = 0.25*L[1][0];
    double L01 = 0.25*L[0][1];
    double L11 = 0.25*L[1][1];
	
	double L00oneMinuseta = L00*oneMinuseta;
	double L00onePluseta  = L00*onePluseta;
	double L01oneMinusxi  = L01*oneMinusxi;
	double L01onePlusxi   = L01*onePlusxi;

	double L10oneMinuseta = L10*oneMinuseta;
	double L10onePluseta  = L10*onePluseta;
	double L11oneMinusxi  = L11*oneMinusxi;
	double L11onePlusxi   = L11*onePlusxi;

	// See Cook, Malkus, Plesha p. 169 for the derivation of these terms
    shp[0][0] = -L00oneMinuseta - L01oneMinusxi;	// N_1,1
    shp[0][1] =  L00oneMinuseta - L01onePlusxi;		// N_2,1
    shp[0][2] =  L00onePluseta  + L01onePlusxi;		// N_3,1
    shp[0][3] = -L00onePluseta  + L01oneMinusxi;	// N_4,1
	
    shp[1][0] = -L10oneMinuseta - L11oneMinusxi;	// N_1,2
    shp[1][1] =  L10oneMinuseta - L11onePlusxi;		// N_2,2
    shp[1][2] =  L10onePluseta  + L11onePlusxi;		// N_3,2
    shp[1][3] = -L10onePluseta  + L11oneMinusxi;	// N_4,2

    return detJ;
}

void 
FourNodeQuad::setPressureLoadAtNodes(void)
{
        pressureLoad.Zero();

	if (pressure == 0.0)
		return;

	const Vector &node1 = theNodes[0]->getCrds();
	const Vector &node2 = theNodes[1]->getCrds();
	const Vector &node3 = theNodes[2]->getCrds();
	const Vector &node4 = theNodes[3]->getCrds();

	double x1 = node1(0);
	double y1 = node1(1);
	double x2 = node2(0);
	double y2 = node2(1);
	double x3 = node3(0);
	double y3 = node3(1);
	double x4 = node4(0);
	double y4 = node4(1);

	double dx12 = x2-x1;
	double dy12 = y2-y1;
	double dx23 = x3-x2;
	double dy23 = y3-y2;
	double dx34 = x4-x3;
	double dy34 = y4-y3;
	double dx41 = x1-x4;
	double dy41 = y1-y4;

	double pressureOver2 = pressure/2.0;

	// Contribution from side 12
	pressureLoad(0) += pressureOver2*dy12;
	pressureLoad(2) += pressureOver2*dy12;
	pressureLoad(1) += pressureOver2*-dx12;
	pressureLoad(3) += pressureOver2*-dx12;

	// Contribution from side 23
	pressureLoad(2) += pressureOver2*dy23;
	pressureLoad(4) += pressureOver2*dy23;
	pressureLoad(3) += pressureOver2*-dx23;
	pressureLoad(5) += pressureOver2*-dx23;

	// Contribution from side 34
	pressureLoad(4) += pressureOver2*dy34;
	pressureLoad(6) += pressureOver2*dy34;
	pressureLoad(5) += pressureOver2*-dx34;
	pressureLoad(7) += pressureOver2*-dx34;

	// Contribution from side 41
	pressureLoad(6) += pressureOver2*dy41;
	pressureLoad(0) += pressureOver2*dy41;
	pressureLoad(7) += pressureOver2*-dx41;
	pressureLoad(1) += pressureOver2*-dx41;

	//pressureLoad = pressureLoad*thickness;
}








