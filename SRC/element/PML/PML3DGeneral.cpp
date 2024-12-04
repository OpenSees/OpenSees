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
#include <stdio.h> 
#include <stdlib.h> 
#include <cmath> 

#include <ID.h> 
#include <Vector.h>
#include <Matrix.h>
#include <Element.h>
#include <Node.h>
#include <Domain.h>
#include <ErrorHandler.h>
#include "PML3DGeneral.h"
#include <Renderer.h>
#include <ElementResponse.h>
#include <Parameter.h>
#include <ElementalLoad.h>

#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <elementAPI.h>

#include <fstream>
#include <Eigen/Dense>



// =======================================================================
// PML3DGeneral element tcl command
// =======================================================================
void* OPS_PML3DGeneral()
{
	// check if the total number of arguments passed is correct
	if (OPS_GetNumRemainingInputArgs() < (15)) {
		opserr << "WARNING insufficient arguments\n";
		opserr << "Want: element PML3DGeneral eleTag? [8 integer nodeTags] $matTag? gamma? beta? eta? m_pml? L_pml? R_pml? x0_pml? y0_pml? z0_pml? nx_pml? ny_pml? nz_pml? \n";
		return 0;
	}

	// reading element tag and node numbers and matTag
	int idata[9];
	int num = 9;
	if (OPS_GetIntInput(&num, idata) < 0) {
		opserr << "WARNING: invalid integer data : could be the tag or the node numbers \n";
		return 0;
	}





    //reading E , nu, rho
    double Elasticity[3];
    num = 3;
    if (OPS_GetDoubleInput(&num, Elasticity) < 0) {
        opserr << "WARNING: invalid double data: could be Elasticity, Poisson's ratio and density\n";
        return 0;
    }

    double _E = Elasticity[0];
    double _nu = Elasticity[1];
    double _rho = Elasticity[2];

	// reading Newmark parameters
	double Newmark[3];
	num = 3;
	if (OPS_GetDoubleInput(&num, Newmark) < 0) {
		opserr << "WARNING: invalid double data: could be Newmark parameters\n";
		return 0;
	}

	// reading material properties
	double dData[9]; 
    num = 9;
	if (OPS_GetDoubleInput(&num, dData) < 0) {
		opserr << "WARNING: invalid PML Material Propertises\n";
		return 0;
	}

	// create a new PML3DGeneral element and add it to the Domain
    double _gamma = Newmark[0];
    double _beta  = Newmark[1];
    double _eta   = Newmark[2]; 

    // double _E     = dData[0];
    // double _nu    = dData[1];
    // double _rho   = dData[2];

    double _m_pml = dData[0];
    double _L_pml = dData[1];
    double _R_pml = dData[2];

    double _x0_pml = dData[3];
    double _y0_pml = dData[4];
    double _z0_pml = dData[5];

    double _nx_pml = dData[6];
    double _ny_pml = dData[7];
    double _nz_pml = dData[8];

	return new PML3DGeneral(idata[0], &idata[1],
                            _E, _nu, _rho,
                            _gamma, _beta, _eta,  
                            _m_pml, _L_pml, _R_pml, 
                            _x0_pml, _y0_pml, _z0_pml, 
                            _nx_pml, _ny_pml, _nz_pml);
}

// =======================================================================
// Static data
// =======================================================================
Vector  PML3DGeneral::resid(PML3DGeneral_NUM_DOF);
double  PML3DGeneral::eta = 0.;
double  PML3DGeneral::beta = 0.;
double  PML3DGeneral::gamma = 0.;
double  PML3DGeneral::dt = 0.;
int     PML3DGeneral::eleCount = 0;
// =======================================================================
// Null constructor
// =======================================================================
PML3DGeneral::PML3DGeneral() 
    :Element(0, ELE_TAG_PML3DGeneral),
	connectedExternalNodes(PML3DGeneral_NUM_NODES),
    ubar(PML3DGeneral_NUM_DOF), ubart(PML3DGeneral_NUM_DOF),
    StiffnessMatrix(PML3DGeneral_NUM_DOF, PML3DGeneral_NUM_DOF),
    ImpedanceMatrix(PML3DGeneral_NUM_DOF, PML3DGeneral_NUM_DOF),
    MassMatrix(PML3DGeneral_NUM_DOF, PML3DGeneral_NUM_DOF),
    DampingMatrix(PML3DGeneral_NUM_DOF, PML3DGeneral_NUM_DOF),
    EffectiveStiffnessMatrix(PML3DGeneral_NUM_DOF, PML3DGeneral_NUM_DOF),
    m_pml(0), L_pml(0), R_pml(0),
    x0_pml(0), y0_pml(0), z0_pml(0),
    nx_pml(0), ny_pml(0), nz_pml(0),
    E(0), nu(0), rho(0),
    Domainptr(0)
{
    // initialize node pointers
    for (int i = 0; i < PML3DGeneral_NUM_NODES; i++) {
        connectedExternalNodes(i) = 0;
        nodePointers[i] = 0;
    }

    // initialize static newmark parameters
    eta = 0;
    beta = 0;
    gamma = 0;

    ubar.Zero();
    ubart.Zero();
    StiffnessMatrix.Zero();
    MassMatrix.Zero();
    DampingMatrix.Zero();
    EffectiveStiffnessMatrix.Zero();
    ImpedanceMatrix.Zero();


}
// =======================================================================
// Full constructor
// =======================================================================
PML3DGeneral::PML3DGeneral(int tag, int* nodeTags,
                                double _E, double _nu, double _rho,
                                double _gamma, double _beta, double _eta,
                                double _m_pml, double _L_pml, double _R_pml, 
                                double _x0_pml, double _y0_pml, double _z0_pml,
                                double _nx_pml, double _ny_pml, double _nz_pml)
                                :Element(tag, ELE_TAG_PML3DGeneral),
                                connectedExternalNodes(PML3DGeneral_NUM_NODES),
                                ubar(PML3DGeneral_NUM_DOF), ubart(PML3DGeneral_NUM_DOF),
                                StiffnessMatrix(PML3DGeneral_NUM_DOF, PML3DGeneral_NUM_DOF),
                                ImpedanceMatrix(PML3DGeneral_NUM_DOF, PML3DGeneral_NUM_DOF),
                                MassMatrix(PML3DGeneral_NUM_DOF, PML3DGeneral_NUM_DOF),
                                DampingMatrix(PML3DGeneral_NUM_DOF, PML3DGeneral_NUM_DOF),
                                EffectiveStiffnessMatrix(PML3DGeneral_NUM_DOF, PML3DGeneral_NUM_DOF),
                                m_pml(_m_pml), L_pml(_L_pml), R_pml(_R_pml),
                                x0_pml(_x0_pml), y0_pml(_y0_pml), z0_pml(_z0_pml),
                                nx_pml(_nx_pml), ny_pml(_ny_pml), nz_pml(_nz_pml),
                                E(_E), nu(_nu), rho(_rho),
                                Domainptr(0)
{
    eleCount++;
    if (eleCount == 1) {
        opserr << "PML3DGeneral::PML3DGeneral()" << endln;
    }
    // initialize node pointers
    for (int i = 0; i < PML3DGeneral_NUM_NODES; i++) {
        connectedExternalNodes(i) = nodeTags[i];
        nodePointers[i] = 0;
    }

    // initialize static newmark parameters
    eta = _eta;
    beta = _beta;
    gamma = _gamma;

    // initialize vectors and matrices
	ubart.Zero();
	ubar.Zero();
    StiffnessMatrix.Zero();
    MassMatrix.Zero();
    DampingMatrix.Zero();
    EffectiveStiffnessMatrix.Zero();
    ImpedanceMatrix.Zero();

    // create the coordinate vectors
    opserr << "PML3DGeneral::PML3DGeneral()" << endln;
    opserr << "element data: "    << endln;
    opserr << "\t tag: "    << tag    << endln;
    opserr << "\t nodeTags: "     << connectedExternalNodes;
    opserr << "\t E: "     << E     << endln;
    opserr << "\t nu: "    << nu    << endln;
    opserr << "\t rho: "   << rho   << endln;
    opserr << "\t eta: "    << eta    << endln;
    opserr << "\t beta: "   << beta   << endln;
    opserr << "\t gamma: "  << gamma  << endln;
    opserr << "\t m_pml: "  << m_pml  << endln;
    opserr << "\t L_pml: "  << L_pml  << endln;
    opserr << "\t R_pml: "  << R_pml  << endln;
    opserr << "\t x0_pml: " << x0_pml << endln;
    opserr << "\t y0_pml: " << y0_pml << endln;
    opserr << "\t z0_pml: " << z0_pml << endln;
    opserr << "\t nx_pml: " << nx_pml << endln;
    opserr << "\t ny_pml: " << ny_pml << endln;
    opserr << "\t nz_pml: " << nz_pml << endln;
}

// =======================================================================
// Destructor
// =======================================================================
PML3DGeneral::~PML3DGeneral()
{
    // does nothing
}

// =======================================================================
// Set domain
// =======================================================================
void PML3DGeneral::setDomain(Domain* theDomain)
{

    if (theDomain == 0) {
        opserr << "PML3DGeneral::setDomain() - FATAL error setting domain: " << this->getTag() << endln;
        exit(-1);
    }
    Domainptr = theDomain;
    // node pointers
    for (int i = 0; i < PML3DGeneral_NUM_NODES; i++) {
        nodePointers[i] = theDomain->getNode(connectedExternalNodes(i));
    }
    this->DomainComponent::setDomain(theDomain);


    // compute K 
    this->ComputeStiffnessMatrix();
    // compute M
    this->ComputeMassMatrix();
    // compute C
    this->ComputePMLMatrix();
    // compute G
    this->ComputePMLMatrix();


    Vector X1 = nodePointers[0]->getCrds();
    Vector X2 = nodePointers[1]->getCrds();
    Vector X3 = nodePointers[2]->getCrds();
    Vector X4 = nodePointers[3]->getCrds();
    Vector X5 = nodePointers[4]->getCrds();
    Vector X6 = nodePointers[5]->getCrds();
    Vector X7 = nodePointers[6]->getCrds();
    Vector X8 = nodePointers[7]->getCrds();

    opserr << "X1: " << X1(0) << " " << X1(1) << " " << X1(2) << endln;
    opserr << "X2: " << X2(0) << " " << X2(1) << " " << X2(2) << endln;
    opserr << "X3: " << X3(0) << " " << X3(1) << " " << X3(2) << endln;
    opserr << "X4: " << X4(0) << " " << X4(1) << " " << X4(2) << endln;
    opserr << "X5: " << X5(0) << " " << X5(1) << " " << X5(2) << endln;
    opserr << "X6: " << X6(0) << " " << X6(1) << " " << X6(2) << endln;
    opserr << "X7: " << X7(0) << " " << X7(1) << " " << X7(2) << endln;
    opserr << "X8: " << X8(0) << " " << X8(1) << " " << X8(2) << endln;

}

// =======================================================================
// Computes the jacobian of the transformation. 
// =======================================================================
Eigen::MatrixXd 
PML3DGeneral::ComputeJacobianMatrix(const double ri, const double si, const double ti) const{
    //Gets the element coordinates in deformed configuration.
    Vector X1 = nodePointers[0]->getCrds();
    Vector X2 = nodePointers[1]->getCrds();
    Vector X3 = nodePointers[2]->getCrds();
    Vector X4 = nodePointers[3]->getCrds();
    Vector X5 = nodePointers[4]->getCrds();
    Vector X6 = nodePointers[5]->getCrds();
    Vector X7 = nodePointers[6]->getCrds();
    Vector X8 = nodePointers[7]->getCrds();

	//Jacobian coefficients:
	double J11 = -1.0/8.0*(1.0 - si)*(1.0 - ti)*X1(0) + 1.0/8.0*(1.0 - si)*(1.0 - ti)*X2(0) + 1.0/8.0*(1.0 + si)*(1.0 - ti)*X3(0) - 1.0/8.0*(1.0 + si)*(1.0 - ti)*X4(0) - 1.0/8.0*(1.0 - si)*(1.0 + ti)*X5(0) + 1.0/8.0*(1.0 - si)*(1.0 + ti)*X6(0) + 1.0/8.0*(1.0 + si)*(1.0 + ti)*X7(0) - 1.0/8.0*(1.0 + si)*(1.0 + ti)*X8(0);
	double J12 = -1.0/8.0*(1.0 - si)*(1.0 - ti)*X1(1) + 1.0/8.0*(1.0 - si)*(1.0 - ti)*X2(1) + 1.0/8.0*(1.0 + si)*(1.0 - ti)*X3(1) - 1.0/8.0*(1.0 + si)*(1.0 - ti)*X4(1) - 1.0/8.0*(1.0 - si)*(1.0 + ti)*X5(1) + 1.0/8.0*(1.0 - si)*(1.0 + ti)*X6(1) + 1.0/8.0*(1.0 + si)*(1.0 + ti)*X7(1) - 1.0/8.0*(1.0 + si)*(1.0 + ti)*X8(1);
	double J13 = -1.0/8.0*(1.0 - si)*(1.0 - ti)*X1(2) + 1.0/8.0*(1.0 - si)*(1.0 - ti)*X2(2) + 1.0/8.0*(1.0 + si)*(1.0 - ti)*X3(2) - 1.0/8.0*(1.0 + si)*(1.0 - ti)*X4(2) - 1.0/8.0*(1.0 - si)*(1.0 + ti)*X5(2) + 1.0/8.0*(1.0 - si)*(1.0 + ti)*X6(2) + 1.0/8.0*(1.0 + si)*(1.0 + ti)*X7(2) - 1.0/8.0*(1.0 + si)*(1.0 + ti)*X8(2); 
	double J21 = -1.0/8.0*(1.0 - ri)*(1.0 - ti)*X1(0) - 1.0/8.0*(1.0 + ri)*(1.0 - ti)*X2(0) + 1.0/8.0*(1.0 + ri)*(1.0 - ti)*X3(0) + 1.0/8.0*(1.0 - ri)*(1.0 - ti)*X4(0) - 1.0/8.0*(1.0 - ri)*(1.0 + ti)*X5(0) - 1.0/8.0*(1.0 + ri)*(1.0 + ti)*X6(0) + 1.0/8.0*(1.0 + ri)*(1.0 + ti)*X7(0) + 1.0/8.0*(1.0 - ri)*(1.0 + ti)*X8(0);
	double J22 = -1.0/8.0*(1.0 - ri)*(1.0 - ti)*X1(1) - 1.0/8.0*(1.0 + ri)*(1.0 - ti)*X2(1) + 1.0/8.0*(1.0 + ri)*(1.0 - ti)*X3(1) + 1.0/8.0*(1.0 - ri)*(1.0 - ti)*X4(1) - 1.0/8.0*(1.0 - ri)*(1.0 + ti)*X5(1) - 1.0/8.0*(1.0 + ri)*(1.0 + ti)*X6(1) + 1.0/8.0*(1.0 + ri)*(1.0 + ti)*X7(1) + 1.0/8.0*(1.0 - ri)*(1.0 + ti)*X8(1);
	double J23 = -1.0/8.0*(1.0 - ri)*(1.0 - ti)*X1(2) - 1.0/8.0*(1.0 + ri)*(1.0 - ti)*X2(2) + 1.0/8.0*(1.0 + ri)*(1.0 - ti)*X3(2) + 1.0/8.0*(1.0 - ri)*(1.0 - ti)*X4(2) - 1.0/8.0*(1.0 - ri)*(1.0 + ti)*X5(2) - 1.0/8.0*(1.0 + ri)*(1.0 + ti)*X6(2) + 1.0/8.0*(1.0 + ri)*(1.0 + ti)*X7(2) + 1.0/8.0*(1.0 - ri)*(1.0 + ti)*X8(2);
	double J31 = -1.0/8.0*(1.0 - ri)*(1.0 - si)*X1(0) - 1.0/8.0*(1.0 + ri)*(1.0 - si)*X2(0) - 1.0/8.0*(1.0 + ri)*(1.0 + si)*X3(0) - 1.0/8.0*(1.0 - ri)*(1.0 + si)*X4(0) + 1.0/8.0*(1.0 - ri)*(1.0 - si)*X5(0) + 1.0/8.0*(1.0 + ri)*(1.0 - si)*X6(0) + 1.0/8.0*(1.0 + ri)*(1.0 + si)*X7(0) + 1.0/8.0*(1.0 - ri)*(1.0 + si)*X8(0);
	double J32 = -1.0/8.0*(1.0 - ri)*(1.0 - si)*X1(1) - 1.0/8.0*(1.0 + ri)*(1.0 - si)*X2(1) - 1.0/8.0*(1.0 + ri)*(1.0 + si)*X3(1) - 1.0/8.0*(1.0 - ri)*(1.0 + si)*X4(1) + 1.0/8.0*(1.0 - ri)*(1.0 - si)*X5(1) + 1.0/8.0*(1.0 + ri)*(1.0 - si)*X6(1) + 1.0/8.0*(1.0 + ri)*(1.0 + si)*X7(1) + 1.0/8.0*(1.0 - ri)*(1.0 + si)*X8(1);
	double J33 = -1.0/8.0*(1.0 - ri)*(1.0 - si)*X1(2) - 1.0/8.0*(1.0 + ri)*(1.0 - si)*X2(2) - 1.0/8.0*(1.0 + ri)*(1.0 + si)*X3(2) - 1.0/8.0*(1.0 - ri)*(1.0 + si)*X4(2) + 1.0/8.0*(1.0 - ri)*(1.0 - si)*X5(2) + 1.0/8.0*(1.0 + ri)*(1.0 - si)*X6(2) + 1.0/8.0*(1.0 + ri)*(1.0 + si)*X7(2) + 1.0/8.0*(1.0 - ri)*(1.0 + si)*X8(2);

	//Jacobia Matrix definition:
	Eigen::MatrixXd Jij(3,3);
	Jij << J11, J12, J13,
		   J21, J22, J23,
		   J31, J32, J33;

	return Jij;
}

// =======================================================================
// Compute Shape Function at Gauss Point
// =======================================================================
Eigen::MatrixXd 
PML3DGeneral::ComputeShapeFunctionMatrix(const double ri, const double si, const double ti) const{
	//Shape function coefficients:
	double H11 = 1.0/8.0*(1.0 - ri)*(1.0 - si)*(1.0 - ti);
	double H22 = 1.0/8.0*(1.0 + ri)*(1.0 - si)*(1.0 - ti);
	double H33 = 1.0/8.0*(1.0 + ri)*(1.0 + si)*(1.0 - ti);
	double H44 = 1.0/8.0*(1.0 - ri)*(1.0 + si)*(1.0 - ti);
	double H55 = 1.0/8.0*(1.0 - ri)*(1.0 - si)*(1.0 + ti);
	double H66 = 1.0/8.0*(1.0 + ri)*(1.0 - si)*(1.0 + ti);
	double H77 = 1.0/8.0*(1.0 + ri)*(1.0 + si)*(1.0 + ti);
	double H88 = 1.0/8.0*(1.0 - ri)*(1.0 + si)*(1.0 + ti);

	//Shape function matrix:
	Eigen::MatrixXd Hij(3,24);
	Hij << H11, 0.0, 0.0, H22, 0.0, 0.0, H33, 0.0, 0.0, H44, 0.0, 0.0, H55, 0.0, 0.0, H66, 0.0, 0.0, H77, 0.0, 0.0, H88, 0.0, 0.0,
		   0.0, H11, 0.0, 0.0, H22, 0.0, 0.0, H33, 0.0, 0.0, H44, 0.0, 0.0, H55, 0.0, 0.0, H66, 0.0, 0.0, H77, 0.0, 0.0, H88, 0.0,
		   0.0, 0.0, H11, 0.0, 0.0, H22, 0.0, 0.0, H33, 0.0, 0.0, H44, 0.0, 0.0, H55, 0.0, 0.0, H66, 0.0, 0.0, H77, 0.0, 0.0, H88;

	return Hij;
}

// =======================================================================
// Evaluates the lumped-mass matrix matrix at a given Gauss point.
// =======================================================================
Eigen::MatrixXd 
PML3DGeneral::ComputeStrainDisplacementMatrix(const double ri, const double si, const double ti, const Eigen::MatrixXd &Jij) const{
	//Inverse jacobian matrix:
	Eigen::MatrixXd J = Jij.inverse();

	//Strain-displacement matrix coefficients:
	double B11 = -1.0/8.0*J(0,2)*(1.0 - ri)*(1.0 - si) - 1.0/8.0*J(0,1)*(1.0 - ri)*(1.0 - ti) - 1.0/8.0*J(0,0)*(1.0 - si)*(1.0 - ti);
	double B21 = -1.0/8.0*J(0,2)*(1.0 + ri)*(1.0 - si) - 1.0/8.0*J(0,1)*(1.0 + ri)*(1.0 - ti) + 1.0/8.0*J(0,0)*(1.0 - si)*(1.0 - ti);
	double B31 = -1.0/8.0*J(0,2)*(1.0 + ri)*(1.0 + si) + 1.0/8.0*J(0,1)*(1.0 + ri)*(1.0 - ti) + 1.0/8.0*J(0,0)*(1.0 + si)*(1.0 - ti);
	double B41 = -1.0/8.0*J(0,2)*(1.0 - ri)*(1.0 + si) + 1.0/8.0*J(0,1)*(1.0 - ri)*(1.0 - ti) - 1.0/8.0*J(0,0)*(1.0 + si)*(1.0 - ti);
	double B51 =  1.0/8.0*J(0,2)*(1.0 - ri)*(1.0 - si) - 1.0/8.0*J(0,1)*(1.0 - ri)*(1.0 + ti) - 1.0/8.0*J(0,0)*(1.0 - si)*(1.0 + ti);
	double B61 =  1.0/8.0*J(0,2)*(1.0 + ri)*(1.0 - si) - 1.0/8.0*J(0,1)*(1.0 + ri)*(1.0 + ti) + 1.0/8.0*J(0,0)*(1.0 - si)*(1.0 + ti);
	double B71 =  1.0/8.0*J(0,2)*(1.0 + ri)*(1.0 + si) + 1.0/8.0*J(0,1)*(1.0 + ri)*(1.0 + ti) + 1.0/8.0*J(0,0)*(1.0 + si)*(1.0 + ti);
	double B81 =  1.0/8.0*J(0,2)*(1.0 - ri)*(1.0 + si) + 1.0/8.0*J(0,1)*(1.0 - ri)*(1.0 + ti) - 1.0/8.0*J(0,0)*(1.0 + si)*(1.0 + ti);

	double B12 = -1.0/8.0*J(1,2)*(1.0 - ri)*(1.0 - si) - 1.0/8.0*J(1,1)*(1.0 - ri)*(1.0 - ti) - 1.0/8.0*J(1,0)*(1.0 - si)*(1.0 - ti);
	double B22 = -1.0/8.0*J(1,2)*(1.0 + ri)*(1.0 - si) - 1.0/8.0*J(1,1)*(1.0 + ri)*(1.0 - ti) + 1.0/8.0*J(1,0)*(1.0 - si)*(1.0 - ti);
	double B32 = -1.0/8.0*J(1,2)*(1.0 + ri)*(1.0 + si) + 1.0/8.0*J(1,1)*(1.0 + ri)*(1.0 - ti) + 1.0/8.0*J(1,0)*(1.0 + si)*(1.0 - ti);
	double B42 = -1.0/8.0*J(1,2)*(1.0 - ri)*(1.0 + si) + 1.0/8.0*J(1,1)*(1.0 - ri)*(1.0 - ti) - 1.0/8.0*J(1,0)*(1.0 + si)*(1.0 - ti);
	double B52 =  1.0/8.0*J(1,2)*(1.0 - ri)*(1.0 - si) - 1.0/8.0*J(1,1)*(1.0 - ri)*(1.0 + ti) - 1.0/8.0*J(1,0)*(1.0 - si)*(1.0 + ti);
	double B62 =  1.0/8.0*J(1,2)*(1.0 + ri)*(1.0 - si) - 1.0/8.0*J(1,1)*(1.0 + ri)*(1.0 + ti) + 1.0/8.0*J(1,0)*(1.0 - si)*(1.0 + ti);
	double B72 =  1.0/8.0*J(1,2)*(1.0 + ri)*(1.0 + si) + 1.0/8.0*J(1,1)*(1.0 + ri)*(1.0 + ti) + 1.0/8.0*J(1,0)*(1.0 + si)*(1.0 + ti);
	double B82 =  1.0/8.0*J(1,2)*(1.0 - ri)*(1.0 + si) + 1.0/8.0*J(1,1)*(1.0 - ri)*(1.0 + ti) - 1.0/8.0*J(1,0)*(1.0 + si)*(1.0 + ti);

	double B13 = -1.0/8.0*J(2,2)*(1.0 - ri)*(1.0 - si) - 1.0/8.0*J(2,1)*(1.0 - ri)*(1.0 - ti) - 1.0/8.0*J(2,0)*(1.0 - si)*(1.0 - ti);
	double B23 = -1.0/8.0*J(2,2)*(1.0 + ri)*(1.0 - si) - 1.0/8.0*J(2,1)*(1.0 + ri)*(1.0 - ti) + 1.0/8.0*J(2,0)*(1.0 - si)*(1.0 - ti);
	double B33 = -1.0/8.0*J(2,2)*(1.0 + ri)*(1.0 + si) + 1.0/8.0*J(2,1)*(1.0 + ri)*(1.0 - ti) + 1.0/8.0*J(2,0)*(1.0 + si)*(1.0 - ti);
	double B43 = -1.0/8.0*J(2,2)*(1.0 - ri)*(1.0 + si) + 1.0/8.0*J(2,1)*(1.0 - ri)*(1.0 - ti) - 1.0/8.0*J(2,0)*(1.0 + si)*(1.0 - ti);
	double B53 =  1.0/8.0*J(2,2)*(1.0 - ri)*(1.0 - si) - 1.0/8.0*J(2,1)*(1.0 - ri)*(1.0 + ti) - 1.0/8.0*J(2,0)*(1.0 - si)*(1.0 + ti);
	double B63 =  1.0/8.0*J(2,2)*(1.0 + ri)*(1.0 - si) - 1.0/8.0*J(2,1)*(1.0 + ri)*(1.0 + ti) + 1.0/8.0*J(2,0)*(1.0 - si)*(1.0 + ti);
	double B73 =  1.0/8.0*J(2,2)*(1.0 + ri)*(1.0 + si) + 1.0/8.0*J(2,1)*(1.0 + ri)*(1.0 + ti) + 1.0/8.0*J(2,0)*(1.0 + si)*(1.0 + ti);
	double B83 =  1.0/8.0*J(2,2)*(1.0 - ri)*(1.0 + si) + 1.0/8.0*J(2,1)*(1.0 - ri)*(1.0 + ti) - 1.0/8.0*J(2,0)*(1.0 + si)*(1.0 + ti);

	//Deformation matrix definition:
	Eigen::MatrixXd Bij(6,24);
	Bij <<  B11, 0.0, 0.0, B21, 0.0, 0.0, B31, 0.0, 0.0, B41, 0.0, 0.0, B51, 0.0, 0.0, B61, 0.0, 0.0, B71, 0.0, 0.0, B81, 0.0, 0.0,
		    0.0, B12, 0.0, 0.0, B22, 0.0, 0.0, B32, 0.0, 0.0, B42, 0.0, 0.0, B52, 0.0, 0.0, B62, 0.0, 0.0, B72, 0.0, 0.0, B82, 0.0,
		    0.0, 0.0, B13, 0.0, 0.0, B23, 0.0, 0.0, B33, 0.0, 0.0, B43, 0.0, 0.0, B53, 0.0, 0.0, B63, 0.0, 0.0, B73, 0.0, 0.0, B83,
		    B12, B11, 0.0, B22, B21, 0.0, B32, B31, 0.0, B42, B41, 0.0, B52, B51, 0.0, B62, B61, 0.0, B72, B71, 0.0, B82, B81, 0.0,
		    0.0, B13, B12, 0.0, B23, B22, 0.0, B33, B32, 0.0, B43, B42, 0.0, B53, B52, 0.0, B63, B62, 0.0, B73, B72, 0.0, B83, B82,
		    B13, 0.0, B11, B23, 0.0, B21, B33, 0.0, B31, B43, 0.0, B41, B53, 0.0, B51, B63, 0.0, B61, B73, 0.0, B71, B83, 0.0, B81;
	return Bij;
}

//=======================================================================
//Evaluates the stretching parameters of PML
//=======================================================================
Eigen::VectorXd
PML3DGeneral::ComputePMLStretchingFactors(const double ri, const double si, const double ti, const double rho, const double mu, const double lambda) const {
	//Gets the element coordinates in deformed configuration. 
    Vector X1 = nodePointers[0]->getCrds();
    Vector X2 = nodePointers[1]->getCrds();
    Vector X3 = nodePointers[2]->getCrds();
    Vector X4 = nodePointers[3]->getCrds();
    Vector X5 = nodePointers[4]->getCrds();
    Vector X6 = nodePointers[5]->getCrds();
    Vector X7 = nodePointers[6]->getCrds();
    Vector X8 = nodePointers[7]->getCrds();

	Eigen::VectorXd X(24);
    X << X1(0), X1(1), X1(2), X2(0), X2(1), X2(2), X3(0), X3(1), X3(2), X4(0), X4(1), X4(2), X5(0), X5(1), X5(2), X6(0), X6(1), X6(2), X7(0), X7(1), X7(2), X8(0), X8(1), X8(2);

	Eigen::MatrixXd Hij = ComputeShapeFunctionMatrix(ri, si, ti);

	Eigen::VectorXd XGauss = Hij*X;

	double x = XGauss(0);
	double y = XGauss(1);
	double z = XGauss(2);

	//P wave velocity
	double V_pml = sqrt((lambda + 2.0*mu)/rho);
	
	//characteristic length
	double b_pml = L_pml/10.0;

	//stretching parameters
	double a0 = (m_pml+1.0)*b_pml/2.0/L_pml*log(1.0/R_pml);
	double b0 = (m_pml+1.0)*V_pml/2.0/L_pml*log(1.0/R_pml);

	double ax = 1.0+a0*pow((x-x0_pml)*nx_pml/L_pml,m_pml);
	double ay = 1.0+a0*pow((y-y0_pml)*ny_pml/L_pml,m_pml);
	double az = 1.0+a0*pow((z-z0_pml)*nz_pml/L_pml,m_pml);

	double bx = b0*pow((x-x0_pml)*nx_pml/L_pml,m_pml);
	double by = b0*pow((y-y0_pml)*ny_pml/L_pml,m_pml);
	double bz = b0*pow((z-z0_pml)*nz_pml/L_pml,m_pml);

	Eigen::VectorXd abc_pml(6);
	abc_pml << ax, ay, az, bx, by, bz;

	return abc_pml;
}

// =======================================================================
//	compute mass matrix
// =======================================================================
void PML3DGeneral::ComputeMassMatrix()
{
    MassMatrix.Zero();

    Eigen::VectorXd wi(8);
    Eigen::MatrixXd xi(8,3);

    wi << 1.0000,  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000;
    xi << -0.577350269189626, -0.577350269189626, -0.577350269189626, 
           0.577350269189626, -0.577350269189626, -0.577350269189626, 
          -0.577350269189626,  0.577350269189626, -0.577350269189626, 
           0.577350269189626,  0.577350269189626, -0.577350269189626, 
          -0.577350269189626, -0.577350269189626,  0.577350269189626, 
           0.577350269189626, -0.577350269189626,  0.577350269189626, 
          -0.577350269189626,  0.577350269189626,  0.577350269189626, 
           0.577350269189626,  0.577350269189626,  0.577350269189626;


    double lambda = E*nu/(1.0 + nu)/(1.0 - 2.0*nu);
    double mu = E/2.0/(1.0 + nu);

	//Numerical integration.
    for(unsigned int i = 0; i < wi.size(); i++){
		double ri = xi(i,0);
		double si = xi(i,1);
		double ti = xi(i,2);

		Eigen::VectorXd abc_pml = ComputePMLStretchingFactors(ri, si, ti, rho, mu, lambda);

		double a = abc_pml(0)*abc_pml(1)*abc_pml(2);

		double H11 = 1.0/8.0*(1.0 - ri)*(1.0 - si)*(1.0 - ti);
		double H22 = 1.0/8.0*(1.0 + ri)*(1.0 - si)*(1.0 - ti);
		double H33 = 1.0/8.0*(1.0 + ri)*(1.0 + si)*(1.0 - ti);
		double H44 = 1.0/8.0*(1.0 - ri)*(1.0 + si)*(1.0 - ti);
		double H55 = 1.0/8.0*(1.0 - ri)*(1.0 - si)*(1.0 + ti);
		double H66 = 1.0/8.0*(1.0 + ri)*(1.0 - si)*(1.0 + ti);
		double H77 = 1.0/8.0*(1.0 + ri)*(1.0 + si)*(1.0 + ti);
		double H88 = 1.0/8.0*(1.0 - ri)*(1.0 + si)*(1.0 + ti);

		Eigen::VectorXd SF(8);
		SF << H11, H22, H33, H44, H55, H66, H77, H88;

		//Jacobian matrix:
		Eigen::MatrixXd Jij = ComputeJacobianMatrix(ri, si, ti);
		double D = fabs(Jij.determinant());

		for (unsigned int j = 0; j < 8 ; j++) {
			for (unsigned int k = 0; k < 8 ; k++) {
				MassMatrix(9*j  ,9*k  ) += rho*a*SF(j)*SF(k)*wi(i)*D;
				MassMatrix(9*j+1,9*k+1) += rho*a*SF(j)*SF(k)*wi(i)*D;
				MassMatrix(9*j+2,9*k+2) += rho*a*SF(j)*SF(k)*wi(i)*D;
				MassMatrix(9*j+3,9*k+3) += -a*(lambda+mu)/mu/(3*lambda+2*mu)*SF(j)*SF(k)*wi(i)*D;
				MassMatrix(9*j+4,9*k+4) += -a*(lambda+mu)/mu/(3*lambda+2*mu)*SF(j)*SF(k)*wi(i)*D;
				MassMatrix(9*j+5,9*k+5) += -a*(lambda+mu)/mu/(3*lambda+2*mu)*SF(j)*SF(k)*wi(i)*D;
				MassMatrix(9*j+6,9*k+6) += -a/mu*SF(j)*SF(k)*wi(i)*D;
				MassMatrix(9*j+7,9*k+7) += -a/mu*SF(j)*SF(k)*wi(i)*D;
				MassMatrix(9*j+8,9*k+8) += -a/mu*SF(j)*SF(k)*wi(i)*D;
				
				MassMatrix(9*j+3,9*k+4) += a*lambda/2/mu/(3*lambda+2*mu)*SF(j)*SF(k)*wi(i)*D;
				MassMatrix(9*j+3,9*k+5) += a*lambda/2/mu/(3*lambda+2*mu)*SF(j)*SF(k)*wi(i)*D;
				MassMatrix(9*j+4,9*k+5) += a*lambda/2/mu/(3*lambda+2*mu)*SF(j)*SF(k)*wi(i)*D;

				MassMatrix(9*j+4,9*k+3) += a*lambda/2/mu/(3*lambda+2*mu)*SF(j)*SF(k)*wi(i)*D;
				MassMatrix(9*j+5,9*k+3) += a*lambda/2/mu/(3*lambda+2*mu)*SF(j)*SF(k)*wi(i)*D;
				MassMatrix(9*j+5,9*k+4) += a*lambda/2/mu/(3*lambda+2*mu)*SF(j)*SF(k)*wi(i)*D;
			}
		}
	}
}

// =======================================================================
//	compute damping matrix
// =======================================================================
void PML3DGeneral::ComputeDampingMatrix()
{
    DampingMatrix.Zero();

    Eigen::VectorXd wi(8);
    Eigen::MatrixXd xi(8,3);

    wi << 1.0000,  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000;
    xi << -0.577350269189626, -0.577350269189626, -0.577350269189626, 
           0.577350269189626, -0.577350269189626, -0.577350269189626, 
          -0.577350269189626,  0.577350269189626, -0.577350269189626, 
           0.577350269189626,  0.577350269189626, -0.577350269189626, 
          -0.577350269189626, -0.577350269189626,  0.577350269189626, 
           0.577350269189626, -0.577350269189626,  0.577350269189626, 
          -0.577350269189626,  0.577350269189626,  0.577350269189626, 
           0.577350269189626,  0.577350269189626,  0.577350269189626;


    double lambda = E*nu/(1.0 + nu)/(1.0 - 2.0*nu);
    double mu = E/2.0/(1.0 + nu);


	//Numerical integration.
    for(unsigned int i = 0; i < wi.size(); i++){
		double ri = xi(i,0);
		double si = xi(i,1);
		double ti = xi(i,2);

		Eigen::VectorXd abc_pml = ComputePMLStretchingFactors(ri, si, ti, rho, mu, lambda);
		double ax = abc_pml(0);
		double ay = abc_pml(1);
		double az = abc_pml(2);
		double bx = abc_pml(3);
		double by = abc_pml(4);
		double bz = abc_pml(5);

		double b  = ax*ay*bz + ax*by*az + bx*ay*az;

		double H11 = 1.0/8.0*(1.0 - ri)*(1.0 - si)*(1.0 - ti);
		double H22 = 1.0/8.0*(1.0 + ri)*(1.0 - si)*(1.0 - ti);
		double H33 = 1.0/8.0*(1.0 + ri)*(1.0 + si)*(1.0 - ti);
		double H44 = 1.0/8.0*(1.0 - ri)*(1.0 + si)*(1.0 - ti);
		double H55 = 1.0/8.0*(1.0 - ri)*(1.0 - si)*(1.0 + ti);
		double H66 = 1.0/8.0*(1.0 + ri)*(1.0 - si)*(1.0 + ti);
		double H77 = 1.0/8.0*(1.0 + ri)*(1.0 + si)*(1.0 + ti);
		double H88 = 1.0/8.0*(1.0 - ri)*(1.0 + si)*(1.0 + ti);

		//Jacobian matrix:
		Eigen::MatrixXd Jij = ComputeJacobianMatrix(ri, si, ti);
		
		double D = fabs(Jij.determinant());

		Eigen::MatrixXd J = Jij.inverse();

		//Strain-displacement matrix coefficients:
		double B11 = -1.0/8.0*J(0,2)*(1.0 - ri)*(1.0 - si) - 1.0/8.0*J(0,1)*(1.0 - ri)*(1.0 - ti) - 1.0/8.0*J(0,0)*(1.0 - si)*(1.0 - ti);
		double B21 = -1.0/8.0*J(0,2)*(1.0 + ri)*(1.0 - si) - 1.0/8.0*J(0,1)*(1.0 + ri)*(1.0 - ti) + 1.0/8.0*J(0,0)*(1.0 - si)*(1.0 - ti);
		double B31 = -1.0/8.0*J(0,2)*(1.0 + ri)*(1.0 + si) + 1.0/8.0*J(0,1)*(1.0 + ri)*(1.0 - ti) + 1.0/8.0*J(0,0)*(1.0 + si)*(1.0 - ti);
		double B41 = -1.0/8.0*J(0,2)*(1.0 - ri)*(1.0 + si) + 1.0/8.0*J(0,1)*(1.0 - ri)*(1.0 - ti) - 1.0/8.0*J(0,0)*(1.0 + si)*(1.0 - ti);
		double B51 =  1.0/8.0*J(0,2)*(1.0 - ri)*(1.0 - si) - 1.0/8.0*J(0,1)*(1.0 - ri)*(1.0 + ti) - 1.0/8.0*J(0,0)*(1.0 - si)*(1.0 + ti);
		double B61 =  1.0/8.0*J(0,2)*(1.0 + ri)*(1.0 - si) - 1.0/8.0*J(0,1)*(1.0 + ri)*(1.0 + ti) + 1.0/8.0*J(0,0)*(1.0 - si)*(1.0 + ti);
		double B71 =  1.0/8.0*J(0,2)*(1.0 + ri)*(1.0 + si) + 1.0/8.0*J(0,1)*(1.0 + ri)*(1.0 + ti) + 1.0/8.0*J(0,0)*(1.0 + si)*(1.0 + ti);
		double B81 =  1.0/8.0*J(0,2)*(1.0 - ri)*(1.0 + si) + 1.0/8.0*J(0,1)*(1.0 - ri)*(1.0 + ti) - 1.0/8.0*J(0,0)*(1.0 + si)*(1.0 + ti);

		double B12 = -1.0/8.0*J(1,2)*(1.0 - ri)*(1.0 - si) - 1.0/8.0*J(1,1)*(1.0 - ri)*(1.0 - ti) - 1.0/8.0*J(1,0)*(1.0 - si)*(1.0 - ti);
		double B22 = -1.0/8.0*J(1,2)*(1.0 + ri)*(1.0 - si) - 1.0/8.0*J(1,1)*(1.0 + ri)*(1.0 - ti) + 1.0/8.0*J(1,0)*(1.0 - si)*(1.0 - ti);
		double B32 = -1.0/8.0*J(1,2)*(1.0 + ri)*(1.0 + si) + 1.0/8.0*J(1,1)*(1.0 + ri)*(1.0 - ti) + 1.0/8.0*J(1,0)*(1.0 + si)*(1.0 - ti);
		double B42 = -1.0/8.0*J(1,2)*(1.0 - ri)*(1.0 + si) + 1.0/8.0*J(1,1)*(1.0 - ri)*(1.0 - ti) - 1.0/8.0*J(1,0)*(1.0 + si)*(1.0 - ti);
		double B52 =  1.0/8.0*J(1,2)*(1.0 - ri)*(1.0 - si) - 1.0/8.0*J(1,1)*(1.0 - ri)*(1.0 + ti) - 1.0/8.0*J(1,0)*(1.0 - si)*(1.0 + ti);
		double B62 =  1.0/8.0*J(1,2)*(1.0 + ri)*(1.0 - si) - 1.0/8.0*J(1,1)*(1.0 + ri)*(1.0 + ti) + 1.0/8.0*J(1,0)*(1.0 - si)*(1.0 + ti);
		double B72 =  1.0/8.0*J(1,2)*(1.0 + ri)*(1.0 + si) + 1.0/8.0*J(1,1)*(1.0 + ri)*(1.0 + ti) + 1.0/8.0*J(1,0)*(1.0 + si)*(1.0 + ti);
		double B82 =  1.0/8.0*J(1,2)*(1.0 - ri)*(1.0 + si) + 1.0/8.0*J(1,1)*(1.0 - ri)*(1.0 + ti) - 1.0/8.0*J(1,0)*(1.0 + si)*(1.0 + ti);

		double B13 = -1.0/8.0*J(2,2)*(1.0 - ri)*(1.0 - si) - 1.0/8.0*J(2,1)*(1.0 - ri)*(1.0 - ti) - 1.0/8.0*J(2,0)*(1.0 - si)*(1.0 - ti);
		double B23 = -1.0/8.0*J(2,2)*(1.0 + ri)*(1.0 - si) - 1.0/8.0*J(2,1)*(1.0 + ri)*(1.0 - ti) + 1.0/8.0*J(2,0)*(1.0 - si)*(1.0 - ti);
		double B33 = -1.0/8.0*J(2,2)*(1.0 + ri)*(1.0 + si) + 1.0/8.0*J(2,1)*(1.0 + ri)*(1.0 - ti) + 1.0/8.0*J(2,0)*(1.0 + si)*(1.0 - ti);
		double B43 = -1.0/8.0*J(2,2)*(1.0 - ri)*(1.0 + si) + 1.0/8.0*J(2,1)*(1.0 - ri)*(1.0 - ti) - 1.0/8.0*J(2,0)*(1.0 + si)*(1.0 - ti);
		double B53 =  1.0/8.0*J(2,2)*(1.0 - ri)*(1.0 - si) - 1.0/8.0*J(2,1)*(1.0 - ri)*(1.0 + ti) - 1.0/8.0*J(2,0)*(1.0 - si)*(1.0 + ti);
		double B63 =  1.0/8.0*J(2,2)*(1.0 + ri)*(1.0 - si) - 1.0/8.0*J(2,1)*(1.0 + ri)*(1.0 + ti) + 1.0/8.0*J(2,0)*(1.0 - si)*(1.0 + ti);
		double B73 =  1.0/8.0*J(2,2)*(1.0 + ri)*(1.0 + si) + 1.0/8.0*J(2,1)*(1.0 + ri)*(1.0 + ti) + 1.0/8.0*J(2,0)*(1.0 + si)*(1.0 + ti);
		double B83 =  1.0/8.0*J(2,2)*(1.0 - ri)*(1.0 + si) + 1.0/8.0*J(2,1)*(1.0 - ri)*(1.0 + ti) - 1.0/8.0*J(2,0)*(1.0 + si)*(1.0 + ti);

		Eigen::VectorXd SF(8);
		SF << H11, H22, H33, H44, H55, H66, H77, H88;

		Eigen::VectorXd dSFdx(8);
		dSFdx << B11, B21, B31, B41, B51, B61, B71, B81;

		Eigen::VectorXd dSFdy(8);
		dSFdy << B12, B22, B32, B42, B52, B62, B72, B82;

		Eigen::VectorXd dSFdz(8);
		dSFdz << B13, B23, B33, B43, B53, B63, B73, B83;

		for (unsigned int j = 0; j < 8 ; j++) {
			for (unsigned int k = 0; k < 8 ; k++) {
				//The state vector: U = [u1, u2, u3, s11, s22, s33, s12, s23, s13]
				DampingMatrix(9*j  ,9*k  ) += rho*b*SF(j)*SF(k)*wi(i)*D;
				DampingMatrix(9*j+1,9*k+1) += rho*b*SF(j)*SF(k)*wi(i)*D;
				DampingMatrix(9*j+2,9*k+2) += rho*b*SF(j)*SF(k)*wi(i)*D;
				DampingMatrix(9*j+3,9*k+3) += -b*(lambda + mu)/mu/(3.0*lambda + 2.0*mu)*SF(j)*SF(k)*wi(i)*D;
				DampingMatrix(9*j+4,9*k+4) += -b*(lambda + mu)/mu/(3.0*lambda + 2.0*mu)*SF(j)*SF(k)*wi(i)*D;
				DampingMatrix(9*j+5,9*k+5) += -b*(lambda + mu)/mu/(3.0*lambda + 2.0*mu)*SF(j)*SF(k)*wi(i)*D;
				DampingMatrix(9*j+6,9*k+6) += -b/mu*SF(j)*SF(k)*wi(i)*D;
				DampingMatrix(9*j+7,9*k+7) += -b/mu*SF(j)*SF(k)*wi(i)*D;
				DampingMatrix(9*j+8,9*k+8) += -b/mu*SF(j)*SF(k)*wi(i)*D;
				
				DampingMatrix(9*j+3,9*k+4) += b*lambda/2.0/mu/(3.0*lambda + 2.0*mu)*SF(j)*SF(k)*wi(i)*D;
				DampingMatrix(9*j+3,9*k+5) += b*lambda/2.0/mu/(3.0*lambda + 2.0*mu)*SF(j)*SF(k)*wi(i)*D;
				DampingMatrix(9*j+4,9*k+5) += b*lambda/2.0/mu/(3.0*lambda + 2.0*mu)*SF(j)*SF(k)*wi(i)*D;

				DampingMatrix(9*j+4,9*k+3) += b*lambda/2.0/mu/(3.0*lambda + 2.0*mu)*SF(j)*SF(k)*wi(i)*D;
				DampingMatrix(9*j+5,9*k+3) += b*lambda/2.0/mu/(3.0*lambda + 2.0*mu)*SF(j)*SF(k)*wi(i)*D;
				DampingMatrix(9*j+5,9*k+4) += b*lambda/2.0/mu/(3.0*lambda + 2.0*mu)*SF(j)*SF(k)*wi(i)*D;

				DampingMatrix(9*j  ,9*k+3) += dSFdx(j)*SF(k)*(ay*az)*wi(i)*D; 
				DampingMatrix(9*j  ,9*k+6) += dSFdy(j)*SF(k)*(ax*az)*wi(i)*D;
				DampingMatrix(9*j  ,9*k+8) += dSFdz(j)*SF(k)*(ax*ay)*wi(i)*D;
				DampingMatrix(9*j+1,9*k+4) += dSFdy(j)*SF(k)*(ax*az)*wi(i)*D;
				DampingMatrix(9*j+1,9*k+6) += dSFdx(j)*SF(k)*(ay*az)*wi(i)*D;
				DampingMatrix(9*j+1,9*k+7) += dSFdz(j)*SF(k)*(ax*ay)*wi(i)*D;
				DampingMatrix(9*j+2,9*k+5) += dSFdz(j)*SF(k)*(ax*ay)*wi(i)*D;
				DampingMatrix(9*j+2,9*k+8) += dSFdx(j)*SF(k)*(ay*az)*wi(i)*D;
				DampingMatrix(9*j+2,9*k+7) += dSFdy(j)*SF(k)*(ax*az)*wi(i)*D;

				DampingMatrix(9*j+3,9*k  ) += dSFdx(k)*SF(j)*(ay*az)*wi(i)*D; 
				DampingMatrix(9*j+6,9*k  ) += dSFdy(k)*SF(j)*(ax*az)*wi(i)*D;
				DampingMatrix(9*j+8,9*k  ) += dSFdz(k)*SF(j)*(ax*ay)*wi(i)*D;
				DampingMatrix(9*j+4,9*k+1) += dSFdy(k)*SF(j)*(ax*az)*wi(i)*D;
				DampingMatrix(9*j+6,9*k+1) += dSFdx(k)*SF(j)*(ay*az)*wi(i)*D;
				DampingMatrix(9*j+7,9*k+1) += dSFdz(k)*SF(j)*(ax*ay)*wi(i)*D;
				DampingMatrix(9*j+5,9*k+2) += dSFdz(k)*SF(j)*(ax*ay)*wi(i)*D;
				DampingMatrix(9*j+8,9*k+2) += dSFdx(k)*SF(j)*(ay*az)*wi(i)*D;
				DampingMatrix(9*j+7,9*k+2) += dSFdy(k)*SF(j)*(ax*az)*wi(i)*D;
			}
		}
	}
	
}

// =======================================================================
//	calculating the PML matrix
// =======================================================================
void PML3DGeneral::ComputePMLMatrix() {
	//Impedance buffer matrix definition:
    ImpedanceMatrix.Zero();

    Eigen::VectorXd wi(8);
    Eigen::MatrixXd xi(8,3);

    wi << 1.0000,  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000;
    xi << -0.577350269189626, -0.577350269189626, -0.577350269189626, 
           0.577350269189626, -0.577350269189626, -0.577350269189626, 
          -0.577350269189626,  0.577350269189626, -0.577350269189626, 
           0.577350269189626,  0.577350269189626, -0.577350269189626, 
          -0.577350269189626, -0.577350269189626,  0.577350269189626, 
           0.577350269189626, -0.577350269189626,  0.577350269189626, 
          -0.577350269189626,  0.577350269189626,  0.577350269189626, 
           0.577350269189626,  0.577350269189626,  0.577350269189626;


    double lambda = E*nu/(1.0 + nu)/(1.0 - 2.0*nu);
    double mu = E/2.0/(1.0 + nu);

	//Numerical integration.
    for(unsigned int i = 0; i < wi.size(); i++){
        //Gets the integration point coordinates.
		double ri = xi(i,0);
		double si = xi(i,1);
		double ti = xi(i,2);

		Eigen::VectorXd abc_pml = ComputePMLStretchingFactors(ri, si, ti, rho, mu, lambda);

		//double ax = abc_pml(0);
		//double ay = abc_pml(1);
		//double az = abc_pml(2);
		double bx = abc_pml(3);
		double by = abc_pml(4);
		double bz = abc_pml(5);

		double d  = bx*by*bz;

		double H11 = 1.0/8.0*(1.0 - ri)*(1.0 - si)*(1.0 - ti);
		double H22 = 1.0/8.0*(1.0 + ri)*(1.0 - si)*(1.0 - ti);
		double H33 = 1.0/8.0*(1.0 + ri)*(1.0 + si)*(1.0 - ti);
		double H44 = 1.0/8.0*(1.0 - ri)*(1.0 + si)*(1.0 - ti);
		double H55 = 1.0/8.0*(1.0 - ri)*(1.0 - si)*(1.0 + ti);
		double H66 = 1.0/8.0*(1.0 + ri)*(1.0 - si)*(1.0 + ti);
		double H77 = 1.0/8.0*(1.0 + ri)*(1.0 + si)*(1.0 + ti);
		double H88 = 1.0/8.0*(1.0 - ri)*(1.0 + si)*(1.0 + ti);

		//Jacobian matrix:
		Eigen::MatrixXd Jij = ComputeJacobianMatrix(ri, si, ti);
		
		double D = fabs(Jij.determinant());

		Eigen::MatrixXd J = Jij.inverse();

		//Strain-displacement matrix coefficients:
		double B11 = -1.0/8.0*J(0,2)*(1.0 - ri)*(1.0 - si) - 1.0/8.0*J(0,1)*(1.0 - ri)*(1.0 - ti) - 1.0/8.0*J(0,0)*(1.0 - si)*(1.0 - ti);
		double B21 = -1.0/8.0*J(0,2)*(1.0 + ri)*(1.0 - si) - 1.0/8.0*J(0,1)*(1.0 + ri)*(1.0 - ti) + 1.0/8.0*J(0,0)*(1.0 - si)*(1.0 - ti);
		double B31 = -1.0/8.0*J(0,2)*(1.0 + ri)*(1.0 + si) + 1.0/8.0*J(0,1)*(1.0 + ri)*(1.0 - ti) + 1.0/8.0*J(0,0)*(1.0 + si)*(1.0 - ti);
		double B41 = -1.0/8.0*J(0,2)*(1.0 - ri)*(1.0 + si) + 1.0/8.0*J(0,1)*(1.0 - ri)*(1.0 - ti) - 1.0/8.0*J(0,0)*(1.0 + si)*(1.0 - ti);
		double B51 =  1.0/8.0*J(0,2)*(1.0 - ri)*(1.0 - si) - 1.0/8.0*J(0,1)*(1.0 - ri)*(1.0 + ti) - 1.0/8.0*J(0,0)*(1.0 - si)*(1.0 + ti);
		double B61 =  1.0/8.0*J(0,2)*(1.0 + ri)*(1.0 - si) - 1.0/8.0*J(0,1)*(1.0 + ri)*(1.0 + ti) + 1.0/8.0*J(0,0)*(1.0 - si)*(1.0 + ti);
		double B71 =  1.0/8.0*J(0,2)*(1.0 + ri)*(1.0 + si) + 1.0/8.0*J(0,1)*(1.0 + ri)*(1.0 + ti) + 1.0/8.0*J(0,0)*(1.0 + si)*(1.0 + ti);
		double B81 =  1.0/8.0*J(0,2)*(1.0 - ri)*(1.0 + si) + 1.0/8.0*J(0,1)*(1.0 - ri)*(1.0 + ti) - 1.0/8.0*J(0,0)*(1.0 + si)*(1.0 + ti);

		double B12 = -1.0/8.0*J(1,2)*(1.0 - ri)*(1.0 - si) - 1.0/8.0*J(1,1)*(1.0 - ri)*(1.0 - ti) - 1.0/8.0*J(1,0)*(1.0 - si)*(1.0 - ti);
		double B22 = -1.0/8.0*J(1,2)*(1.0 + ri)*(1.0 - si) - 1.0/8.0*J(1,1)*(1.0 + ri)*(1.0 - ti) + 1.0/8.0*J(1,0)*(1.0 - si)*(1.0 - ti);
		double B32 = -1.0/8.0*J(1,2)*(1.0 + ri)*(1.0 + si) + 1.0/8.0*J(1,1)*(1.0 + ri)*(1.0 - ti) + 1.0/8.0*J(1,0)*(1.0 + si)*(1.0 - ti);
		double B42 = -1.0/8.0*J(1,2)*(1.0 - ri)*(1.0 + si) + 1.0/8.0*J(1,1)*(1.0 - ri)*(1.0 - ti) - 1.0/8.0*J(1,0)*(1.0 + si)*(1.0 - ti);
		double B52 =  1.0/8.0*J(1,2)*(1.0 - ri)*(1.0 - si) - 1.0/8.0*J(1,1)*(1.0 - ri)*(1.0 + ti) - 1.0/8.0*J(1,0)*(1.0 - si)*(1.0 + ti);
		double B62 =  1.0/8.0*J(1,2)*(1.0 + ri)*(1.0 - si) - 1.0/8.0*J(1,1)*(1.0 + ri)*(1.0 + ti) + 1.0/8.0*J(1,0)*(1.0 - si)*(1.0 + ti);
		double B72 =  1.0/8.0*J(1,2)*(1.0 + ri)*(1.0 + si) + 1.0/8.0*J(1,1)*(1.0 + ri)*(1.0 + ti) + 1.0/8.0*J(1,0)*(1.0 + si)*(1.0 + ti);
		double B82 =  1.0/8.0*J(1,2)*(1.0 - ri)*(1.0 + si) + 1.0/8.0*J(1,1)*(1.0 - ri)*(1.0 + ti) - 1.0/8.0*J(1,0)*(1.0 + si)*(1.0 + ti);

		double B13 = -1.0/8.0*J(2,2)*(1.0 - ri)*(1.0 - si) - 1.0/8.0*J(2,1)*(1.0 - ri)*(1.0 - ti) - 1.0/8.0*J(2,0)*(1.0 - si)*(1.0 - ti);
		double B23 = -1.0/8.0*J(2,2)*(1.0 + ri)*(1.0 - si) - 1.0/8.0*J(2,1)*(1.0 + ri)*(1.0 - ti) + 1.0/8.0*J(2,0)*(1.0 - si)*(1.0 - ti);
		double B33 = -1.0/8.0*J(2,2)*(1.0 + ri)*(1.0 + si) + 1.0/8.0*J(2,1)*(1.0 + ri)*(1.0 - ti) + 1.0/8.0*J(2,0)*(1.0 + si)*(1.0 - ti);
		double B43 = -1.0/8.0*J(2,2)*(1.0 - ri)*(1.0 + si) + 1.0/8.0*J(2,1)*(1.0 - ri)*(1.0 - ti) - 1.0/8.0*J(2,0)*(1.0 + si)*(1.0 - ti);
		double B53 =  1.0/8.0*J(2,2)*(1.0 - ri)*(1.0 - si) - 1.0/8.0*J(2,1)*(1.0 - ri)*(1.0 + ti) - 1.0/8.0*J(2,0)*(1.0 - si)*(1.0 + ti);
		double B63 =  1.0/8.0*J(2,2)*(1.0 + ri)*(1.0 - si) - 1.0/8.0*J(2,1)*(1.0 + ri)*(1.0 + ti) + 1.0/8.0*J(2,0)*(1.0 - si)*(1.0 + ti);
		double B73 =  1.0/8.0*J(2,2)*(1.0 + ri)*(1.0 + si) + 1.0/8.0*J(2,1)*(1.0 + ri)*(1.0 + ti) + 1.0/8.0*J(2,0)*(1.0 + si)*(1.0 + ti);
		double B83 =  1.0/8.0*J(2,2)*(1.0 - ri)*(1.0 + si) + 1.0/8.0*J(2,1)*(1.0 - ri)*(1.0 + ti) - 1.0/8.0*J(2,0)*(1.0 + si)*(1.0 + ti);

		Eigen::VectorXd SF(8);
		SF << H11, H22, H33, H44, H55, H66, H77, H88;

		Eigen::VectorXd dSFdx(8);
		dSFdx << B11, B21, B31, B41, B51, B61, B71, B81;

		Eigen::VectorXd dSFdy(8);
		dSFdy << B12, B22, B32, B42, B52, B62, B72, B82;

		Eigen::VectorXd dSFdz(8);
		dSFdz << B13, B23, B33, B43, B53, B63, B73, B83;

		for (unsigned int j = 0; j < 8 ; j++) {
			for (unsigned int k = 0; k < 8 ; k++) {
				//The state vector: U = [u1, u2, u3, s11, s22, s33, s12, s23, s13]
				ImpedanceMatrix(9*j  ,9*k  ) += rho*d*SF(j)*SF(k)*wi(i)*D;
				ImpedanceMatrix(9*j+1,9*k+1) += rho*d*SF(j)*SF(k)*wi(i)*D;
				ImpedanceMatrix(9*j+2,9*k+2) += rho*d*SF(j)*SF(k)*wi(i)*D;

				ImpedanceMatrix(9*j+3,9*k+3) += -d*(lambda + mu)/mu/(3.0*lambda + 2.0*mu)*SF(j)*SF(k)*wi(i)*D;
				ImpedanceMatrix(9*j+4,9*k+4) += -d*(lambda + mu)/mu/(3.0*lambda + 2.0*mu)*SF(j)*SF(k)*wi(i)*D;
				ImpedanceMatrix(9*j+5,9*k+5) += -d*(lambda + mu)/mu/(3.0*lambda + 2.0*mu)*SF(j)*SF(k)*wi(i)*D;
				ImpedanceMatrix(9*j+6,9*k+6) += -d/mu*SF(j)*SF(k)*wi(i)*D;
				ImpedanceMatrix(9*j+7,9*k+7) += -d/mu*SF(j)*SF(k)*wi(i)*D;
				ImpedanceMatrix(9*j+8,9*k+8) += -d/mu*SF(j)*SF(k)*wi(i)*D;
				
				ImpedanceMatrix(9*j+3,9*k+4) += d*lambda/2.0/mu/(3.0*lambda + 2.0*mu)*SF(j)*SF(k)*wi(i)*D;
				ImpedanceMatrix(9*j+3,9*k+5) += d*lambda/2.0/mu/(3.0*lambda + 2.0*mu)*SF(j)*SF(k)*wi(i)*D;
				ImpedanceMatrix(9*j+4,9*k+5) += d*lambda/2.0/mu/(3.0*lambda + 2.0*mu)*SF(j)*SF(k)*wi(i)*D;

				ImpedanceMatrix(9*j+4,9*k+3) += d*lambda/2.0/mu/(3.0*lambda + 2.0*mu)*SF(j)*SF(k)*wi(i)*D;
				ImpedanceMatrix(9*j+5,9*k+3) += d*lambda/2.0/mu/(3.0*lambda + 2.0*mu)*SF(j)*SF(k)*wi(i)*D;
				ImpedanceMatrix(9*j+5,9*k+4) += d*lambda/2.0/mu/(3.0*lambda + 2.0*mu)*SF(j)*SF(k)*wi(i)*D;

				ImpedanceMatrix(9*j  ,9*k+3) += dSFdx(j)*SF(k)*(by*bz)*wi(i)*D; 
				ImpedanceMatrix(9*j  ,9*k+6) += dSFdy(j)*SF(k)*(bx*bz)*wi(i)*D;
				ImpedanceMatrix(9*j  ,9*k+8) += dSFdz(j)*SF(k)*(bx*by)*wi(i)*D;
				ImpedanceMatrix(9*j+1,9*k+4) += dSFdy(j)*SF(k)*(bx*bz)*wi(i)*D;
				ImpedanceMatrix(9*j+1,9*k+6) += dSFdx(j)*SF(k)*(by*bz)*wi(i)*D;
				ImpedanceMatrix(9*j+1,9*k+7) += dSFdz(j)*SF(k)*(bx*by)*wi(i)*D;
				ImpedanceMatrix(9*j+2,9*k+5) += dSFdz(j)*SF(k)*(bx*by)*wi(i)*D;
				ImpedanceMatrix(9*j+2,9*k+8) += dSFdx(j)*SF(k)*(by*bz)*wi(i)*D;
				ImpedanceMatrix(9*j+2,9*k+7) += dSFdy(j)*SF(k)*(bx*bz)*wi(i)*D;

				ImpedanceMatrix(9*j+3,9*k  ) += dSFdx(k)*SF(j)*(by*bz)*wi(i)*D; 
				ImpedanceMatrix(9*j+6,9*k  ) += dSFdy(k)*SF(j)*(bx*bz)*wi(i)*D;
				ImpedanceMatrix(9*j+8,9*k  ) += dSFdz(k)*SF(j)*(bx*by)*wi(i)*D;
				ImpedanceMatrix(9*j+4,9*k+1) += dSFdy(k)*SF(j)*(bx*bz)*wi(i)*D;
				ImpedanceMatrix(9*j+6,9*k+1) += dSFdx(k)*SF(j)*(by*bz)*wi(i)*D;
				ImpedanceMatrix(9*j+7,9*k+1) += dSFdz(k)*SF(j)*(bx*by)*wi(i)*D;
				ImpedanceMatrix(9*j+5,9*k+2) += dSFdz(k)*SF(j)*(bx*by)*wi(i)*D;
				ImpedanceMatrix(9*j+8,9*k+2) += dSFdx(k)*SF(j)*(by*bz)*wi(i)*D;
				ImpedanceMatrix(9*j+7,9*k+2) += dSFdy(k)*SF(j)*(bx*bz)*wi(i)*D;
			}
		}
	}
}

//======================================================================
// compute stiffness matrix
//======================================================================
void PML3DGeneral::ComputeStiffnessMatrix() {

    //Stiffness matrix definition:
    StiffnessMatrix.Zero();

    Eigen::VectorXd wi(8);
    Eigen::MatrixXd xi(8,3);

    wi << 1.0000,  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000;
    xi << -0.577350269189626, -0.577350269189626, -0.577350269189626, 
           0.577350269189626, -0.577350269189626, -0.577350269189626, 
          -0.577350269189626,  0.577350269189626, -0.577350269189626, 
           0.577350269189626,  0.577350269189626, -0.577350269189626, 
          -0.577350269189626, -0.577350269189626,  0.577350269189626, 
           0.577350269189626, -0.577350269189626,  0.577350269189626, 
          -0.577350269189626,  0.577350269189626,  0.577350269189626, 
           0.577350269189626,  0.577350269189626,  0.577350269189626;


    double lambda = E*nu/(1.0 + nu)/(1.0 - 2.0*nu);
    double mu = E/2.0/(1.0 + nu);


	//Numerical integration.
    for(unsigned int i = 0; i < wi.size(); i++){
		//Gets material properties:
		double ri = xi(i,0);
		double si = xi(i,1);
		double ti = xi(i,2);

		Eigen::VectorXd abc_pml = ComputePMLStretchingFactors(ri, si, ti, rho, mu, lambda);

		double ax = abc_pml(0);
		double ay = abc_pml(1);
		double az = abc_pml(2);
		double bx = abc_pml(3);
		double by = abc_pml(4);
		double bz = abc_pml(5);

		double c  = ax*by*bz + bx*by*az + bx*ay*bz;

		double H11 = 1.0/8.0*(1.0 - ri)*(1.0 - si)*(1.0 - ti);
		double H22 = 1.0/8.0*(1.0 + ri)*(1.0 - si)*(1.0 - ti);
		double H33 = 1.0/8.0*(1.0 + ri)*(1.0 + si)*(1.0 - ti);
		double H44 = 1.0/8.0*(1.0 - ri)*(1.0 + si)*(1.0 - ti);
		double H55 = 1.0/8.0*(1.0 - ri)*(1.0 - si)*(1.0 + ti);
		double H66 = 1.0/8.0*(1.0 + ri)*(1.0 - si)*(1.0 + ti);
		double H77 = 1.0/8.0*(1.0 + ri)*(1.0 + si)*(1.0 + ti);
		double H88 = 1.0/8.0*(1.0 - ri)*(1.0 + si)*(1.0 + ti);

		//Jacobian matrix:
		Eigen::MatrixXd Jij = ComputeJacobianMatrix(ri, si, ti);
		
		double D = fabs(Jij.determinant());

		Eigen::MatrixXd J = Jij.inverse();

		//Strain-displacement matrix coefficients:
		double B11 = -1.0/8.0*J(0,2)*(1.0 - ri)*(1.0 - si) - 1.0/8.0*J(0,1)*(1.0 - ri)*(1.0 - ti) - 1.0/8.0*J(0,0)*(1.0 - si)*(1.0 - ti);
		double B21 = -1.0/8.0*J(0,2)*(1.0 + ri)*(1.0 - si) - 1.0/8.0*J(0,1)*(1.0 + ri)*(1.0 - ti) + 1.0/8.0*J(0,0)*(1.0 - si)*(1.0 - ti);
		double B31 = -1.0/8.0*J(0,2)*(1.0 + ri)*(1.0 + si) + 1.0/8.0*J(0,1)*(1.0 + ri)*(1.0 - ti) + 1.0/8.0*J(0,0)*(1.0 + si)*(1.0 - ti);
		double B41 = -1.0/8.0*J(0,2)*(1.0 - ri)*(1.0 + si) + 1.0/8.0*J(0,1)*(1.0 - ri)*(1.0 - ti) - 1.0/8.0*J(0,0)*(1.0 + si)*(1.0 - ti);
		double B51 =  1.0/8.0*J(0,2)*(1.0 - ri)*(1.0 - si) - 1.0/8.0*J(0,1)*(1.0 - ri)*(1.0 + ti) - 1.0/8.0*J(0,0)*(1.0 - si)*(1.0 + ti);
		double B61 =  1.0/8.0*J(0,2)*(1.0 + ri)*(1.0 - si) - 1.0/8.0*J(0,1)*(1.0 + ri)*(1.0 + ti) + 1.0/8.0*J(0,0)*(1.0 - si)*(1.0 + ti);
		double B71 =  1.0/8.0*J(0,2)*(1.0 + ri)*(1.0 + si) + 1.0/8.0*J(0,1)*(1.0 + ri)*(1.0 + ti) + 1.0/8.0*J(0,0)*(1.0 + si)*(1.0 + ti);
		double B81 =  1.0/8.0*J(0,2)*(1.0 - ri)*(1.0 + si) + 1.0/8.0*J(0,1)*(1.0 - ri)*(1.0 + ti) - 1.0/8.0*J(0,0)*(1.0 + si)*(1.0 + ti);

		double B12 = -1.0/8.0*J(1,2)*(1.0 - ri)*(1.0 - si) - 1.0/8.0*J(1,1)*(1.0 - ri)*(1.0 - ti) - 1.0/8.0*J(1,0)*(1.0 - si)*(1.0 - ti);
		double B22 = -1.0/8.0*J(1,2)*(1.0 + ri)*(1.0 - si) - 1.0/8.0*J(1,1)*(1.0 + ri)*(1.0 - ti) + 1.0/8.0*J(1,0)*(1.0 - si)*(1.0 - ti);
		double B32 = -1.0/8.0*J(1,2)*(1.0 + ri)*(1.0 + si) + 1.0/8.0*J(1,1)*(1.0 + ri)*(1.0 - ti) + 1.0/8.0*J(1,0)*(1.0 + si)*(1.0 - ti);
		double B42 = -1.0/8.0*J(1,2)*(1.0 - ri)*(1.0 + si) + 1.0/8.0*J(1,1)*(1.0 - ri)*(1.0 - ti) - 1.0/8.0*J(1,0)*(1.0 + si)*(1.0 - ti);
		double B52 =  1.0/8.0*J(1,2)*(1.0 - ri)*(1.0 - si) - 1.0/8.0*J(1,1)*(1.0 - ri)*(1.0 + ti) - 1.0/8.0*J(1,0)*(1.0 - si)*(1.0 + ti);
		double B62 =  1.0/8.0*J(1,2)*(1.0 + ri)*(1.0 - si) - 1.0/8.0*J(1,1)*(1.0 + ri)*(1.0 + ti) + 1.0/8.0*J(1,0)*(1.0 - si)*(1.0 + ti);
		double B72 =  1.0/8.0*J(1,2)*(1.0 + ri)*(1.0 + si) + 1.0/8.0*J(1,1)*(1.0 + ri)*(1.0 + ti) + 1.0/8.0*J(1,0)*(1.0 + si)*(1.0 + ti);
		double B82 =  1.0/8.0*J(1,2)*(1.0 - ri)*(1.0 + si) + 1.0/8.0*J(1,1)*(1.0 - ri)*(1.0 + ti) - 1.0/8.0*J(1,0)*(1.0 + si)*(1.0 + ti);

		double B13 = -1.0/8.0*J(2,2)*(1.0 - ri)*(1.0 - si) - 1.0/8.0*J(2,1)*(1.0 - ri)*(1.0 - ti) - 1.0/8.0*J(2,0)*(1.0 - si)*(1.0 - ti);
		double B23 = -1.0/8.0*J(2,2)*(1.0 + ri)*(1.0 - si) - 1.0/8.0*J(2,1)*(1.0 + ri)*(1.0 - ti) + 1.0/8.0*J(2,0)*(1.0 - si)*(1.0 - ti);
		double B33 = -1.0/8.0*J(2,2)*(1.0 + ri)*(1.0 + si) + 1.0/8.0*J(2,1)*(1.0 + ri)*(1.0 - ti) + 1.0/8.0*J(2,0)*(1.0 + si)*(1.0 - ti);
		double B43 = -1.0/8.0*J(2,2)*(1.0 - ri)*(1.0 + si) + 1.0/8.0*J(2,1)*(1.0 - ri)*(1.0 - ti) - 1.0/8.0*J(2,0)*(1.0 + si)*(1.0 - ti);
		double B53 =  1.0/8.0*J(2,2)*(1.0 - ri)*(1.0 - si) - 1.0/8.0*J(2,1)*(1.0 - ri)*(1.0 + ti) - 1.0/8.0*J(2,0)*(1.0 - si)*(1.0 + ti);
		double B63 =  1.0/8.0*J(2,2)*(1.0 + ri)*(1.0 - si) - 1.0/8.0*J(2,1)*(1.0 + ri)*(1.0 + ti) + 1.0/8.0*J(2,0)*(1.0 - si)*(1.0 + ti);
		double B73 =  1.0/8.0*J(2,2)*(1.0 + ri)*(1.0 + si) + 1.0/8.0*J(2,1)*(1.0 + ri)*(1.0 + ti) + 1.0/8.0*J(2,0)*(1.0 + si)*(1.0 + ti);
		double B83 =  1.0/8.0*J(2,2)*(1.0 - ri)*(1.0 + si) + 1.0/8.0*J(2,1)*(1.0 - ri)*(1.0 + ti) - 1.0/8.0*J(2,0)*(1.0 + si)*(1.0 + ti);

		Eigen::VectorXd SF(8);
		SF << H11, H22, H33, H44, H55, H66, H77, H88;

		Eigen::VectorXd dSFdx(8);
		dSFdx << B11, B21, B31, B41, B51, B61, B71, B81;

		Eigen::VectorXd dSFdy(8);
		dSFdy << B12, B22, B32, B42, B52, B62, B72, B82;

		Eigen::VectorXd dSFdz(8);
		dSFdz << B13, B23, B33, B43, B53, B63, B73, B83;

		for (unsigned int j = 0; j < 8 ; j++) {
			for (unsigned int k = 0; k < 8 ; k++) {
				//The state vector: U = [u1, u2, u3, s11, s22, s33, s12, s23, s13]
				StiffnessMatrix(9*j  ,9*k  ) += rho*c*SF(j)*SF(k)*wi(i)*D;
				StiffnessMatrix(9*j+1,9*k+1) += rho*c*SF(j)*SF(k)*wi(i)*D;
				StiffnessMatrix(9*j+2,9*k+2) += rho*c*SF(j)*SF(k)*wi(i)*D;
				StiffnessMatrix(9*j+3,9*k+3) += -c*(lambda + mu)/mu/(3.0*lambda + 2.0*mu)*SF(j)*SF(k)*wi(i)*D;
				StiffnessMatrix(9*j+4,9*k+4) += -c*(lambda + mu)/mu/(3.0*lambda + 2.0*mu)*SF(j)*SF(k)*wi(i)*D;
				StiffnessMatrix(9*j+5,9*k+5) += -c*(lambda + mu)/mu/(3.0*lambda + 2.0*mu)*SF(j)*SF(k)*wi(i)*D;
				StiffnessMatrix(9*j+6,9*k+6) += -c/mu*SF(j)*SF(k)*wi(i)*D;
				StiffnessMatrix(9*j+7,9*k+7) += -c/mu*SF(j)*SF(k)*wi(i)*D;
				StiffnessMatrix(9*j+8,9*k+8) += -c/mu*SF(j)*SF(k)*wi(i)*D;
				
				StiffnessMatrix(9*j+3,9*k+4) += c*lambda/2.0/mu/(3.0*lambda + 2.0*mu)*SF(j)*SF(k)*wi(i)*D;
				StiffnessMatrix(9*j+3,9*k+5) += c*lambda/2.0/mu/(3.0*lambda + 2.0*mu)*SF(j)*SF(k)*wi(i)*D;
				StiffnessMatrix(9*j+4,9*k+5) += c*lambda/2.0/mu/(3.0*lambda + 2.0*mu)*SF(j)*SF(k)*wi(i)*D;

				StiffnessMatrix(9*j+4,9*k+3) += c*lambda/2.0/mu/(3.0*lambda + 2.0*mu)*SF(j)*SF(k)*wi(i)*D;
				StiffnessMatrix(9*j+5,9*k+3) += c*lambda/2.0/mu/(3.0*lambda + 2.0*mu)*SF(j)*SF(k)*wi(i)*D;
				StiffnessMatrix(9*j+5,9*k+4) += c*lambda/2.0/mu/(3.0*lambda + 2.0*mu)*SF(j)*SF(k)*wi(i)*D;

				StiffnessMatrix(9*j  ,9*k+3) += dSFdx(j)*SF(k)*(ay*bz+az*by)*wi(i)*D; 
				StiffnessMatrix(9*j  ,9*k+6) += dSFdy(j)*SF(k)*(ax*bz+az*bx)*wi(i)*D;
				StiffnessMatrix(9*j  ,9*k+8) += dSFdz(j)*SF(k)*(ax*by+ay*bx)*wi(i)*D;
				StiffnessMatrix(9*j+1,9*k+4) += dSFdy(j)*SF(k)*(ax*bz+az*bx)*wi(i)*D;
				StiffnessMatrix(9*j+1,9*k+6) += dSFdx(j)*SF(k)*(ay*bz+az*by)*wi(i)*D;
				StiffnessMatrix(9*j+1,9*k+7) += dSFdz(j)*SF(k)*(ax*by+ay*bx)*wi(i)*D;
				StiffnessMatrix(9*j+2,9*k+5) += dSFdz(j)*SF(k)*(ax*by+ay*bx)*wi(i)*D;
				StiffnessMatrix(9*j+2,9*k+8) += dSFdx(j)*SF(k)*(ay*bz+az*by)*wi(i)*D;
				StiffnessMatrix(9*j+2,9*k+7) += dSFdy(j)*SF(k)*(ax*bz+az*bx)*wi(i)*D;

				StiffnessMatrix(9*j+3,9*k  ) += dSFdx(k)*SF(j)*(ay*bz+az*by)*wi(i)*D; 
				StiffnessMatrix(9*j+6,9*k  ) += dSFdy(k)*SF(j)*(ax*bz+az*bx)*wi(i)*D;
				StiffnessMatrix(9*j+8,9*k  ) += dSFdz(k)*SF(j)*(ax*by+ay*bx)*wi(i)*D;
				StiffnessMatrix(9*j+4,9*k+1) += dSFdy(k)*SF(j)*(ax*bz+az*bx)*wi(i)*D;
				StiffnessMatrix(9*j+6,9*k+1) += dSFdx(k)*SF(j)*(ay*bz+az*by)*wi(i)*D;
				StiffnessMatrix(9*j+7,9*k+1) += dSFdz(k)*SF(j)*(ax*by+ay*bx)*wi(i)*D;
				StiffnessMatrix(9*j+5,9*k+2) += dSFdz(k)*SF(j)*(ax*by+ay*bx)*wi(i)*D;
				StiffnessMatrix(9*j+8,9*k+2) += dSFdx(k)*SF(j)*(ay*bz+az*by)*wi(i)*D;
				StiffnessMatrix(9*j+7,9*k+2) += dSFdy(k)*SF(j)*(ax*bz+az*bx)*wi(i)*D;
			}
		}
	}
}

// =======================================================================
//	return stiffness matrix 
// =======================================================================
const Matrix& PML3DGeneral::getTangentStiff()
{
    // calculate the effective stiffness matrix
    double dt = Domainptr->getDT();
    double cg = (eta*dt)/beta;
    EffectiveStiffnessMatrix.addMatrix(0.0, StiffnessMatrix, 1.0);
    EffectiveStiffnessMatrix.addMatrix(1.0, ImpedanceMatrix, cg);
    return EffectiveStiffnessMatrix;
}

// =======================================================================
//	return initial stiffness matrix 
// =======================================================================
const Matrix& PML3DGeneral::getInitialStiff()
{
	return this->getTangentStiff();
}
// =======================================================================
//	return damping matrix
// =======================================================================
const Matrix& PML3DGeneral::getDamp()
{
    return DampingMatrix;
}

// =======================================================================
//	return mass matrix
// =======================================================================
const Matrix& PML3DGeneral::getMass()
{
    return MassMatrix;
}

// =======================================================================
// Ressisting force
// =======================================================================
//get residual
const Vector& PML3DGeneral::getResistingForce()
{
	//
	// perform: R = K * u
	//
    static Vector disp(PML3DGeneral_NUM_DOF);
	int loc = 0;
	for (int i = 0; i < PML3DGeneral_NUM_NODES; i++) {
		const Vector& uNode = nodePointers[i]->getTrialDisp();
		for (int j = 0; j < 9; j++) {
			disp(loc) = uNode(j);
            loc++;
        }
	}
	resid.addMatrixVector(0.0, StiffnessMatrix, disp, 1.0);
	return resid;
}

// =======================================================================
// get residual with inertia terms
// =======================================================================
const Vector&
PML3DGeneral::getResistingForceIncInertia()
{
    // get the current velocity and acceleration
    static Vector disp(PML3DGeneral_NUM_DOF);
    static Vector vel(PML3DGeneral_NUM_DOF);
    static Vector accel(PML3DGeneral_NUM_DOF);
    int loc = 0;
    for (int i = 0; i < PML3DGeneral_NUM_NODES; i++) {
        const Vector& uNode = nodePointers[i]->getTrialDisp();
        const Vector& velNode = nodePointers[i]->getTrialVel();
        const Vector& accelNode = nodePointers[i]->getTrialAccel();
        for (int j = 0; j < 9; j++) {
            disp(loc) = uNode(j);
            vel(loc) = velNode(j);
            accel(loc) = accelNode(j);
            loc++;
        }
    }

    // calculate residual
    resid.addMatrixVector(0.0, StiffnessMatrix, disp, 1.0);
    resid.addMatrixVector(1.0, DampingMatrix, vel, 1.0);
    resid.addMatrixVector(1.0, MassMatrix, accel, 1.0);
    resid.addMatrixVector(1.0, ImpedanceMatrix, ubar, 1.0);

    return resid;        
}

// =======================================================================
// update
// =======================================================================
int PML3DGeneral::update(void)
{
    // get the current time
    dt = Domainptr->getDT();
    int loc = 0;
    double c1 = dt;
    double c2 = dt * dt * 0.5;
    double c3 = dt*dt*dt*((1.0/6.0)-eta);
    double c4 = dt*dt*dt*eta;
    for (int i = 0; i < PML3DGeneral_NUM_NODES; i++) {
        const Vector& uNode = nodePointers[i]->getDisp();
        const Vector& vNode = nodePointers[i]->getVel();
        const Vector& aNode = nodePointers[i]->getAccel();
        const Vector& atpdt = nodePointers[i]->getTrialAccel();
        for (int j = 0; j < 9; j++) {
            ubar(loc) = ubart(loc) + uNode(j)*c1 + vNode(j)*c2 + aNode(j)*c3 + atpdt(j)*c4; 
            loc++;
        }
    }
	return 0;
}

// =======================================================================
// commit state
// =======================================================================
int  PML3DGeneral::commitState()
{
	int success = 0;
	if ((success = this->Element::commitState()) != 0) {
		opserr << "PML3DGeneral::commitState () - failed in base class";
	}
	// set ubart to ubar
	for (int i = 0; i < PML3DGeneral_NUM_DOF; i++) {
		ubart(i) = ubar(i);
	}
	return success;
}

// =======================================================================
// revert to last commit 
// =======================================================================
int  PML3DGeneral::revertToLastCommit()
{
	int success = 0;

	// set ubar to ubart
	for (int i = 0; i < PML3DGeneral_NUM_DOF; i++) {
		ubar(i) = ubart(i);
	}

	return success;
}


// =======================================================================
// revert to start
// =======================================================================
int  PML3DGeneral::revertToStart()
{
	int success = 0;
	ubar.Zero();
	ubart.Zero();
	return success;
}


// =======================================================================
// get the number of external nodes
// =======================================================================
int  PML3DGeneral::getNumExternalNodes() const
{
	return PML3DGeneral_NUM_NODES;
}

// =======================================================================
// return connected external nodes
// =======================================================================
const ID& PML3DGeneral::getExternalNodes()
{
	return connectedExternalNodes;
}

// =======================================================================
// return node pointers
// =======================================================================
Node** PML3DGeneral::getNodePtrs(void)
{
	return nodePointers;
}

// =======================================================================
// return number of dofs
// =======================================================================
int  PML3DGeneral::getNumDOF()
{
	return PML3DGeneral_NUM_DOF;
}

// =======================================================================
// add load
// =======================================================================
int PML3DGeneral::addLoad(ElementalLoad* theLoad, double loadFactor)
{
	return -1;
}

// =======================================================================
// add zero load
// =======================================================================
void  PML3DGeneral::zeroLoad()
{
	return;
}

// =======================================================================
// senself
// =======================================================================
int  PML3DGeneral::sendSelf(int commitTag,
	Channel& theChannel)
{
    return 0;
}

// =======================================================================
// recvself
// =======================================================================
int  PML3DGeneral::recvSelf(int commitTag,
	Channel& theChannel,
	FEM_ObjectBroker& theBroker)
{
    return 0;
}


// =======================================================================
// display
// =======================================================================
int PML3DGeneral::displaySelf(Renderer& theViewer, int displayMode, float fact, const char** modes, int numMode)
{
    return 0;
}

// =======================================================================
// setresponse
// =======================================================================
Response* PML3DGeneral::setResponse(const char** argv, int argc, OPS_Stream& output)
{
	Response* theResponse = 0;
	return theResponse;
}

// =======================================================================
// getresponse
// =======================================================================
int PML3DGeneral::getResponse(int responseID, Information& eleInfo)
{
	return -1;
}

// =======================================================================
// set parameter
// =======================================================================
int PML3DGeneral::setParameter(const char** argv, int argc, Parameter& param)
{
	int res = -1;
	return res;
}

// =======================================================================
// update parameter
// =======================================================================
int PML3DGeneral::updateParameter(int parameterID, Information& info)
{
	int res = -1;
	return res;
}
// =======================================================================
// print
// =======================================================================
void  PML3DGeneral::Print(OPS_Stream &s, int flag) {
	return;
}
