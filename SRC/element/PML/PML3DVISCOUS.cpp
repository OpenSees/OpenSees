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

// Written by: Amin Pakzad, Pedro Arduino (parduino@uw.edu)
//
// Eight node PML3DVISCOUS element .. a c++ wrapper to fortran routine 
// provided by Wenyang Zhang (zwyll@ucla.edu), University of California, Los Angeles
//
// University of Washington, UC. Los Angeles, U.C. Berkeley, 12, 2020


#include "PML3DVISCOUS.h"

#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 

#include <OPS_Globals.h>
#include <ID.h> 
#include <Vector.h>
#include <Matrix.h>
#include <Element.h>
#include <Node.h>
#include <Domain.h>
#include <ErrorHandler.h>
#include <Renderer.h>
#include <ElementResponse.h>
#include <Parameter.h>
#include <ElementalLoad.h>
#include <NDMaterial.h>
#include <ElasticIsotropicMaterial.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <elementAPI.h>
#include <cstring>



void OPS_PML3DVISCOUS_Box(double* boxParams, double* eleCenter, double* Xref, double* Normal) {
	// calculate the Xref and Normal for the Box meshType
	// boxParams: XC, YC, ZC, L, W, H
	// XC, YC, ZC: center of the box at the surface
	// L, W, H: length, width, height of the box
	Xref[0] = 0.0;
	Xref[1] = 0.0;
	Xref[2] = 0.0;
	Normal[0] = 0.0;
	Normal[1] = 0.0;
	Normal[2] = 0.0;
	double xref, yref, zref, nx, ny, nz;
	xref = 0.0; yref = 0.0; zref = 0.0;
	nx = 0.0; ny = 0.0; nz = 0.0;

	// calculate the Normal
	double widthx  = boxParams[3]/2.;
	double widthy  = boxParams[4]/2.;
	double depth   = boxParams[5];
	double Xmesh = boxParams[0];
	double Ymesh = boxParams[1];
	double Zmesh = boxParams[2];

	double x1 = eleCenter[0]; 
	double x2 = eleCenter[1];
	double x3 = eleCenter[2];

	int EleType_arg = 0;

	if (x2 < Ymesh - widthy) {
		if (x1 < Xmesh - widthx) {
			if (x3 < Zmesh - depth) {
				EleType_arg = 15;
				nx = -1.0;
				ny = -1.0;
				nz = -1.0;
				xref = Xmesh - widthx;
				yref = Ymesh - widthy;
				zref = Zmesh - depth;
			}
			else {
				EleType_arg = 6;
				nx = -1.0;
				ny = -1.0;
				nz = 0.0;
				xref = Xmesh - widthx;
				yref = Ymesh - widthy;
				zref = Zmesh;
			}
		}
		else if (x1 < Xmesh + widthx) {
			if (x3 < Zmesh - depth) {
				EleType_arg = 11;
				nx = 0.0;
				ny = -1.0;
				nz = -1.0;
				xref = Xmesh;
				yref = Ymesh - widthy;
				zref = Zmesh - depth;
			}
			else {
				EleType_arg = 2;
				nx = 0.0;
				ny = -1.0;
				nz = 0.0;
				xref = Xmesh;
				yref = Ymesh - widthy;
				zref = Zmesh;
			}
		}
		else {
			if (x3 < Zmesh - depth) {
				EleType_arg = 16;
				nx = 1.0;
				ny = -1.0;
				nz = -1.0;
				xref = Xmesh + widthx;
				yref = Ymesh - widthy;
				zref = Zmesh - depth;
			}
			else {
				EleType_arg = 7;
				nx = 1.0;
				ny = -1.0;
				nz = 0.0;
				xref = Xmesh + widthx;
				yref = Ymesh - widthy;
				zref = Zmesh;
			}
		} 
	} else if (x2 < Ymesh + widthy) {
		if (x1 < Xmesh - widthx) {
			if (x3 < Zmesh - depth) {
				EleType_arg = 14;
				nx = -1.0;
				ny = 0.0;
				nz = -1.0;
				xref = Xmesh - widthx;
				yref = Ymesh;
				zref = Zmesh - depth;
			}
			else {
				EleType_arg = 5;
				nx = -1.0;
				ny = 0.0;
				nz = 0.0;
				xref = Xmesh - widthx;
				yref = Ymesh;
				zref = Zmesh;
			}
		}
		else if (x1 < Xmesh + widthx) {
			if (x3 < Zmesh - depth) {
				EleType_arg = 10;
				nx = 0.0;
				ny = 0.0;
				nz = -1.0;
				xref = Xmesh;
				yref = Ymesh;
				zref = Zmesh - depth;
			}
			else {
				EleType_arg = 1;
				nx = 0.0;
				ny = 0.0;
				nz = 0.0;
				xref = Xmesh;
				yref = Ymesh;
				zref = Zmesh;
			}
		}
		else {
			if (x3 < Zmesh - depth) {
				EleType_arg = 12;
				nx = 1.0;
				ny = 0.0;
				nz = -1.0;
				xref = Xmesh + widthx;
				yref = Ymesh;
				zref = Zmesh - depth;
			}
			else {
				EleType_arg = 3;
				nx = 1.0;
				ny = 0.0;
				nz = 0.0;
				xref = Xmesh + widthx;
				yref = Ymesh;
				zref = Zmesh;
			}
		}
	} else {
		if (x1 < Xmesh - widthx) {
			if (x3 < Zmesh - depth) {
				EleType_arg = 18;
				nx = -1.0;
				ny = 1.0;
				nz = -1.0;
				xref = Xmesh - widthx;
				yref = Ymesh + widthy;
				zref = Zmesh - depth;
			}
			else {
				EleType_arg = 9;
				nx = -1.0;
				ny = 1.0;
				nz = 0.0;
				xref = Xmesh - widthx;
				yref = Ymesh + widthy;
				zref = Zmesh;
			}
		}
		else if (x1 < Xmesh + widthx) {
			if (x3 < Zmesh - depth) {
				EleType_arg = 13;
				nx = 0.0;
				ny = 1.0;
				nz = -1.0;
				xref = Xmesh;
				yref = Ymesh + widthy;
				zref = Zmesh - depth;
			}
			else {
				EleType_arg = 4;
				nx = 0.0;
				ny = 1.0;
				nz = 0.0;
				xref = Xmesh;
				yref = Ymesh + widthy;
				zref = Zmesh;
			}
		}
		else {
			if (x3 < Zmesh - depth) {
				EleType_arg = 17;
				nx = 1.0;
				ny = 1.0;
				nz = -1.0;
				xref = Xmesh + widthx;
				yref = Ymesh + widthy;
				zref = Zmesh - depth;
			}
			else {
				EleType_arg = 8;
				nx = 1.0;
				ny = 1.0;
				nz = 0.0;
				xref = Xmesh + widthx;
				yref = Ymesh + widthy;
				zref = Zmesh;
			}
		}
	}

    // int tag = 1;
	// opserr << "Element Tag: " << tag << "\n";
	// opserr << "Element Type(" << tag << "): " << EleType_arg << "\n";
	// // opserr << "Element Center(" << tag << "): " << elementCenter[0] << " " << elementCenter[1] << " " << elementCenter[2] << "\n";
	// opserr << "Element Reference(" << tag << "): " << xref << " " << yref << " " << zref << "\n";
	// opserr << "Element Normal(" << tag << "): " << nx << " " << ny << " " << nz << "\n";

	// return the calculated Xref and Normal
	Xref[0] = xref;
	Xref[1] = yref;
	Xref[2] = zref;
	Normal[0] = nx;
	Normal[1] = ny;
	Normal[2] = nz;
}


// =======================================================================
// PML3DVISCOUS element tcl command
// =======================================================================
void* OPS_PML3DVISCOUS()
{
	// check if the total number of arguments passed is correct
	if (OPS_GetNumRemainingInputArgs() < (11)) {
		opserr << "WARNING insufficient arguments\n";
		opserr << "Want: element PML3DVISCOUS eleTag? Node1? Node2? Node3? Node4? Node5? Node6? Node7? Node 8? matTag? $PMLThickness\n";
		return 0;
	}

	// reading element tag and node numbers 
	int idata[9];
	int num = 9;
	if (OPS_GetIntInput(&num, idata) < 0) {
		opserr << "WARNING: invalid integer data : could be the tag or the node numbers \n";
		return 0;
	}
	// THIS SHOULD BE DONE IN setDomain()
	// calculate the number center of the element by averaging the node coordinates
	double elementCenter[3] = {0.0, 0.0, 0.0};
	/*
	// OPS_GetNodeCrd(int* nodeTag, int* sizeData, double* data);
	for (int i = 1; i < 9; i++) {
		int nodeTag = idata[i];
		double nodeCoord[3] = {0.0, 0.0, 0.0};
		num = 3;
		// REWRITE TO DO THIS IN setDomain of the element, not in the OPS_command
		if (OPS_GetNodeCrd(&nodeTag, &num, nodeCoord) < 0) {
			opserr << "WARNING: invalid node tag: " << nodeTag << endln;
			return 0;
		}
		for (int j = 0; j < 3; j++) {
			elementCenter[j] += nodeCoord[j];
		}
	}
	for (int j = 0; j < 3; j++) {
		elementCenter[j] /= 8;
	}
	*/

	// get the material tag
	int matTag;
	num = 1;
	if (OPS_GetIntInput(&num, &matTag) < 0) {
		opserr << "WARNING: invalid integer data: could be the material tag\n";
		return 0;
	}

	NDMaterial* mat = OPS_getNDMaterial(matTag);
	if (mat == 0) {
		opserr << "WARNING material not found\n";
		opserr << "material tag: " << matTag;
		opserr << "\nPML3DVISCOUS element: " << idata[0] << endln;
	}
	// check if the material is "ElasticIsotropicMaterial"
	if (strcmp(mat->getClassType(), "ElasticIsotropicMaterial") != 0) {
		opserr << "Error: PML3DVISCOUS element only supports ElasticIsotropicMaterial\n";
		opserr << "\tMaterial provided is of type: " << mat->getClassType() << endln;
		return 0;
	}

	double E = 0.0, nu = 0.0, rho = 0.0;

	const char* argv_E[] = {"E"};
	const char* argv_nu[] = {"nu"};

	Parameter paramE(1);
	Parameter paramNu(2);

	// This is the intended way to access material properties
	if (mat->setParameter(argv_E, 1, paramE) >= 0) {
		E = paramE.getValue();
	} else {
		opserr << "Warning: Material doesn't support E parameter\n";
	}

	if (mat->setParameter(argv_nu, 1, paramNu) >= 0) {
		nu = paramNu.getValue();
	} else {
		opserr << "Warning: Material doesn't support nu parameter\n";
	}

	// Always available in NDMaterial base class
	rho = mat->getRho();

	// opserr << "PML3DVISCOUS element: E=" << E << ", nu=" << nu << ", rho=" << rho << endln;



	// get the PML thickness
	double PMLThickness;
	num = 1;
	if (OPS_GetDoubleInput(&num, &PMLThickness) < 0) {
		opserr << "WARNING: invalid double data: could be the PMLThickness\n";
		return 0;
	}



	// reading the meshType string and it related parameters
	num = 1;
	if (OPS_GetNumRemainingInputArgs() < 1) {
		opserr << "Error: need meshType\n";
		opserr << "meshType could be: \"General\", \"Box\", \"Sphere\", \"Cylinder\"\n";
		return 0;
	}
	double Xref[3]={0.0, 0.0, 0.0}; // double Normal[]
	double Normal[3]={0.0, 0.0, 0.0}; // double Normal[]

	const char* meshType = OPS_GetString();
	if (strcmp(meshType, "General") == 0) {
		// reading the meshType string and it related parameters
		// this meshType needs to read the following parameters
		// Xref, Yref, Zref, N1, N2, N3
		num = 6;
		if (OPS_GetNumRemainingInputArgs() < 6) {
			opserr << "Error: need meshType parameters\n";
			opserr << "meshType parameters: Xref, Yref, Zref, N1, N2, N3\n";
			return 0;
		}
		num = 3;
		if (OPS_GetDoubleInput(&num, Xref) < 0) {
			opserr << "WARNING: invalid double data: could be Xref, Yref, Zref\n";
			return 0;
		}

		num = 3;
		if (OPS_GetDoubleInput(&num, Normal) < 0) {
			opserr << "WARNING: invalid double data: could be N1, N2, N3\n";
			return 0;
		}	
	}
	else if (strcmp(meshType, "Box") == 0 || strcmp(meshType, "box") == 0) {
		// reading the meshType string and it related parameters
		// this meshType needs to read the following parameters
		// XC, YC, ZC, L, W, H
		// XC, YC, ZC: center of the box at the surface
		// L, W, H: length, width, height of the box
		num = 6;
		if (OPS_GetNumRemainingInputArgs() < 6) {
			opserr << "Error: need meshType parameters\n";
			opserr << "meshType parameters: XC, YC, ZC, L, W, H\n";
			return 0;
		}
		double boxParams[6];
		num = 6;
		if (OPS_GetDoubleInput(&num, boxParams) < 0) {
			opserr << "WARNING: invalid double data: could be XC, YC, ZC, L, W, H\n";
			return 0;
		}
		// call the helper function to calculate the Xref and Normal
		OPS_PML3DVISCOUS_Box(boxParams, elementCenter, Xref, Normal);
	}
	else if (strcmp(meshType, "Sphere") == 0 || strcmp(meshType, "sphere") == 0) {
		// reading the meshType string and it related parameters
		// this meshType needs to read the following parameters
		// XC, YC, ZC, R
		// XC, YC, ZC: center of the sphere at the surface
		// R: radius of the sphere
		num = 4;
		if (OPS_GetNumRemainingInputArgs() < 4) {
			opserr << "Error: need meshType parameters\n";
			opserr << "meshType parameters: XC, YC, ZC, R\n";
			return 0;
		}
		double sphereParams[4];
		num = 4;
		if (OPS_GetDoubleInput(&num, sphereParams) < 0) {
			opserr << "WARNING: invalid double data: could be XC, YC, ZC, R\n";
			return 0;
		}
		// call the helper function to calculate the Xref and Normal
		// not implemented yet
		opserr << "Error: Sphere meshType not implemented yet\n";
		return 0;
	}
	else if (strcmp(meshType, "Cylinder") == 0 || strcmp(meshType, "cylinder") == 0) {
		// reading the meshType string and it related parameters
		// this meshType needs to read the following parameters
		// XC, YC, ZC, R, H, dir
		// XC, YC, ZC: center of the cylinder at the surface
		// R: radius of the cylinder
		// H: height of the cylinder
		// dir: direction of the cylinder
		num = 6;
		if (OPS_GetNumRemainingInputArgs() < 6) {
			opserr << "Error: need meshType parameters\n";
			opserr << "meshType parameters: XC, YC, ZC, R, H, dir\n";
			return 0;
		}
		double cylinderParams[6];
		num = 6;
		if (OPS_GetDoubleInput(&num, cylinderParams) < 0) {
			opserr << "WARNING: invalid double data: could be XC, YC, ZC, R, H, dir\n";
			return 0;
		}
		// call the helper function to calculate the Xref and Normal
		// not implemented yet
		opserr << "Error: Cylinder meshType not implemented yet\n";
		return 0;
	}
	else {
		opserr << "Error: meshType not recognized\n";
		opserr << "meshType could be: \"General\", \"Box\", \"Sphere\", \"Cylinder\"\n";
		return 0;
	}
	// cp_ref = SQRT(E *(1.d0-xnu)/rho/(1.d0+xnu)/(1.d0-2.d0*xnu))
	double Cp = sqrt(E * (1.0 - nu) / rho / (1.0 + nu) / (1.0 - 2.0 * nu));
	double m_coeff = 2.0;
	double R = 1.0e-8;
	double PML_b = PMLThickness / 1.0;
	double alpha_0 = 0.0;
	double beta_0 = 0.0;
	bool   ExplicitAlphaBeta = false;
	double gamma = 0.5;
	double beta = 0.25;
	double eta = 1.0/12.0;
	double keisi = 1.0/48.0;
	
	// now iterate overe the remaining arguments to read optional parameters
	while (OPS_GetNumRemainingInputArgs() > 0) {
		const char* option = OPS_GetString();
		if (strcmp(option, "-Newmark") == 0 || strcmp(option, "-newmark") == 0) {
			// reading Newmark parameters
			double Newmark[4];
			num = 4;
			if (OPS_GetDoubleInput(&num, Newmark) < 0) {
				opserr << "WARNING: invalid double data: could be Newmark parameters\n";
				opserr << "Newmark parameters: gamma, beta, dt, alpha_f\n";
				return 0;
			}
			gamma = Newmark[0];
			beta = Newmark[1];
			eta = Newmark[2];
			keisi = Newmark[3];
		}
		else if (strcmp(option, "-Cp") == 0 || strcmp(option, "-cp") == 0) {
			// reading the PML parameters
			num = 1;
			if (OPS_GetDoubleInput(&num, &Cp) < 0) {
				opserr << "WARNING: invalid double data: could be Cp\n";
				return 0;
			}
		}
		else if (strcmp(option, "-m") == 0 || strcmp(option, "-M") == 0) {
			// reading the PML parameters
			num = 1;
			if (OPS_GetDoubleInput(&num, &m_coeff) < 0) {
				opserr << "WARNING: invalid double data: could be m\n";
				return 0;
			}
		}
		else if (strcmp(option, "-R") == 0 || strcmp(option, "-r") == 0) {
			// reading the PML parameters
			num = 1;
			if (OPS_GetDoubleInput(&num, &R) < 0) {
				opserr << "WARNING: invalid double data: could be R\n";
				return 0;
			}
		}
		else if (strcmp(option, "-alphabeta") == 0) {
			// reading the PML parameters
			double alphabeta[2];
			num = 2;
			if (OPS_GetDoubleInput(&num, alphabeta) < 0) {
				opserr << "WARNING: invalid double data: could be alphabeta\n";
				return 0;
			}
			alpha_0 = alphabeta[0];
			beta_0 = alphabeta[1];
			ExplicitAlphaBeta = true;
		}
		else {
			opserr << "WARNING: unknown option: " << option << endln;
			return 0;
		}
	}

		
	if (!ExplicitAlphaBeta) {
		
		alpha_0 = ((m_coeff + 1) * PML_b) / (2.0 * PML_b)*log10(1.0 / R);
		beta_0  = ((m_coeff + 1) * Cp)    / (2.0 * PML_b)*log10(1.0 / R);
	} else {
		opserr << "Warning: PML Element(" << idata[0] << ") is using user defined alpha and beta\n"; 
	}

	// opserr << "PML3DVISCOUS element, tag: " << idata[0] << endln;
	// opserr << "Node 1: " << idata[1] << endln;
	// opserr << "Node 2: " << idata[2] << endln;
	// opserr << "Node 3: " << idata[3] << endln;
	// opserr << "Node 4: " << idata[4] << endln;
	// opserr << "Node 5: " << idata[5] << endln;
	// opserr << "Node 6: " << idata[6] << endln;
	// opserr << "Node 7: " << idata[7] << endln;
	// opserr << "Node 8: " << idata[8] << endln;
	// opserr << "Element center: " << elementCenter[0] << ", " << elementCenter[1] << ", " << elementCenter[2] << endln;
	// opserr << "PML Thickness: " << PMLThickness << endln;
	// opserr << "Material tag: " << matTag << endln;
	// opserr << "Material type: " << mat->getClassType() << endln;
	// opserr << "Material information: " << endln;
	// opserr << "\t E: " << E << endln;
	// opserr << "\t nu: " << nu << endln;	
	// opserr << "\t rho: " << rho << endln;
	// opserr << "Mesh Type: " << meshType << endln;
	// opserr << "Xref: " << Xref[0] << ", " << Xref[1] << ", " << Xref[2] << endln;
	// opserr << "Normal: " << Normal[0] << ", " << Normal[1] << ", " << Normal[2] << endln;
	// opserr << "Cp: " << Cp << endln;
	// opserr << "m: " << m_coeff << endln;
	// opserr << "R: " << R << endln;
	// opserr << "alpha_0: " << alpha_0 << endln;
	// opserr << "beta_0: " << beta_0 << endln;
	// opserr << "ExplicitAlphaBeta: " << ExplicitAlphaBeta << endln;
	

	return new PML3DVISCOUS(idata[0], &idata[1], 
							mat, PMLThickness,
							Xref, Normal, 
							alpha_0, beta_0, ExplicitAlphaBeta,
							Cp, m_coeff, R,
							gamma, beta, eta, keisi); 
}


// =======================================================================
// static data
// =======================================================================
double  PML3DVISCOUS::gamma = 0.;
double  PML3DVISCOUS::beta = 0.;
double  PML3DVISCOUS::eta = 0.;
double  PML3DVISCOUS::keisi = 0.;
double  PML3DVISCOUS::dt = 0.;
int     PML3DVISCOUS::eleCount = 0;
int     PML3DVISCOUS::ComputedEleTag = -1;

double PML3DVISCOUS::M[PML3DVISCOUS_NUM_DOF * PML3DVISCOUS_NUM_DOF] = {0};
double PML3DVISCOUS::C[PML3DVISCOUS_NUM_DOF * PML3DVISCOUS_NUM_DOF] = {0};
double PML3DVISCOUS::K[PML3DVISCOUS_NUM_DOF * PML3DVISCOUS_NUM_DOF] = {0};
double PML3DVISCOUS::G[PML3DVISCOUS_NUM_DOF * PML3DVISCOUS_NUM_DOF] = {0};
double PML3DVISCOUS::H[PML3DVISCOUS_NUM_DOF * PML3DVISCOUS_NUM_DOF] = {0};
double PML3DVISCOUS::Keff[PML3DVISCOUS_NUM_DOF * PML3DVISCOUS_NUM_DOF] = {0};


// Matrix  PML3DVISCOUS::tangent(PML3DVISCOUS_NUM_DOF, PML3DVISCOUS_NUM_DOF);
// Matrix  PML3DVISCOUS::mass(PML3DVISCOUS_NUM_DOF, PML3DVISCOUS_NUM_DOF);
// Matrix  PML3DVISCOUS::damping(PML3DVISCOUS_NUM_DOF, PML3DVISCOUS_NUM_DOF);
Matrix  PML3DVISCOUS::tangent(K, PML3DVISCOUS_NUM_DOF, PML3DVISCOUS_NUM_DOF);
Matrix  PML3DVISCOUS::mass(M, PML3DVISCOUS_NUM_DOF, PML3DVISCOUS_NUM_DOF);
Matrix  PML3DVISCOUS::damping(C, PML3DVISCOUS_NUM_DOF, PML3DVISCOUS_NUM_DOF);
Matrix  PML3DVISCOUS::Gmat(G, PML3DVISCOUS_NUM_DOF, PML3DVISCOUS_NUM_DOF);
Matrix  PML3DVISCOUS::Hmat(H, PML3DVISCOUS_NUM_DOF, PML3DVISCOUS_NUM_DOF);
Matrix  PML3DVISCOUS::keffmat(Keff, PML3DVISCOUS_NUM_DOF, PML3DVISCOUS_NUM_DOF);
Vector  PML3DVISCOUS::resid(PML3DVISCOUS_NUM_DOF);
// int     PML3DVISCOUS::numberOfElements = 0;


// =======================================================================
// null constructor
// =======================================================================
PML3DVISCOUS::PML3DVISCOUS()
	:Element(0, ELE_TAG_PML3DVISCOUS),
	connectedExternalNodes(PML3DVISCOUS_NUM_NODES),
	ubar(PML3DVISCOUS_NUM_DOF),
	ubart(PML3DVISCOUS_NUM_DOF),
    ubarbar(PML3DVISCOUS_NUM_DOF),
    ubarbart(PML3DVISCOUS_NUM_DOF)
{
	for (int i = 0; i < PML3DVISCOUS_NUM_NODES; i++) {
		nodePointers[i] = 0;
	}
	dt = 0;
	ubar.Zero();
	ubart.Zero();
    ubarbar.Zero();
    ubarbart.Zero();
	updateflag = 0;
	update_dt = 0;
	gamma = 0;
	beta = 0;
	eta = 0;
	keisi = 0;
	tangent.Zero();
	mass.Zero();
	damping.Zero();
	Gmat.Zero();
	Hmat.Zero();
	keffmat.Zero();
	theMaterial = 0;
	PML_L = 0;
	cp_ref = 0;
	m_coeff = 0;
	R_coeff = 0;
	alpha0 = 0;
	beta0 = 0;
	ExplicitAlphaBeta = false;
	E = 0.0;
	nu = 0.0;
	rho = 0.0;


}

// =======================================================================
// Full constructor
// =======================================================================
PML3DVISCOUS::PML3DVISCOUS(int tag, int* nodeTags,
						   NDMaterial* theMat, 
						   double PMLThickness,
						   double* Xref, double* Normal,
						   double alpha_0, double beta_0, bool explicitAB,
						   double Cp, double m_coeff, double R,
						   double gammaN, double betaN, double etaN, double keisiN
	):Element(tag, ELE_TAG_PML3DVISCOUS),
	connectedExternalNodes(PML3DVISCOUS_NUM_NODES),
	ubar(PML3DVISCOUS_NUM_DOF),
	ubart(PML3DVISCOUS_NUM_DOF),
    ubarbar(PML3DVISCOUS_NUM_DOF),
    ubarbart(PML3DVISCOUS_NUM_DOF),
	theMaterial(theMat), PML_L(PMLThickness),
	cp_ref(Cp), m_coeff(m_coeff), R_coeff(R),
	alpha0(alpha_0), beta0(beta_0), ExplicitAlphaBeta(explicitAB)
{
	eleCount++;
	if (eleCount == 1) {
		opserr << "Perfectly Matched Layer 3D (PMLVISCOUS) element -  Written: W. Zhang, E. Taciroglu, A. Pakzad, P. Arduino, UCLA, U.Washington\n ";
		tangent.Zero();
		mass.Zero();
		damping.Zero();
		Gmat.Zero();
		Hmat.Zero();
		keffmat.Zero();
	}
	// initialize node pointers
	for (int i = 0; i < PML3DVISCOUS_NUM_NODES; i++) {
		connectedExternalNodes(i) = nodeTags[i];
		nodePointers[i] = 0;
	}

	// initialize Newmark parameters
	gamma = gammaN;
	beta  = betaN;
	eta   = etaN;
	keisi = keisiN;

	// initialize the ubar and ubart vectors to zero
	ubart.Zero();
	ubar.Zero();
    ubarbar.Zero();
    ubarbart.Zero();
	updateflag = 0;
	update_dt = 0;


	xref = Xref[0];
	yref = Xref[1];
	zref = Xref[2];

	nx = Normal[0];
	ny = Normal[1];
	nz = Normal[2];

	Parameter paramE(1);
	Parameter paramNu(2);
	const char* argv_E[] = {"E"};
	const char* argv_nu[] = {"nu"};
		if (strcmp(theMaterial->getClassType(), "ElasticIsotropicMaterial") == 0) {
		int res = 0;
		res = theMaterial->setParameter(argv_nu, 1, paramNu);
		res = theMaterial->setParameter(argv_E, 1, paramE);
        if (res >= 0) {
			nu = paramNu.getValue();
			E = paramE.getValue();
		} else {
			opserr << "Error: PML3DVISCOUS element only supports ElasticIsotropicMaterial\n";
			opserr << "\tMaterial provided is of type: " << theMaterial->getClassType() << endln;
			return;
		}
	}
	rho = theMaterial->getRho();
	// opserr << "PML3DVISCOUS element(Domain): E=" << E << ", nu=" << nu << ", rho=" << theMaterial->getRho() << endln;

	// print out the information
	// opserr << "PML3DVISCOUS element, tag: " << this->getTag() << endln;
	// opserr << "Node 1: " << connectedExternalNodes(0) << endln;
	// opserr << "Node 2: " << connectedExternalNodes(1) << endln;
	// opserr << "Node 3: " << connectedExternalNodes(2) << endln;
	// opserr << "Node 4: " << connectedExternalNodes(3) << endln;
	// opserr << "Node 5: " << connectedExternalNodes(4) << endln;
	// opserr << "Node 6: " << connectedExternalNodes(5) << endln;
	// opserr << "Node 7: " << connectedExternalNodes(6) << endln;
	// opserr << "Node 8: " << connectedExternalNodes(7) << endln;
	// opserr << "PML Thickness: " << PML_L << endln;
	// opserr << "Material tag: " << theMaterial->getTag() << endln;
	// opserr << "Material type: " << theMaterial->getClassType() << endln;
	// opserr << "Material information: " << endln;
	// opserr << "\t E: " << theMaterial->getElasticModulus() << endln;
	// opserr << "\t nu: " << theMaterial->getPoissonsRatio() << endln;
	// opserr << "\t rho: " << theMaterial->getRho() << endln;
	// opserr << "Xref: " << xref << ", " << yref << ", " << zref << endln;
	// opserr << "Normal: " << nx << ", " << ny << ", " << nz << endln;
	// opserr << "Cp: " << cp_ref << endln;
	// opserr << "m: " << m_coeff << endln;
	// opserr << "R: " << R_coeff << endln;
	// opserr << "alpha_0: " << alpha0 << endln;
	// opserr << "beta_0: " << beta0 << endln;
	// if (ExplicitAlphaBeta) {
	// 	opserr << "ExplicitAlphaBeta: " << "True" << endln;
	// } else {
	// 	opserr << "ExplicitAlphaBeta: " << "False" << endln;
	// }
	// opserr << "Newmark parameters: " << endln;
	// opserr << "\t gamma: " << gamma << endln;
	// opserr << "\t beta: " << beta << endln;
	// opserr << "\t eta: " << eta << endln;
	// opserr << "\t keisi: " << keisi << endln;
	// opserr << "Rayleigh Damping parameters: " << endln;
	// opserr << "\t alphaM: " << alphaM << endln;
	// opserr << "\t betaK: " << betaK << endln;
	// opserr << "\t betaK0: " << betaK0 << endln;
	// opserr << "\t betaKc: " << betaKc << endln;

}

// =======================================================================
//  destructor
// ======================================================================= 
PML3DVISCOUS::~PML3DVISCOUS()
{

}

// =======================================================================
// Set Domain
// =======================================================================
void  PML3DVISCOUS::setDomain(Domain* theDomain)
{

	Domainptr = theDomain;

	// node pointers
	for (int i = 0; i < PML3DVISCOUS_NUM_NODES; i++)
		nodePointers[i] = theDomain->getNode(connectedExternalNodes(i));

	this->DomainComponent::setDomain(theDomain);
	for (int i = 0; i < PML3DVISCOUS_NUM_NODES; i++) {
		const Vector& loc = nodePointers[i]->getCrds();
		coords[i * 3] = loc(0);
		coords[i * 3 + 1] = loc(1);
		coords[i * 3 + 2] = loc(2);
	}
}
// =======================================================================
// Calculate the matricies
// =======================================================================
void PML3DVISCOUS::calculateMatrices()
{
	// create the coordinate vectors
	int NDOFEL = PML3DVISCOUS_NUM_DOF;
	int NPROPS = 13;
	int MCRD = 3;
	int NNODE = 8;
	int LFLAGS = 12;
	for (int i = 0; i < PML3DVISCOUS_NUM_DOF*PML3DVISCOUS_NUM_DOF; i++) {
		C[i] = 0.0;
		K[i] = 0.0;
		M[i] = 0.0;
		G[i] = 0.0;
		H[i] = 0.0;
	}

	double betarayleigh = (betaK0 > betaK) ? betaK0 : betaK;
	betarayleigh = (betaKc > betarayleigh) ? betaKc : betarayleigh;
	double props[16];
	props[0] = E;
	props[1] = nu;
	props[2] = rho;
	props[3] = 6.0;
	props[4] = PML_L;
	props[5] = xref;
	props[6] = yref;
	props[7] = zref;
	props[8] = nx;
	props[9] = ny;
	props[10] = nz;
	props[11] = alpha0;
	props[12] = beta0;
	props[13] = alphaM;	
	props[14] = betarayleigh;
	props[15] = m_coeff;
	
	pml3d_(M, C, K, G, H, &NDOFEL, props, coords, &MCRD, &NNODE, &LFLAGS);
	ComputedEleTag = this->getTag();
}



// =======================================================================
// update
// =======================================================================
int PML3DVISCOUS::update(void)
{
	dt = Domainptr->getDT();
	// opserr << "dt = " << dt << "\n";	
	// get u, v, a from nodes and calculate the ubar vector
	int loc = 0;
	double c1 = dt;
	double c2 = dt * dt * 0.5;
	double c3 = dt*dt*dt*((1.0/6.0)-eta);
	double c4 = dt*dt*dt*eta;
	for (int i = 0; i < PML3DVISCOUS_NUM_NODES; i++) {
		const Vector& uNode = nodePointers[i]->getDisp();
		const Vector& vNode = nodePointers[i]->getVel();
		const Vector& aNode = nodePointers[i]->getAccel();
		const Vector& atpdt = nodePointers[i]->getTrialAccel();
		for (int j = 0; j < 9; j++) {
			ubar(loc) = ubart(loc) + uNode(j)*c1 + vNode(j)*c2 + aNode(j)*c3 + atpdt(j)*c4; 
			loc++;
		}
	}

    loc = 0;
    c1 = dt;
    c2 = dt * dt * 0.5;
    c3 = dt * dt * dt /6.0;
    c4 = dt * dt * dt * dt * (1.0/24.0 - keisi);
    int c5 = dt * dt * dt * dt * keisi;
    for (int i = 0; i < PML3DVISCOUS_NUM_NODES; i++) {
        const Vector& uNode = nodePointers[i]->getDisp();
        const Vector& vNode = nodePointers[i]->getVel();
        const Vector& aNode = nodePointers[i]->getAccel();
        const Vector& atpdt = nodePointers[i]->getTrialAccel();
        for (int j = 0; j < 9; j++) {
            ubarbar(loc) = ubarbart(loc) + ubart(loc)*c1 + uNode(j)*c2 + vNode(j)*c3 + aNode(j)*c4 + atpdt(j)*c5; 
            loc++;
        }
    }

	return 0;
}

// =======================================================================
//	return stiffness matrix 
// =======================================================================
const Matrix& PML3DVISCOUS::getTangentStiff()
{
	if (ComputedEleTag != this->getTag()) {
		this->calculateMatrices();
	} 

	// check if the dt is changed to update the tangent stiffness matrix
	double cg = eta*dt/beta;
    double ch = dt * dt * keisi/beta;
	//keff = k + cg*g( k and g are symmetric matrices)
	for (int i = 0; i < PML3DVISCOUS_NUM_DOF*PML3DVISCOUS_NUM_DOF; i++) {
		Keff[i] = K[i] + cg*G[i] + ch*H[i];
	}
	return keffmat;
}

// =======================================================================
//	return initial stiffness matrix 
// =======================================================================
const Matrix& PML3DVISCOUS::getInitialStiff()
{
	return this->getTangentStiff();
}

// =======================================================================
//	return mass matrix
// =======================================================================
const Matrix& PML3DVISCOUS::getMass()
{
	if (ComputedEleTag != this->getTag()) {
		this->calculateMatrices();
	}
	return mass;
}

// =======================================================================
//	return damping matrix
// =======================================================================
const Matrix& PML3DVISCOUS::getDamp()
{
	if (ComputedEleTag != this->getTag()) {
		this->calculateMatrices();
	}
	return damping;
}

// =======================================================================
// Ressisting force
// =======================================================================
//get residual
const Vector& PML3DVISCOUS::getResistingForce()
{
	if (ComputedEleTag != this->getTag()) {
		this->calculateMatrices();
	}
	int numNodalDOF = 9;
	static Vector theVector(PML3DVISCOUS_NUM_DOF);


	//
	// perform: R = K * u
	//

	int loc = 0;
	for (int i = 0; i < PML3DVISCOUS_NUM_NODES; i++) {
		const Vector& uNode = nodePointers[i]->getTrialDisp();
		for (int j = 0; j < numNodalDOF; j++)
			theVector(loc++) = uNode(j);
	}
	resid.addMatrixVector(0.0, tangent, theVector, 1.0);
	return resid;
}


// =======================================================================
//
// =======================================================================
//get residual with inertia terms
const Vector&
PML3DVISCOUS::getResistingForceIncInertia()
{
	if (ComputedEleTag != this->getTag()) {
		this->calculateMatrices();
	}
    // R += K*u
	static Vector theVector(PML3DVISCOUS_NUM_DOF);

	int loc = 0;
	for (int i = 0; i < PML3DVISCOUS_NUM_NODES; i++) {
		const Vector& uNode = nodePointers[i]->getTrialDisp();
		for (int j = 0; j < 9; j++)
			theVector(loc++) = uNode(j);
	}
	resid.addMatrixVector(0.0, tangent, theVector, 1.0);



	// R += M*a
	loc = 0;
	Node** theNodes = this->getNodePtrs();
	for (int i = 0; i < PML3DVISCOUS_NUM_NODES; i++) {
		const Vector& acc = theNodes[i]->getTrialAccel();
		for (int j = 0; j < 9; j++) {
			theVector(loc++) = acc(j);
		}
	}
	resid.addMatrixVector(1.0, mass, theVector, 1.0);

	// R += C*v
	loc = 0;
	for (int i = 0; i < PML3DVISCOUS_NUM_NODES; i++) {
		const Vector& vel = theNodes[i]->getTrialVel();
		for (int j = 0; j < 9; j++) {
			theVector(loc++) = vel[j];
		}
	}
	resid.addMatrixVector(1.0, damping, theVector, 1.0);


	// R += G*ubar
	resid.addMatrixVector(1.0, Gmat, ubar, 1.0);

    // R += H*ubarbar
    resid.addMatrixVector(1.0, Hmat, ubarbar, 1.0);
	return resid;
}

// =======================================================================
// get the number of external nodes
// =======================================================================
int  PML3DVISCOUS::getNumExternalNodes() const
{
	return PML3DVISCOUS_NUM_NODES;
}

// =======================================================================
// return connected external nodes
// =======================================================================
const ID& PML3DVISCOUS::getExternalNodes()
{
	return connectedExternalNodes;
}

// =======================================================================
// return node pointers
// =======================================================================
Node** PML3DVISCOUS::getNodePtrs(void)
{
	return nodePointers;
}

// =======================================================================
// return number of dofs
// =======================================================================
int  PML3DVISCOUS::getNumDOF()
{
	return PML3DVISCOUS_NUM_DOF;
}

// =======================================================================
// commit state
// =======================================================================
int  PML3DVISCOUS::commitState()
{
	int success = 0;
	if ((success = this->Element::commitState()) != 0) {
		opserr << "PML3DVISCOUS::commitState () - failed in base class";
	}

	// set ubart to ubar
	for (int i = 0; i < PML3DVISCOUS_NUM_DOF; i++) {
		ubart(i) = ubar(i);
	}

	updateflag = 0;
	return success;
}

// =======================================================================
// revert to last commit 
// =======================================================================
int  PML3DVISCOUS::revertToLastCommit()
{
	int success = 0;

	// set ubar to ubart
	for (int i = 0; i < PML3DVISCOUS_NUM_DOF; i++) {
		ubar(i) = ubart(i);
        ubarbar(i) = ubarbart(i);
	}

	return success;
}

// =======================================================================
// revert to start
// =======================================================================
int  PML3DVISCOUS::revertToStart()
{
	int success = 0;

	// set ubar and ubart to zero
	for (int i = 0; i < PML3DVISCOUS_NUM_DOF; i++) {
		ubar(i) = 0.0;
		ubart(i) = 0.0;
        ubarbar(i) = 0.0;
        ubarbart(i) = 0.0;
	}

	return success;
}

// =======================================================================
// add load
// =======================================================================
int PML3DVISCOUS::addLoad(ElementalLoad* theLoad, double loadFactor)
{
	return -1;
}

// =======================================================================
// add zero load
// =======================================================================
void  PML3DVISCOUS::zeroLoad()
{
	return;
}

// =======================================================================
// senself
// =======================================================================
int  PML3DVISCOUS::sendSelf(int commitTag,
	Channel& theChannel)
{
	return 0;
	// int res = 0;

	// // note: we don't check for dataTag == 0 for Element
	// // objects as that is taken care of in a commit by the Domain
	// // object - don't want to have to do the check if sending data
	// int dataTag = this->getDbTag();

	// // PML3DVISCOUS packs its data into a Vector and sends this to theChannel
	// // along with its dbTag and the commitTag passed in the arguments
	// static Vector data(PML3DVISCOUS_NUM_PROPS + 4);
	// data(0) = this->getTag();

	// for (int ii = 1; ii <= PML3DVISCOUS_NUM_PROPS; ii++) {
	// 	data(ii) = props[ii - 1];
	// }
	// data(PML3DVISCOUS_NUM_PROPS+1) = eta;
	// data(PML3DVISCOUS_NUM_PROPS+2) = beta;
	// data(PML3DVISCOUS_NUM_PROPS+3) = gamma;

	// res += theChannel.sendVector(dataTag, commitTag, data);
	// if (res < 0) {
	// 	opserr << "WARNING PML3DVISCOUS::sendSelf() - " << this->getTag() << " failed to send Vector\n";
	// 	return res;
	// }


	// // PML3DVISCOUS then sends the tags of its four nodes
	// res += theChannel.sendID(dataTag, commitTag, connectedExternalNodes);
	// if (res < 0) {
	// 	opserr << "WARNING PML3DVISCOUS::sendSelf() - " << this->getTag() << " failed to send ID\n";
	// 	return res;
	// }

	// return res;
}

// =======================================================================
// recvself
// =======================================================================
int  PML3DVISCOUS::recvSelf(int commitTag,
	Channel& theChannel,
	FEM_ObjectBroker& theBroker)
{
	return 0;
	// int res = 0;

	// int dataTag = this->getDbTag();

	// // PML3DVISCOUS creates a Vector, receives the Vector and then sets the 
	// // internal data with the data in the Vector
	// static Vector data(PML3DVISCOUS_NUM_PROPS + 4);
	// res += theChannel.recvVector(dataTag, commitTag, data);
	// if (res < 0) {
	// 	opserr << "WARNING PML3DVISCOUS::recvSelf() - failed to receive Vector\n";
	// 	return res;
	// }

	// this->setTag((int)data(0));


	// for (int ii = 1; ii <= PML3DVISCOUS_NUM_PROPS; ii++) {
	// 	props[ii - 1] = data(ii);
	// }

	// eta   = data(PML3DVISCOUS_NUM_PROPS+1);
	// beta  = data(PML3DVISCOUS_NUM_PROPS+2);
	// gamma = data(PML3DVISCOUS_NUM_PROPS+3);

	// // PML3DVISCOUS now receives the tags of its four external nodes
	// res += theChannel.recvID(dataTag, commitTag, connectedExternalNodes);
	// if (res < 0) {
	// 	opserr << "WARNING PML3DVISCOUS::recvSelf() - " << this->getTag() << " failed to receive ID\n";
	// 	return res;
	// }

	// return res;
}


// =======================================================================
// display
// =======================================================================
int PML3DVISCOUS::displaySelf(Renderer& theViewer, int displayMode, float fact, const char** modes, int numMode)
{
	// Get the end point display coords
	static Vector v1(3);
	static Vector v2(3);
	static Vector v3(3);
	static Vector v4(3);
	static Vector v5(3);
	static Vector v6(3);
	static Vector v7(3);
	static Vector v8(3);
	nodePointers[0]->getDisplayCrds(v1, fact, displayMode);
	nodePointers[1]->getDisplayCrds(v2, fact, displayMode);
	nodePointers[2]->getDisplayCrds(v3, fact, displayMode);
	nodePointers[3]->getDisplayCrds(v4, fact, displayMode);
	nodePointers[4]->getDisplayCrds(v5, fact, displayMode);
	nodePointers[5]->getDisplayCrds(v6, fact, displayMode);
	nodePointers[6]->getDisplayCrds(v7, fact, displayMode);
	nodePointers[7]->getDisplayCrds(v8, fact, displayMode);

	// place values in coords matrix
	static Matrix Coords(8, 3);
	for (int i = 0; i < 3; i++) {
		Coords(0, i) = v1(i);
		Coords(1, i) = v2(i);
		Coords(2, i) = v3(i);
		Coords(3, i) = v4(i);
		Coords(4, i) = v5(i);
		Coords(5, i) = v6(i);
		Coords(6, i) = v7(i);
		Coords(7, i) = v8(i);
	}

	// fill RGB vector
	static Vector values(8);
	for (int i = 0; i < 8; i++)
		values(i) = 1.0;

	// draw the cube
	return theViewer.drawCube(Coords, values, this->getTag());
}

// =======================================================================
// setresponse
// =======================================================================
Response* PML3DVISCOUS::setResponse(const char** argv, int argc, OPS_Stream& output)
{
	Response* theResponse = 0;

	// char outputData[32];

	// output.tag("ElementOutput");
	// output.attr("eleType", "PML3DVISCOUS");
	// output.attr("eleTag", this->getTag());
	// for (int i = 1; i <= 8; i++) {
	// 	sprintf(outputData, "node%d", i);
	// 	output.attr(outputData, nodePointers[i - 1]->getTag());
	// }

	// if (strcmp(argv[0], "force") == 0 || strcmp(argv[0], "forces") == 0) {

	// 	for (int i = 1; i <= 8; i++) {
	// 		sprintf(outputData, "P1_%d", i);
	// 		output.tag("ResponseType", outputData);
	// 		sprintf(outputData, "P2_%d", i);
	// 		output.tag("ResponseType", outputData);
	// 		sprintf(outputData, "P3_%d", i);
	// 		output.tag("ResponseType", outputData);
	// 	}

	// 	theResponse = new ElementResponse(this, 1, resid);
	// }
	// output.endTag(); // ElementOutput
	return theResponse;
}

// =======================================================================
// getresponse
// =======================================================================
int PML3DVISCOUS::getResponse(int responseID, Information& eleInfo)
{
	// static Vector stresses(48);

	// if (responseID == 1)
	// 	return eleInfo.setVector(this->getResistingForce());

	return -1;
}

// =======================================================================
// set parameter
// =======================================================================
int PML3DVISCOUS::setParameter(const char** argv, int argc, Parameter& param)
{
	int res = -1;
	return res;
}

// =======================================================================
// update parameter
// =======================================================================
int PML3DVISCOUS::updateParameter(int parameterID, Information& info)
{
	int res = -1;
	return res;
}

// =======================================================================
// print
// =======================================================================
void  PML3DVISCOUS::Print(OPS_Stream &s, int flag) {
	
  	if (flag == OPS_PRINT_CURRENTSTATE) {
		s << "Element: " << this->getTag() << endln;
		s << "type: PML3DVISCOUS \n";
		s << "Nodes: " << connectedExternalNodes;
		s << "eta: " << eta << " beta: " << beta << " gamma: " << gamma << endln;
		s << endln;
		s << "Resisting Force (no inertia): " << this->getResistingForce();
  	}
    
 	 if (flag == OPS_PRINT_PRINTMODEL_JSON) {
		s << "\t\t\t{";
		s << "\"name\": " << this->getTag() << ", ";
		s << "\"type\": \"PML3DVISCOUS\", ";
		s << "\"nodes\": [" << connectedExternalNodes(0) << ", ";
		for (int i = 1; i < 7; i++)
		s << connectedExternalNodes(i) << ", ";
		s << connectedExternalNodes(7) << "], ";
  	}
	return;
}
