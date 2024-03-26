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

// Code developed by: Cristian V. Miculaș  (github user name: cvmiculas)
// Element conceptualization: Cristian V. Miculaș (cristian.miculas@uc.pt), Ricardo J. Costa (rjcosta@dec.uc.pt) and Luís Simões da Silva (luisss@dec.uc.pt)
// Affiliation: Civil Engineering Department, ISISE, University of Coimbra, Portugal
// Acknowledgements: This work has been supported in part by national funds through FCT – Foundation for Science and Technology, Portugal, under grant agreement SFRH/BD/138151/2018 awarded to Cristian V. Miculaş

// Created: 04.02.2024
// Revised: 

// Description:
// This file contains the class implementation for the 3D beam-column-joint element object designed for the 3D innovative plug-and-play steel tubular joint configuration proposed within the INNO3DJOINTS project (https://ec.europa.eu/info/funding-tenders/opportunities/portal/screen/how-to-participate/org-details/960532413/project/749959/program/31061225/details).
// This element has 5 external nodes (6 DOFs/node) and 4 internal nodes (1 DOF/node), resulting in a total of 34 DOFs.
// This element can be viewed as a 2D plate in the 3D space.
// This element has 32 componenets (0D elements), each allowing for a different uniaxial material tag.

// References:
//
// 1. C.V. Miculaş, Innovative plug-and-play joints for hybrid tubular constructions (Ph.D. thesis), University of Coimbra, Portugal, 2023, https://estudogeral.uc.pt/handle/10316/110990
//
// 2. C. V. Miculaş, R. J. Costa, L. S. da Silva, R. Simões, H. Craveiro, T. Tankova, 3D macro-element for innovative plug-and-play joints, J. Constructional Steel Research 214 (2024), https://doi.org/10.1016/j.jcsr.2023.108436
//
// 3. C.V. Miculaş, R.J. Costa, L. Simões da Silva, R. Simões, H. Craveiro, T. Tankova, Macro-modelling of the three-dimensional interaction between the faces of a steel tubular column joint, in: F. Di Trapani, C. Demartino, G.C. Marano, G. Monti (Eds.), Proceedings of the 2022 Eurasian OpenSees Days, Springer Nature Switzerland, Cham, 2023, pp. 408–422, http://dx.doi.org/10.1007/978-3-031-30125-4_37
//


#include <Inno3DPnPJoint.h>

#include <Domain.h>
#include <Node.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Renderer.h>
#include <MatrixUtil.h>

#include <UniaxialMaterial.h>
#include <string.h>
#include <math.h>
#include <ElementResponse.h>
#include <elementAPI.h>
#include <stdlib.h>

#ifdef _USRDLL
#define OPS_Export extern "C" _declspec(dllexport)
#elif _MACOSX
#define OPS_Export extern "C" __attribute__((visibility("default")))
#else
#define OPS_Export extern "C"
#endif

static int numInno3DPnPJoint = 0;

// class wide matrices
Matrix Inno3DPnPJoint::Transf(30,30);
Matrix Inno3DPnPJoint::Tran(3,3);

void*OPS_Inno3DPnPJoint(void)
{
	// opserr << "OPS_Inno3DPnPJoint: START" << endln;
	
	// // check to be in 3D
	// if (OPS_GetNDM() != 3 || OPS_GetNDF() != 6) {
	// opserr << "WARNING -- model dimensions and/or nodal DOF not compatible with Inno3DPnPJoint element\n";
	// return 0;
    // }
	
	// print out a message about who wrote this element & any copyright info wanted
	if (numInno3DPnPJoint == 0)
	{
		// opserr << "Inno3DJointND element - Written by Cristian V. Miculas @ISISE - UC, Portugal.\n" << endln;
		numInno3DPnPJoint++;
	}

	Element *theInno3DPnPJoint = 0;

	int numRemainingArgs = OPS_GetNumRemainingInputArgs();
  
	if (numRemainingArgs == 0)
	{
		theInno3DPnPJoint = new Inno3DPnPJoint();
		return theInno3DPnPJoint;
	}

	if (numRemainingArgs != 38)
	{
		opserr << "WARNING error insufficient arguments." << endln;
		opserr << "Want: element Inno3DPnPJoint eleTag? Node1? Node2? Node3? Node4? Node5? Spring01? Spring02? ... Spring32?. " << endln;
		numInno3DPnPJoint++;
	}
	
	// opserr << "numRemainingArgs " << numRemainingArgs << endln;
	
	// get the id and end nodes
	int iData[37]; // vector with 38 positions for ele tag, nodes tag and material tags
	
	// variables 1-6 (eleTag + 5xNodes)
	int numData;

	numData = 6;
	if (OPS_GetIntInput(&numData, iData) != 0)
	{
		opserr << "WARNING error invalid element data.\n";
		return 0;
	}
	int eleTag = iData[0];
	
	
	// variables 6-32 (38xSprings)
	
	// - - - - - C O M P O N E N T   1 - - - - -
	numData = 1;
	if (OPS_GetIntInput(&numData, &iData[6]) != 0)
	{
		opserr << "WARNING error reading the stiffness: " << eleTag << endln;
		return 0;
	}

	int matID_1 = iData[6];
	UniaxialMaterial *theK1 = OPS_GetUniaxialMaterial(matID_1);
	if (theK1 == 0)
	{
		opserr << "WARNING stiffness not found for spring: " << matID_1 << endln;
		return 0;
	}


	// - - - - - C O M P O N E N T   2 - - - - -
	numData = 1;
	if (OPS_GetIntInput(&numData, &iData[7]) != 0)
	{
		opserr << "WARNING error reading the stiffness: " << eleTag << endln;
		return 0;
	}

	int matID_2 = iData[7];
	UniaxialMaterial *theK2 = OPS_GetUniaxialMaterial(matID_2);
	if (theK2 == 0)
	{
		opserr << "WARNING stiffness not found for spring: " << matID_2 << endln;
		return 0;
	}


	// - - - - - C O M P O N E N T   3 - - - - -
	numData = 1;
	if (OPS_GetIntInput(&numData, &iData[8]) != 0)
	{
		opserr << "WARNING error reading the stiffness: " << eleTag << endln;
		return 0;
	}

	int matID_3 = iData[8];
	UniaxialMaterial *theK3 = OPS_GetUniaxialMaterial(matID_3);
	if (theK3 == 0)
	{
		opserr << "WARNING stiffness not found for spring: " << matID_3 << endln;
		return 0;
	}


	// - - - - - C O M P O N E N T   4 - - - - -
	numData = 1;
	if (OPS_GetIntInput(&numData, &iData[9]) != 0)
	{
		opserr << "WARNING error reading the stiffness: " << eleTag << endln;
		return 0;
	}

	int matID_4 = iData[9];
	UniaxialMaterial *theK4 = OPS_GetUniaxialMaterial(matID_4);
	if (theK4 == 0)
	{
		opserr << "WARNING stiffness not found for spring: " << matID_4 << endln;
		return 0;
	}


	// - - - - - C O M P O N E N T   5 - - - - -
	numData = 1;
	if (OPS_GetIntInput(&numData, &iData[10]) != 0)
	{
		opserr << "WARNING error reading the stiffness: " << eleTag << endln;
		return 0;
	}

	int matID_5 = iData[10];
	UniaxialMaterial *theK5 = OPS_GetUniaxialMaterial(matID_5);
	if (theK5 == 0)
	{
		opserr << "WARNING stiffness not found for spring: " << matID_5 << endln;
		return 0;
	}


	// - - - - - C O M P O N E N T   6 - - - - -
	numData = 1;
	if (OPS_GetIntInput(&numData, &iData[11]) != 0)
	{
		opserr << "WARNING error reading the stiffness: " << eleTag << endln;
		return 0;
	}

	int matID_6 = iData[11];
	UniaxialMaterial *theK6 = OPS_GetUniaxialMaterial(matID_6);
	if (theK6 == 0)
	{
		opserr << "WARNING stiffness not found for spring: " << matID_6 << endln;
		return 0;
	}


	// - - - - - C O M P O N E N T   7 - - - - -
	numData = 1;
	if (OPS_GetIntInput(&numData, &iData[12]) != 0)
	{
		opserr << "WARNING error reading the stiffness: " << eleTag << endln;
		return 0;
	}

	int matID_7 = iData[12];
	UniaxialMaterial *theK7 = OPS_GetUniaxialMaterial(matID_7);
	if (theK7 == 0)
	{
		opserr << "WARNING stiffness not found for spring: " << matID_7 << endln;
		return 0;
	}


	// - - - - - C O M P O N E N T   8 - - - - -
	numData = 1;
	if (OPS_GetIntInput(&numData, &iData[13]) != 0)
	{
		opserr << "WARNING error reading the stiffness: " << eleTag << endln;
		return 0;
	}

	int matID_8 = iData[13];
	UniaxialMaterial *theK8 = OPS_GetUniaxialMaterial(matID_8);
	if (theK8 == 0)
	{
		opserr << "WARNING stiffness not found for spring: " << matID_8 << endln;
		return 0;
	}


	// - - - - - C O M P O N E N T   9 - - - - -
	numData = 1;
	if (OPS_GetIntInput(&numData, &iData[14]) != 0)
	{
		opserr << "WARNING error reading the stiffness: " << eleTag << endln;
		return 0;
	}

	int matID_9 = iData[14];
	UniaxialMaterial *theK9 = OPS_GetUniaxialMaterial(matID_9);
	if (theK9 == 0)
	{
		opserr << "WARNING stiffness not found for spring: " << matID_9 << endln;
		return 0;
	}


	// - - - - - C O M P O N E N T   10 - - - - -
	numData = 1;
	if (OPS_GetIntInput(&numData, &iData[15]) != 0)
	{
		opserr << "WARNING error reading the stiffness: " << eleTag << endln;
		return 0;
	}

	int matID_10 = iData[15];
	UniaxialMaterial *theK10 = OPS_GetUniaxialMaterial(matID_10);
	if (theK10 == 0)
	{
		opserr << "WARNING stiffness not found for spring: " << matID_10 << endln;
		return 0;
	}


	// - - - - - C O M P O N E N T   11 - - - - -
	numData = 1;
	if (OPS_GetIntInput(&numData, &iData[16]) != 0)
	{
		opserr << "WARNING error reading the stiffness: " << eleTag << endln;
		return 0;
	}

	int matID_11 = iData[16];
	UniaxialMaterial *theK11 = OPS_GetUniaxialMaterial(matID_11);
	if (theK11 == 0)
	{
		opserr << "WARNING stiffness not found for spring: " << matID_11 << endln;
		return 0;
	}


	// - - - - - C O M P O N E N T   12 - - - - -
	numData = 1;
	if (OPS_GetIntInput(&numData, &iData[17]) != 0)
	{
		opserr << "WARNING error reading the stiffness: " << eleTag << endln;
		return 0;
	}

	int matID_12 = iData[17];
	UniaxialMaterial *theK12 = OPS_GetUniaxialMaterial(matID_12);
	if (theK12 == 0)
	{
		opserr << "WARNING stiffness not found for spring: " << matID_12 << endln;
		return 0;
	}


	// - - - - - C O M P O N E N T   13 - - - - -
	numData = 1;
	if (OPS_GetIntInput(&numData, &iData[18]) != 0)
	{
		opserr << "WARNING error reading the stiffness: " << eleTag << endln;
		return 0;
	}

	int matID_13 = iData[18];
	UniaxialMaterial *theK13 = OPS_GetUniaxialMaterial(matID_13);
	if (theK13 == 0)
	{
		opserr << "WARNING stiffness not found for spring: " << matID_13 << endln;
		return 0;
	}


	// - - - - - C O M P O N E N T   14 - - - - -
	numData = 1;
	if (OPS_GetIntInput(&numData, &iData[19]) != 0)
	{
		opserr << "WARNING error reading the stiffness: " << eleTag << endln;
		return 0;
	}

	int matID_14 = iData[19];
	UniaxialMaterial *theK14 = OPS_GetUniaxialMaterial(matID_14);
	if (theK14 == 0)
	{
		opserr << "WARNING stiffness not found for spring: " << matID_14 << endln;
		return 0;
	}


	// - - - - - C O M P O N E N T   15 - - - - -
	numData = 1;
	if (OPS_GetIntInput(&numData, &iData[20]) != 0)
	{
		opserr << "WARNING error reading the stiffness: " << eleTag << endln;
		return 0;
	}

	int matID_15 = iData[20];
	UniaxialMaterial *theK15 = OPS_GetUniaxialMaterial(matID_15);
	if (theK15 == 0)
	{
		opserr << "WARNING stiffness not found for spring: " << matID_15 << endln;
		return 0;
	}


	// - - - - - C O M P O N E N T   16 - - - - -
	numData = 1;
	if (OPS_GetIntInput(&numData, &iData[21]) != 0)
	{
		opserr << "WARNING error reading the stiffness: " << eleTag << endln;
		return 0;
	}

	int matID_16 = iData[21];
	UniaxialMaterial *theK16 = OPS_GetUniaxialMaterial(matID_16);
	if (theK16 == 0)
	{
		opserr << "WARNING stiffness not found for spring: " << matID_16 << endln;
		return 0;
	}


	// - - - - - C O M P O N E N T   17 - - - - -
	numData = 1;
	if (OPS_GetIntInput(&numData, &iData[22]) != 0)
	{
		opserr << "WARNING error reading the stiffness: " << eleTag << endln;
		return 0;
	}

	int matID_17 = iData[22];
	UniaxialMaterial *theK17 = OPS_GetUniaxialMaterial(matID_17);
	if (theK17 == 0)
	{
		opserr << "WARNING stiffness not found for spring: " << matID_17 << endln;
		return 0;
	}


	// - - - - - C O M P O N E N T   18 - - - - -
	numData = 1;
	if (OPS_GetIntInput(&numData, &iData[23]) != 0)
	{
		opserr << "WARNING error reading the stiffness: " << eleTag << endln;
		return 0;
	}

	int matID_18 = iData[23];
	UniaxialMaterial *theK18 = OPS_GetUniaxialMaterial(matID_18);
	if (theK18 == 0)
	{
		opserr << "WARNING stiffness not found for spring: " << matID_18 << endln;
		return 0;
	}


	// - - - - - C O M P O N E N T   19 - - - - -
	numData = 1;
	if (OPS_GetIntInput(&numData, &iData[24]) != 0)
	{
		opserr << "WARNING error reading the stiffness: " << eleTag << endln;
		return 0;
	}

	int matID_19 = iData[24];
	UniaxialMaterial *theK19 = OPS_GetUniaxialMaterial(matID_19);
	if (theK19 == 0)
	{
		opserr << "WARNING stiffness not found for spring: " << matID_19 << endln;
		return 0;
	}


	// - - - - - C O M P O N E N T   20 - - - - -
	numData = 1;
	if (OPS_GetIntInput(&numData, &iData[25]) != 0)
	{
		opserr << "WARNING error reading the stiffness: " << eleTag << endln;
		return 0;
	}

	int matID_20 = iData[25];
	UniaxialMaterial *theK20 = OPS_GetUniaxialMaterial(matID_20);
	if (theK20 == 0)
	{
		opserr << "WARNING stiffness not found for spring: " << matID_20 << endln;
		return 0;
	}


	// - - - - - C O M P O N E N T   21 - - - - -
	numData = 1;
	if (OPS_GetIntInput(&numData, &iData[26]) != 0)
	{
		opserr << "WARNING error reading the stiffness: " << eleTag << endln;
		return 0;
	}

	int matID_21 = iData[26];
	UniaxialMaterial *theK21 = OPS_GetUniaxialMaterial(matID_21);
	if (theK21 == 0)
	{
		opserr << "WARNING stiffness not found for spring: " << matID_21 << endln;
		return 0;
	}


	// - - - - - C O M P O N E N T   22 - - - - -
	numData = 1;
	if (OPS_GetIntInput(&numData, &iData[27]) != 0)
	{
		opserr << "WARNING error reading the stiffness: " << eleTag << endln;
		return 0;
	}

	int matID_22 = iData[27];
	UniaxialMaterial *theK22 = OPS_GetUniaxialMaterial(matID_22);
	if (theK22 == 0)
	{
		opserr << "WARNING stiffness not found for spring: " << matID_22 << endln;
		return 0;
	}


	// - - - - - C O M P O N E N T   23 - - - - -
	numData = 1;
	if (OPS_GetIntInput(&numData, &iData[28]) != 0)
	{
		opserr << "WARNING error reading the stiffness: " << eleTag << endln;
		return 0;
	}

	int matID_23 = iData[28];
	UniaxialMaterial *theK23 = OPS_GetUniaxialMaterial(matID_23);
	if (theK23 == 0)
	{
		opserr << "WARNING stiffness not found for spring: " << matID_23 << endln;
		return 0;
	}


	// - - - - - C O M P O N E N T   24 - - - - -
	numData = 1;
	if (OPS_GetIntInput(&numData, &iData[29]) != 0)
	{
		opserr << "WARNING error reading the stiffness: " << eleTag << endln;
		return 0;
	}

	int matID_24 = iData[29];
	UniaxialMaterial *theK24 = OPS_GetUniaxialMaterial(matID_24);
	if (theK24 == 0)
	{
		opserr << "WARNING stiffness not found for spring: " << matID_24 << endln;
		return 0;
	}


	// - - - - - C O M P O N E N T   25 - - - - -
	numData = 1;
	if (OPS_GetIntInput(&numData, &iData[30]) != 0)
	{
		opserr << "WARNING error reading the stiffness: " << eleTag << endln;
		return 0;
	}

	int matID_25 = iData[30];
	UniaxialMaterial *theK25 = OPS_GetUniaxialMaterial(matID_25);
	if (theK25 == 0)
	{
		opserr << "WARNING stiffness not found for spring: " << matID_25 << endln;
		return 0;
	}


	// - - - - - C O M P O N E N T   26 - - - - -
	numData = 1;
	if (OPS_GetIntInput(&numData, &iData[31]) != 0)
	{
		opserr << "WARNING error reading the stiffness: " << eleTag << endln;
		return 0;
	}

	int matID_26 = iData[31];
	UniaxialMaterial *theK26 = OPS_GetUniaxialMaterial(matID_26);
	if (theK26 == 0)
	{
		opserr << "WARNING stiffness not found for spring: " << matID_26 << endln;
		return 0;
	}


	// - - - - - C O M P O N E N T   27 - - - - -
	numData = 1;
	if (OPS_GetIntInput(&numData, &iData[32]) != 0)
	{
		opserr << "WARNING error reading the stiffness: " << eleTag << endln;
		return 0;
	}

	int matID_27 = iData[32];
	UniaxialMaterial *theK27 = OPS_GetUniaxialMaterial(matID_27);
	if (theK27 == 0)
	{
		opserr << "WARNING stiffness not found for spring: " << matID_27 << endln;
		return 0;
	}


	// - - - - - C O M P O N E N T   28 - - - - -
	numData = 1;
	if (OPS_GetIntInput(&numData, &iData[33]) != 0)
	{
		opserr << "WARNING error reading the stiffness: " << eleTag << endln;
		return 0;
	}

	int matID_28 = iData[33];
	UniaxialMaterial *theK28 = OPS_GetUniaxialMaterial(matID_28);
	if (theK28 == 0)
	{
		opserr << "WARNING stiffness not found for spring: " << matID_28 << endln;
		return 0;
	}


	// - - - - - C O M P O N E N T   29 - - - - -
	numData = 1;
	if (OPS_GetIntInput(&numData, &iData[34]) != 0)
	{
		opserr << "WARNING error reading the stiffness: " << eleTag << endln;
		return 0;
	}

	int matID_29 = iData[34];
	UniaxialMaterial *theK29 = OPS_GetUniaxialMaterial(matID_29);
	if (theK29 == 0)
	{
		opserr << "WARNING stiffness not found for spring: " << matID_29 << endln;
		return 0;
	}


	// - - - - - C O M P O N E N T   30 - - - - -
	numData = 1;
	if (OPS_GetIntInput(&numData, &iData[35]) != 0)
	{
		opserr << "WARNING error reading the stiffness: " << eleTag << endln;
		return 0;
	}

	int matID_30 = iData[35];
	UniaxialMaterial *theK30 = OPS_GetUniaxialMaterial(matID_30);
	if (theK30 == 0)
	{
		opserr << "WARNING stiffness not found for spring: " << matID_30 << endln;
		return 0;
	}


	// - - - - - C O M P O N E N T   31 - - - - -
	numData = 1;
	if (OPS_GetIntInput(&numData, &iData[36]) != 0)
	{
		opserr << "WARNING error reading the stiffness: " << eleTag << endln;
		return 0;
	}

	int matID_31 = iData[36];
	UniaxialMaterial *theK31 = OPS_GetUniaxialMaterial(matID_31);
	if (theK31 == 0)
	{
		opserr << "WARNING stiffness not found for spring: " << matID_31 << endln;
		return 0;
	}


	// - - - - - C O M P O N E N T   32 - - - - -
	numData = 1;
	if (OPS_GetIntInput(&numData, &iData[37]) != 0)
	{
		opserr << "WARNING error reading the stiffness: " << eleTag << endln;
		return 0;
	}

	int matID_32 = iData[37];
	UniaxialMaterial *theK32 = OPS_GetUniaxialMaterial(matID_32);
	if (theK32 == 0)
	{
		opserr << "WARNING stiffness not found for spring: " << matID_32 << endln;
		return 0;
	}


	theInno3DPnPJoint =  new Inno3DPnPJoint(eleTag, iData[1], iData[2], iData[3], iData[4], iData[5], *theK1, *theK2, *theK3, *theK4, *theK5, *theK6, *theK7, *theK8, *theK9, *theK10, *theK11, *theK12, *theK13, *theK14, *theK15, *theK16, *theK17, *theK18, *theK19, *theK20, *theK21, *theK22, *theK23, *theK24, *theK25, *theK26, *theK27, *theK28, *theK29, *theK30, *theK31, *theK32);



	if (theInno3DPnPJoint == 0)
	{
		delete theK1;
		delete theK2;
		delete theK3;
		delete theK4;
		delete theK5;
		delete theK6;
		delete theK7;
		delete theK8;
		delete theK9;
		delete theK10;
		delete theK11;
		delete theK12;
		delete theK13;
		delete theK14;
		delete theK15;
		delete theK16;
		delete theK17;
		delete theK18;
		delete theK19;
		delete theK20;
		delete theK21;
		delete theK22;
		delete theK23;
		delete theK24;
		delete theK25;
		delete theK26;
		delete theK27;
		delete theK28;
		delete theK29;
		delete theK30;
		delete theK31;
		delete theK32;
		return 0;
	}

	// opserr << "OPS_Inno3DPnPJoint: END" << endln;
	return theInno3DPnPJoint;	
}

// full constructors:
Inno3DPnPJoint::Inno3DPnPJoint(int tag, int Nd1, int Nd2, int Nd3, int Nd4, int Nd5,
				    UniaxialMaterial& theMat1,
				    UniaxialMaterial& theMat2,
				    UniaxialMaterial& theMat3,
				    UniaxialMaterial& theMat4,
				    UniaxialMaterial& theMat5,
				    UniaxialMaterial& theMat6,
				    UniaxialMaterial& theMat7,
				    UniaxialMaterial& theMat8,
				    UniaxialMaterial& theMat9,
				    UniaxialMaterial& theMat10,
				    UniaxialMaterial& theMat11,
					UniaxialMaterial& theMat12,
					UniaxialMaterial& theMat13,
					UniaxialMaterial& theMat14,
					UniaxialMaterial& theMat15,
					UniaxialMaterial& theMat16,
					UniaxialMaterial& theMat17,
					UniaxialMaterial& theMat18,
					UniaxialMaterial& theMat19,
					UniaxialMaterial& theMat20,
					UniaxialMaterial& theMat21,
					UniaxialMaterial& theMat22,
					UniaxialMaterial& theMat23,
					UniaxialMaterial& theMat24,
					UniaxialMaterial& theMat25,
					UniaxialMaterial& theMat26,
					UniaxialMaterial& theMat27,
					UniaxialMaterial& theMat28,
					UniaxialMaterial& theMat29,
					UniaxialMaterial& theMat30,
					UniaxialMaterial& theMat31,
					UniaxialMaterial& theMat32):
	Element(tag, ELE_TAG_Inno3DPnPJoint), ExternalNodes(5),
	nodeDbTag(0), dofDbTag(0), dcX (0.0), dcZ (0.0),
	Uecommit(30), UeIntcommit(4), UeprCommit(30), UeprIntCommit(4),
	matA(32,34), dg_df(4,32), dDef_du(32,4), K(30,30), R(34)
{
	// opserr << "Inno3DPnPJoint full constructor: START" << endln;
	
	// ensure the connectedExternalNode ID is of correct size & set values
	if (ExternalNodes.Size() != 5)
    opserr << "ERROR: Inno3DPnPJoint::Inno3DPnPJoint() " << tag << "failed to create an ID of size 5. " << endln;

	ExternalNodes(0) = Nd1;
	ExternalNodes(1) = Nd2;
	ExternalNodes(2) = Nd3;
	ExternalNodes(3) = Nd4;
	ExternalNodes(4) = Nd5;
	
	// opserr << "node 0: " << ExternalNodes(0) << endln;
	// opserr << "node 1: " << ExternalNodes(1) << endln;
	// opserr << "node 2: " << ExternalNodes(2) << endln;
	// opserr << "node 3: " << ExternalNodes(3) << endln;
	// opserr << "node 4: " << ExternalNodes(4) << endln;
	
	// opserr << "nodeDbTag : " << nodeDbTag << endln;
	// opserr << "dofDbTag : " << dofDbTag << endln;
	// opserr << "ExternalNodes : " << ExternalNodes << endln;

	nodePtr[0] = 0;
	nodePtr[1] = 0;
	nodePtr[2] = 0;
	nodePtr[3] = 0;
	nodePtr[4] = 0;

	MaterialPtr = new UniaxialMaterial*[32];
	for (int x = 0; x <32; x++)
	{
		MaterialPtr[x] = 0;
	}

	// get a copy of the material and check we obtained a valid copy
	MaterialPtr[0] = theMat1.getCopy();
	if (!MaterialPtr[0])
	{
		opserr << "ERROR: Inno3DPnPJoint::Constructor failed to get a copy of material 1. " << endln;
		exit(-1);
	}

	MaterialPtr[1] = theMat2.getCopy();
	if (!MaterialPtr[1])
	{
		opserr << "ERROR: Inno3DPnPJoint::Constructor failed to get a copy of material 2. " << endln;
		exit(-1);
	}

	MaterialPtr[2] = theMat3.getCopy();
	if (!MaterialPtr[2])
	{
		opserr << "ERROR: Inno3DPnPJoint::Constructor failed to get a copy of material 3. " << endln;
		exit(-1);
	}

	MaterialPtr[3] = theMat4.getCopy();
	if (!MaterialPtr[3])
	{
		opserr << "ERROR: Inno3DPnPJoint::Constructor failed to get a copy of material 4. " << endln;
		exit(-1);
	}

	MaterialPtr[4] = theMat5.getCopy();
	if (!MaterialPtr[4])
	{
		opserr << "ERROR: Inno3DPnPJoint::Constructor failed to get a copy of material 5. " << endln;
		exit(-1);
	}

	MaterialPtr[5] = theMat6.getCopy();
	if (!MaterialPtr[5])
	{
		opserr << "ERROR: Inno3DPnPJoint::Constructor failed to get a copy of material 6. " << endln;
		exit(-1);
	}

	MaterialPtr[6] = theMat7.getCopy();
	if (!MaterialPtr[6])
	{
		opserr << "ERROR: Inno3DPnPJoint::Constructor failed to get a copy of material 7. " << endln;
		exit(-1);
	}

	MaterialPtr[7] = theMat8.getCopy();
	if (!MaterialPtr[7])
	{
		opserr << "ERROR: Inno3DPnPJoint::Constructor failed to get a copy of material 8. " << endln;
		exit(-1);
	}

	MaterialPtr[8] = theMat9.getCopy();
	if (!MaterialPtr[8])
	{
		opserr << "ERROR: Inno3DPnPJoint::Constructor failed to get a copy of material 9. " << endln;
		exit(-1);
	}

	MaterialPtr[9] = theMat10.getCopy();
	if (!MaterialPtr[9])
	{
		opserr << "ERROR: Inno3DPnPJoint::Constructor failed to get a copy of material 10. " << endln;
		exit(-1);
	}

	MaterialPtr[10] = theMat11.getCopy();
	if (!MaterialPtr[10])
	{
		opserr << "ERROR: Inno3DPnPJoint::Constructor failed to get a copy of material 11. " << endln;
		exit(-1);
	}

	MaterialPtr[11] = theMat12.getCopy();
	if (!MaterialPtr[11])
	{
		opserr << "ERROR: Inno3DPnPJoint::Constructor failed to get a copy of material 12. " << endln;
		exit(-1);
	}

	MaterialPtr[12] = theMat13.getCopy();
	if (!MaterialPtr[12])
	{
		opserr << "ERROR: Inno3DPnPJoint::Constructor failed to get a copy of material 13. " << endln;
		exit(-1);
	}

	MaterialPtr[13] = theMat14.getCopy();
	if (!MaterialPtr[13])
	{
		opserr << "ERROR: Inno3DPnPJoint::Constructor failed to get a copy of material 14. " << endln;
		exit(-1);
	}

	MaterialPtr[14] = theMat15.getCopy();
	if (!MaterialPtr[14])
	{
		opserr << "ERROR: Inno3DPnPJoint::Constructor failed to get a copy of material 15. " << endln;
		exit(-1);
	}

	MaterialPtr[15] = theMat16.getCopy();
	if (!MaterialPtr[15])
	{
		opserr << "ERROR: Inno3DPnPJoint::Constructor failed to get a copy of material 16. " << endln;
		exit(-1);
	}

	MaterialPtr[16] = theMat17.getCopy();
	if (!MaterialPtr[16])
	{
		opserr << "ERROR: Inno3DPnPJoint::Constructor failed to get a copy of material 17. " << endln;
		exit(-1);
	}

	MaterialPtr[17] = theMat18.getCopy();
	if (!MaterialPtr[17])
	{
		opserr << "ERROR: Inno3DPnPJoint::Constructor failed to get a copy of material 18. " << endln;
		exit(-1);
	}

	MaterialPtr[18] = theMat19.getCopy();
	if (!MaterialPtr[18])
	{
		opserr << "ERROR: Inno3DPnPJoint::Constructor failed to get a copy of material 19. " << endln;
		exit(-1);
	}

	MaterialPtr[19] = theMat20.getCopy();
	if (!MaterialPtr[19])
	{
		opserr << "ERROR: Inno3DPnPJoint::Constructor failed to get a copy of material 20. " << endln;
		exit(-1);
	}

	MaterialPtr[20] = theMat21.getCopy();
	if (!MaterialPtr[20])
	{
		opserr << "ERROR: Inno3DPnPJoint::Constructor failed to get a copy of material 21. " << endln;
		exit(-1);
	}

	MaterialPtr[21] = theMat22.getCopy();
	if (!MaterialPtr[21])
	{
		opserr << "ERROR: Inno3DPnPJoint::Constructor failed to get a copy of material 22. " << endln;
		exit(-1);
	}

	MaterialPtr[22] = theMat23.getCopy();
	if (!MaterialPtr[22])
	{
		opserr << "ERROR: Inno3DPnPJoint::Constructor failed to get a copy of material 23. " << endln;
		exit(-1);
	}

	MaterialPtr[23] = theMat24.getCopy();
	if (!MaterialPtr[23])
	{
		opserr << "ERROR: Inno3DPnPJoint::Constructor failed to get a copy of material 24. " << endln;
		exit(-1);
	}

	MaterialPtr[24] = theMat25.getCopy();
	if (!MaterialPtr[24])
	{
		opserr << "ERROR: Inno3DPnPJoint::Constructor failed to get a copy of material 25. " << endln;
		exit(-1);
	}

	MaterialPtr[25] = theMat26.getCopy();
	if (!MaterialPtr[25])
	{
		opserr << "ERROR: Inno3DPnPJoint::Constructor failed to get a copy of material 26. " << endln;
		exit(-1);
	}

	MaterialPtr[26] = theMat27.getCopy();
	if (!MaterialPtr[26])
	{
		opserr << "ERROR: Inno3DPnPJoint::Constructor failed to get a copy of material 27. " << endln;
		exit(-1);
	}

	MaterialPtr[27] = theMat28.getCopy();
	if (!MaterialPtr[27])
	{
		opserr << "ERROR: Inno3DPnPJoint::Constructor failed to get a copy of material 28. " << endln;
		exit(-1);
	}

	MaterialPtr[28] = theMat29.getCopy();
	if (!MaterialPtr[28])
	{
		opserr << "ERROR: Inno3DPnPJoint::Constructor failed to get a copy of material 29. " << endln;
		exit(-1);
	}

	MaterialPtr[29] = theMat30.getCopy();
	if (!MaterialPtr[29])
	{
		opserr << "ERROR: Inno3DPnPJoint::Constructor failed to get a copy of material 30. " << endln;
		exit(-1);
	}

	MaterialPtr[30] = theMat31.getCopy();
	if (!MaterialPtr[30])
	{
		opserr << "ERROR: Inno3DPnPJoint::Constructor failed to get a copy of material 31. " << endln;
		exit(-1);
	}

	MaterialPtr[31] = theMat32.getCopy();
	if (!MaterialPtr[31])
	{
		opserr << "ERROR: Inno3DPnPJoint::Constructor failed to get a copy of material 32. " << endln;
		exit(-1);
	}

	Uecommit.Zero();
	UeIntcommit.Zero();
	UeprCommit.Zero();
	UeprIntCommit.Zero();

	matA.Zero(); 
	dg_df.Zero();
	dDef_du.Zero();

	K.Zero();
	R.Zero();

	
	// opserr << "Inno3DPnPJoint full constructor: END" << endln;
}


// default constructor:
Inno3DPnPJoint::Inno3DPnPJoint():
	Element(0, ELE_TAG_Inno3DPnPJoint), ExternalNodes(5),
	nodeDbTag(0), dofDbTag(0), dcX (0.0), dcZ (0.0),
	Uecommit(30), UeIntcommit(4), UeprCommit(30), UeprIntCommit(4),
	matA(32,34), dg_df(4,32), dDef_du(32,4), K(30,30), R(34)
{
	// opserr << "Inno3DPnPJoint default constructor: START" << endln;
	
	ExternalNodes(0) = 0;
	ExternalNodes(1) = 0;
	ExternalNodes(2) = 0;
	ExternalNodes(3) = 0;
	ExternalNodes(4) = 0;

	// opserr << "node 1: " << ExternalNodes(0) << endln;
	// opserr << "node 2: " << ExternalNodes(1) << endln;
	// opserr << "node 3: " << ExternalNodes(2) << endln;
	// opserr << "node 4: " << ExternalNodes(3) << endln;
	// opserr << "node 5: " << ExternalNodes(4) << endln;
	
	nodePtr[0] = 0; 
	nodePtr[1] = 0;
	nodePtr[2] = 0;
	nodePtr[3] = 0;
	nodePtr[4] = 0;
	
	// does nothing (invoked by FEM_ObjectBroker)	
	for (int i = 0; i <32; i++)
	{
		MaterialPtr[i] = 0;
	}
	
	// opserr << "Inno3DPnPJoint default constructor: END" << endln;
}


// destructor - provided to clean up any memory
Inno3DPnPJoint::~Inno3DPnPJoint()
{
	// opserr << "Destructor: START" << endln;
	
	for (int i =0; i<32; i++)
	{
	    if (MaterialPtr[i] != 0)
		{
			delete MaterialPtr[i];
		}
	}

	if (MaterialPtr)
	{
		delete [] MaterialPtr;
	}

	// opserr << "Destructor: END" << endln;
}


// public methods
int Inno3DPnPJoint::getNumExternalNodes(void) const
{
	// opserr << "getNumExternalNodes: START/END" << endln;
	return 5;
}


const ID & Inno3DPnPJoint::getExternalNodes(void) 
{
	// opserr << "getExternalNodes: START/END" << endln;
	return this->ExternalNodes;
}


Node **Inno3DPnPJoint::getNodePtrs(void) 
{
	// opserr << "getNodePtrs: START/END" << endln;
	return nodePtr;
}


int Inno3DPnPJoint::getNumDOF(void) 
{
	// opserr << "getNumDOF: START/END" << endln;
    return 30;
}


void Inno3DPnPJoint::setDomain(Domain *theDomain)
{
	// opserr << "setDomain: START" << endln;
	
	// CHECK 0: Domain is not null - invoked when object removed from a domain
    if (theDomain == 0)
	{
		opserr << "ERROR: Inno3DPnPJoint::setDomain -- Domain is null. " << endln;
		return;
	}
	
	// set node pointers
	for (int i = 0; i < 5; i++)
	{
		nodePtr[i] = theDomain -> getNode(ExternalNodes(i));
		if (nodePtr[i] == 0)
		{
			opserr << "ERROR: Inno3DPnPJoint::setDomain. Node pointer is NULL. Node " << this->ExternalNodes(i) << " does not exit in the domain." << endln;
			exit(-1); // don't go any further
		}
    }
	
	// call the base class method
	this->DomainComponent::setDomain(theDomain);

	// CHECK 1: NUMBER OF DOFs
	int dofNd1 = nodePtr[0]->getNumberDOF();
	int dofNd2 = nodePtr[1]->getNumberDOF();
	int dofNd3 = nodePtr[2]->getNumberDOF();
	int dofNd4 = nodePtr[3]->getNumberDOF();
	int dofNd5 = nodePtr[4]->getNumberDOF();
	
	// check for proper nodal degrees of freedom
	if ((dofNd1 != 6) || (dofNd2 != 6) || (dofNd3 != 6) || (dofNd4 != 6) || (dofNd5 != 6))
	{
		opserr << "ERROR: Inno3DPnPJoint::setDomain -- number of DOF associated with the nodes is incorrect." << endln;
		exit(-1); // don't go any further
	}
	
	
	// get coordinates for other checks
	const Vector &node1Crd = nodePtr[0]->getCrds();
	const Vector &node2Crd = nodePtr[1]->getCrds();
	const Vector &node3Crd = nodePtr[2]->getCrds();
	const Vector &node4Crd = nodePtr[3]->getCrds();
	const Vector &node5Crd = nodePtr[4]->getCrds();
	
	Vector Node1(node1Crd);
	Vector Node2(node2Crd);
	Vector Node3(node3Crd);
	Vector Node4(node4Crd);
	Vector Node5(node5Crd);
	
	
	// CHECK 2: COPLANARITY OF NODES
	// create plane with Node 1, Node 2 and Node 3
	double a1 = Node2(0) - Node1(0);
	double b1 = Node2(1) - Node1(1);
	double c1 = Node2(2) - Node1(2);

	double a2 = Node3(0) - Node1(0);
	double b2 = Node3(1) - Node1(1);
	double c2 = Node3(2) - Node1(2);

	double a = b1*c2 - b2*c1;
	double b = a2*c1 - a1*c2;
	double c = a1*b2 - b1*a2;

	// equation of plane is: a*x + b*y + c*z = 0 #
	double d = (-a*Node1(0) - b*Node1(1) - c*Node1(2));

	double errorDomain = 1e-2; // depending on node coord, value 0 cannot be achieved.
	
	// checking if the Node 4 belongs to plane of Node 1, Node 2 and Node 3
	if(fabs(a * Node4(0) + b * Node4(1) + c * Node4(2) + d ) >= errorDomain)
	{
		opserr << "ERROR: Inno3DPnPJoint::setDomain -- Node 4 does NOT belong to plane created by Node 1, Node 2 and Node 3. Check node coordinate definition. \n" << endln;
		exit(-1); // don't go any further
	}	
	
	// checking if the Node 5 belongs to plane of Node 1, Node 2 and Node 3
	if(fabs(a * Node5(0) + b * Node5(1) + c * Node5(2) + d) >= errorDomain)
	{
		opserr << "ERROR: Inno3DPnPJoint::setDomain -- Node 5 does NOT belong to plane created by Node 1, Node 2 and Node 3. Check node coordinate definition. \n" << endln;
		exit(-1); // don't go any further
	}
	
	
	// CHECK 3: PERPENDICULARITY OF NODES' AXIS
	Vector vector13_G = Node1 - Node3;
	Vector vector24_G = Node2 - Node4;
	
	double dotProduct = 0;
	for (int i = 0; i < Node1.Size(); i++)
	{
		dotProduct = dotProduct + vector13_G(i) * vector24_G(i);
    }
	
	if(fabs(dotProduct) >= errorDomain)
	{
		opserr << "ERROR: Inno3DPnPJoint::setDomain -- vector of Node 1 & Node 3 not perpendicular to vector of Node 2 & Node 4. Check node coordinate definition. \n" << endln;
		exit(-1); // don't go any further
	}
	
	
	// CHECK 4: COLINEARITY OF NODES' AXIS
	// colinearity Node 1, Node 5 and Node 3
	Vector vector51_G = Node5 - Node1;
	Vector vector53_G = Node5 - Node3;
	
	Vector crossProd_153(3);
	crossProd_153.Zero();
	
	crossProd_153[0] = vector51_G[1] * vector53_G[2] - vector51_G[2] * vector53_G[1];
	crossProd_153[1] = vector51_G[2] * vector53_G[0] - vector51_G[0] * vector53_G[2];
	crossProd_153[2] = vector51_G[0] * vector53_G[1] - vector51_G[1] * vector53_G[0];

	if(crossProd_153 != 0)
	{
		opserr << "ERROR: Inno3DPnPJoint::setDomain -- Node 1, Node 5 and Node 3 are not colinear. Check node coordinate definition. \n" << endln;
		exit(-1); // don't go any further
	}

	// colinearity Node 2, Node 5 and Node 4
	Vector vector52_G = Node5 - Node2;
	Vector vector54_G = Node5 - Node4;
	
	Vector crossProd_254(3);
	crossProd_254.Zero();
	
	crossProd_254[0] = vector52_G[1] * vector54_G[2] - vector52_G[2] * vector54_G[1];
	crossProd_254[1] = vector52_G[2] * vector54_G[0] - vector52_G[0] * vector54_G[2];
	crossProd_254[2] = vector52_G[0] * vector54_G[1] - vector52_G[1] * vector54_G[0];

	if(crossProd_254 != 0)
	{
		opserr << "ERROR: Inno3DPnPJoint::setDomain -- Node 2, Node 5 and Node 4 are not colinear. Check node coordinate definition. \n" << endln;
		exit(-1); // don't go any further
	}	

	// determine the element width and height. It should be noted that the element width is determined by nodal coordinates 2 and 4 whereas the element height is determined by the nodal coordinates 1 and 3
	// general order anticlockwise arrangement, but can be clockwise too.	
	
	// CHECK 5: (CENTRALITY) POSITION MIDDLE COLUMN'S FACE and SIZE
	double norm_51 = fabs(vector51_G.Norm());
	double norm_53 = fabs(vector53_G.Norm());
	
	double norm_52 = fabs(vector52_G.Norm());
	double norm_54 = fabs(vector54_G.Norm());
	
	// position
	if( (norm_51 != norm_53) || (norm_52 != norm_54) )
	{
		opserr << "ERROR: Inno3DPnPJoint::setDomain -- nodes are not located at the center of the column's face. Check node coordinate definition. \n" << endln;
		exit(-1); // don't go any further
	}
	// size (if < 1e-3, too small)
	else if( (norm_51 <= 1e-3) || (norm_52 <= 1e-3) )
	{
		opserr << "ERROR: Inno3DPnPJoint::setDomain -- length or width <= 1e-3, division by zero occurs. Increase joint size." << endln;
		exit(-1); // don't go any further
	}
	
	// // check cross-section shape (could be commented out)
	// if (norm_51 == norm_52)
	// {
		// opserr << "msg: Inno3DPnPJoint::setDomain -- column type: square." << endln;
	// }
	// else
	// {
		// opserr << "msg: Inno3DPnPJoint::setDomain -- column type: rectangle." << endln;
	// }
	
	// column width and height (in LOCAL plane X-Z)
	dcZ = fabs(vector13_G.Norm());
	dcX = fabs(vector24_G.Norm());
	
	getmatA();

	getdDef_du();
	getdg_df();
	
	formTransfMat();
	
	// opserr << "setDomain: END" << endln;
}   


int Inno3DPnPJoint::commitState(void)
{
	// opserr << "commitState: START" << endln;
	
	// store committed external nodal displacements
	Uecommit = UeprCommit;
	
	// store committed internal nodal displacements
	UeIntcommit = UeprIntCommit;

	// store material history data.
	int mcs = 0;
		for (int j=0; j<32; j++)
		{
			if (MaterialPtr[j] != 0)
			{
				mcs = MaterialPtr[j]->commitState();
			}
			
			if (mcs != 0)
			{
				break;
			}
		}

	// opserr << "commitState: END" << endln;
	return mcs;
}


int Inno3DPnPJoint::revertToLastCommit(void)
{
	// opserr << "revertToLastCommit: START" << endln;
	
	int mcs = 0;
	for (int j=0; j<32; j++)
	{
		if (MaterialPtr[j] != 0)
		{
			mcs = MaterialPtr[j]->revertToLastCommit();
		}
		
		if (mcs != 0)
		{
			break;
		}
	}
	
	UeprCommit = Uecommit;
	UeprIntCommit = UeIntcommit;   
	
	this->update();
	
	// opserr << "revertToLastCommit: END" << endln;
	return mcs;
}


int Inno3DPnPJoint::revertToStart(void)
{
	// opserr << "revertToStart: START" << endln;
	
	int mcs = 0;
	for (int j=0; j<32; j++)
	{
		if (MaterialPtr[j] != 0){
			mcs = MaterialPtr[j]->revertToStart();
		}
		
		if (mcs != 0)
		{
			break;
		}
	}
	
	// opserr << "revertToStart: END" << endln;
	return mcs;
}


int Inno3DPnPJoint::update(void)
{   
	// opserr << "update: START" << endln;
	
	Vector Ue(34);
	Ue.Zero();

    // determine committed displacements given trial displacements
	getGlobalDispls(Ue);
	// opserr << "Global Displs (all) " << Ue << endln;
			
	// update displacements for the external nodes
	UeprCommit.Extract(Ue,0,1.0);
	// opserr << "Global Displs (ext) " << UeprCommit << endln;

	// update displacement for the internal nodes
	UeprIntCommit.Extract(Ue,30,1.0);
	// opserr << "Global Displs (int) " << UeprIntCommit << endln;

	// opserr << "update: END" << endln;
	return 0;
}


const Matrix & Inno3DPnPJoint::getTangentStiff(void)
{
	// opserr << "getTangentStiff: START" << endln;
	
	Vector kSpring(32);
	kSpring.Zero();
	for (int i = 0; i < 32; i++)
	{
		kSpring[i] = 0;
		if (MaterialPtr[i] != NULL)
		{
			kSpring[i] = MaterialPtr[i]->getTangent();
		}
	}
	
	// opserr << "getTangentStiff: " << kSpring << endln;
	
	formK(kSpring);
	
	// print K_joint (condensed 36x36)
	// opserr << "getTangentStiff: " << kSpring << endln;
	
	// opserr << "getTangentStiff: END" << endln;
	return K;
}


const Matrix & Inno3DPnPJoint::getInitialStiff(void)
{
	// opserr << "getInitialStiff: START/END" << endln;
	
	return getTangentStiff();
}


const Vector & Inno3DPnPJoint::getResistingForce(void)
{
	// opserr << "getResistingForce: START/END" << endln;
	// opserr << "getResistingForce: R " << R << endln;
	return R;
}


void Inno3DPnPJoint::getGlobalDispls(Vector &dg) 
{
	// opserr << "getGlobalDispls: START" << endln;
	
	// local variables that will be used in this method
	int converge = 0;
	int linesearch = 0;
	int totalCount = 0;
	int dtConverge = 0;
	int incCount = 0;
	int count = 0;
	int maxTotalCount = 1000;
	int maxCount = 1000;
	double loadStep = 0.0;
	double dLoadStep = 1.0;
	double stepSize = 0.0;
	
	Vector uExtOld(30);   uExtOld.Zero();
	Vector uExt(30);      uExt.Zero();
	Vector duExt(30);     duExt.Zero();
	Vector uIntOld(4);    uIntOld.Zero(); 
	Vector uInt(4);       uInt.Zero();
	Vector duInt(4);      duInt.Zero(); 
	Vector duIntTemp(4);  duIntTemp.Zero();
	Vector intEq(4);      intEq.Zero();
	Vector intEqLast(4);  intEqLast.Zero();
	Vector Uepr(30);      Uepr.Zero();
	Vector UeprInt(4);    UeprInt.Zero();
	
	
	Vector disp1 = nodePtr[0]->getTrialDisp(); 
	Vector disp2 = nodePtr[1]->getTrialDisp();
	Vector disp3 = nodePtr[2]->getTrialDisp();
	Vector disp4 = nodePtr[3]->getTrialDisp();
	Vector disp5 = nodePtr[4]->getTrialDisp();

	// vector containing external (GLOBAL) disp for each dof
	Vector Ut_G(30);
	Ut_G.Zero();
	for (int i = 0; i < 6; i++)
    {
		Ut_G(i)     = disp1(i);
		Ut_G(i+6)   = disp2(i);
		Ut_G(i+12)  = disp3(i);
		Ut_G(i+18)  = disp4(i);
		Ut_G(i+24)  = disp5(i);
    }
	// opserr << "getGlobalDispls Ut_G: " 		<<  Ut_G 		<< endln;

	Vector Ut(30);
	Ut.Zero();
	
	// transformation GLOBAL -> LOCAL
	Ut.addMatrixVector(0.0, Transf, Ut_G, 1.0);
	
	// opserr << "coordinates Ut : " <<  Ut << endln;

	Uepr = Uecommit;
	uExtOld = Uepr;
	uExt = uExtOld;
	
	UeprInt = UeprIntCommit;
	uIntOld = UeprInt;
	uInt = uIntOld;
	
	duExt = Ut - Uepr;

	double tol = 1e-12;
	
	double tolIntEq = tol;
	double normIntEq = tolIntEq;
	
	double tolIntEqdU = tol;
	double normIntEqdU = tolIntEqdU;
	
	double ctolIntEq = tol;
	double ctolIntEqdU = tol;
	
	double toluInt = (tol>tol*uInt.Norm()) ? tol:tol*uInt.Norm();
	
	//opserr << "getGlobalDispls uInt "		<< uInt;
	//opserr << "getGlobalDispls uInt.Norm " 	<< uInt.Norm() << endln;
	//opserr << "getGlobalDispls toluInt "	<< toluInt;
	
	double normDuInt = toluInt;
	
	Vector u(34);
	u.Zero();

	double engrLast = 0.0;
	double engr = 0.0;

	Vector fSpring(32);  		 fSpring.Zero();
	Vector kSpring(32);  		 kSpring.Zero();
	Matrix dintEq_du(4,4);   	 dintEq_du.Zero();
	Matrix df_dDef(32,32);   	 df_dDef.Zero(); //create a Diag matr of springs
	Matrix tempintEq_du (4,32);  tempintEq_du.Zero();

	while ((loadStep < 1.0) && (totalCount < maxTotalCount))
	{
		count = 0;
		converge = 0;
		dtConverge = 0;
		
		while ((!converge) && (count < maxCount))
		{
			// opserr << "getGlobalDispls: converge: "  <<  converge << endln;
			// opserr << "getGlobalDispls: dLoadStep: " <<  dLoadStep << endln;
			
			if (dLoadStep <= 1e-3)
			{
				dLoadStep = dLoadStep;
			}
			totalCount = totalCount + 1;
			count = count + 1;
			
			// opserr << "getGlobalDispls: totalCount: "  <<  totalCount << endln;
			// opserr << "getGlobalDispls: count: "  <<  count << endln;

			for (int ic = 0; ic < 30; ic++ )
			{
				u(ic) = uExt(ic) + duExt(ic);
			}
			u(30) = uInt(0);
			u(31) = uInt(1);
			u(32) = uInt(2);
			u(33) = uInt(3);
		
			// opserr << "getGlobalDispls uInt " << uInt << endln;

			fSpring.Zero();
			kSpring.Zero();
	
			getMatResponse(u,fSpring,kSpring);
			
			// performs internal equilibrium
			// plane Z-X
			intEq(0) = -fSpring(2)  + fSpring(24) + fSpring(28) + fSpring(29);

			intEq(1) = -fSpring(6)  + fSpring(25) + fSpring(29) + fSpring(31);

			intEq(2) =  fSpring(14) - fSpring(26) - fSpring(30) - fSpring(31);

			intEq(3) =  fSpring(18) - fSpring(27) - fSpring(28) - fSpring(30);

			df_dDef.Zero();
			matDiag(kSpring, df_dDef); // create diag matrix from kSpring
			
			tempintEq_du.Zero();
			dintEq_du.Zero();

			tempintEq_du.addMatrixProduct(0.0,dg_df,df_dDef,1.0);
			dintEq_du.addMatrixProduct(0.0,tempintEq_du,dDef_du,1.0);

			normIntEq = intEq.Norm();
			
			// opserr << "getGlobalDispls normIntEq " 	 << normIntEq << endln; 
			
			normIntEqdU = 0.0;
			for (int jc = 0; jc<4 ; jc++)
			{
				normIntEqdU += intEq(jc)*duInt(jc);
				// opserr << "normIntEqdU 1: " <<  normIntEqdU << endln;
			}
			normIntEqdU = fabs(normIntEqdU);
			// opserr << "normIntEqdU 2: " <<  normIntEqdU << endln;
	
			if (totalCount == 1)
			{
				tolIntEq = (tol>tol*normIntEq) ? tol:tol*normIntEq;
				tolIntEqdU = tol;
			}
			else if (totalCount == 4)
			{
				tolIntEqdU = (tol>tol*normIntEqdU) ? tol:tol*normIntEqdU;
			}
			ctolIntEq = (tolIntEq*dLoadStep > tol) ? tolIntEq*dLoadStep:tol;
			ctolIntEqdU = (tolIntEqdU*dLoadStep > tol) ? tolIntEqdU*dLoadStep:tol;


			// check for convergence starts
			if ((normIntEq < ctolIntEq) || ((normIntEqdU < ctolIntEqdU) && (count >1)) || (normDuInt < toluInt) || (dLoadStep < 1e-3))
			{
				converge = 1;
				// opserr << "convergence starts" << endln;
				
				loadStep = loadStep + dLoadStep;
				if (fabs(1.0 - loadStep) < tol)
				{
					loadStep = 1.0;
				}
			}
			else
			{
				// duInt = -dintEq_du/intEq
				dintEq_du.Solve(intEq,duInt);
				duInt *= -1;
	
				normDuInt = duInt.Norm();
				if (!linesearch)
				{
					uInt = uInt + duInt;
				}
				else
				{
					engrLast = 0.0;
					engr = 0.0;
	
					for (int jd = 0; jd<4 ; jd++)
					{
						engrLast += duInt(jd)*intEqLast(jd);
						engr += duInt(jd)*intEq(jd);
					}
	
					if (fabs(engr) > tol*engrLast)
					{
						duIntTemp = duInt;
						duIntTemp *= -1;
						// lineSearch algorithm requirement
						stepSize = getStepSize(engrLast,engr,uExt,duExt,uInt,duIntTemp,tol);
							
						if (fabs(stepSize) > 0.001)
						{
							uInt = uInt + stepSize*duInt;
						}
						else
						{
							uInt = uInt + duInt;
						}
					}
					else
					{
						uInt = uInt + duInt;
					}
					intEqLast = intEq;
				}
			}
    }

	if (!converge && loadStep < 1.0)
	{
		incCount = 0;
		maxCount = 1000;
		if (!linesearch)
		{
			linesearch = 1;
			uInt = uIntOld;
			duInt.Zero();
		}
		else
		{
			uInt = uIntOld;
			duInt.Zero();
			duExt = duExt*0.1;
			dLoadStep = dLoadStep*0.1;
		}
	}
	else if (loadStep < 1.0)
	{
		maxCount = 500;
		incCount = incCount + 1;
		normDuInt = toluInt;
		if ((incCount < maxCount) || dtConverge)
		{
			uExt = uExt + duExt;
			if (loadStep + dLoadStep > 1.0)
			{
				duExt = duExt*(1.0 - loadStep)/dLoadStep;
				dLoadStep = 1.0 - loadStep;
				incCount = 9;
			}
		}
		else
		{
			incCount = 0;
			uExt = uExt + duExt;
			dLoadStep = dLoadStep*10;
			if (loadStep + dLoadStep > 1.0)
			{
				uExt = uExt + duExt*(1.0 - loadStep)/dLoadStep;
				dLoadStep = 1.0 - loadStep;
				incCount = 9;
			}
		}
	}
	
    }

    // determination of stiffness matrix and the residual force vector for the element
    formR(fSpring);
	formK(kSpring);

    dg.Zero();
    // commmited external and internal displacement update
    for (int ig = 0; ig < 30; ig++ )
	{
		if (ig<30)
		{
			dg(ig) = Ut(ig);
			// opserr << "dg(ig) " << ig << " " <<  dg(ig) << endln;
		}               
	}
	dg(30) = uInt(0);
	dg(31) = uInt(1);
	dg(32) = uInt(2);
	dg(33) = uInt(3);

	
	// opserr << "getGlobalDispls: - dg: " << dg << endln;
	// opserr << "getGlobalDispls - fSpring: " << fSpring << endln;
	
	// opserr << "getGlobalDispls: END" << endln;
}


void Inno3DPnPJoint::getMatResponse(Vector U, Vector &fS, Vector &kS)
{
	// opserr << "getMatResponse: START" << endln;
	
    // obtains the material response from the material class
    Vector defSpring(32);
    defSpring.Zero();
    fS.Zero();
    kS.Zero();
	
	// calcualte constitutive matrix
    defSpring.addMatrixVector(0.0, matA, U, 1.0);
	
	// opserr << "defSpring " << defSpring << endln;
	
    for (int j=0; j<32; j++)
    {
		MaterialPtr[j]->setTrialStrain(defSpring(j));
		// opserr <<  MaterialPtr[j]->setTrialStrain(defSpring(j)) << endln;
		
		kS(j) = MaterialPtr[j]->getTangent();
		// opserr << kS(j) << endln;
				
		fS(j) = MaterialPtr[j]->getStress();
		// opserr << fS(j) << endln;
	}

	// opserr << "getMatResponse - fS: " << fS << endln;
	// opserr << "getMatResponse - kS: " << kS << endln; //vector
	// opserr << "getMatResponse - defSpring: " << defSpring << endln;
	
	// opserr << "getMatResponse: END" << endln;
}


void Inno3DPnPJoint::formTransfMat()
{
	// opserr << "formTransfMat: START" << endln;
	
	Transf.Zero();
	Tran.Zero();
	
	// opserr << "Inno3DPnPJoint::formTransfMat Transf: " << Transf << endln;
	// opserr << "Inno3DPnPJoint::formTransfMat Tran: " << Tran << endln;
	
	const Vector &node1Crd = nodePtr[0]->getCrds();
    const Vector &node2Crd = nodePtr[1]->getCrds(); 

	// opserr << "Inno3DPnPJoint::formTransfMat xNorm: " << node1Crd << endln;
	// opserr << "Inno3DPnPJoint::formTransfMat yNorm: " << node2Crd << endln;
	
	Vector Node1(node1Crd);
	Vector Node2(node2Crd);
	
	// opserr << "Inno3DPnPJoint::formTransfMat xNorm: " << Node1 << endln;
	// opserr << "Inno3DPnPJoint::formTransfMat yNorm: " << Node2 << endln;
	
	Vector crossProd_12(3);
	crossProd_12.Zero();
	
	// opserr << "Inno3DPnPJoint::formTransfMat crossProd_12: " << crossProd_12 << endln;

	crossProd_12[0] = Node1[1] * Node2[2] - Node1[2] * Node2[1];
	crossProd_12[1] = Node1[2] * Node2[0] - Node1[0] * Node2[2];
	crossProd_12[2] = Node1[0] * Node2[1] - Node1[1] * Node2[0];
	
	double xNorm = Node2.Norm();
	double yNorm = crossProd_12.Norm();
	double zNorm = Node1.Norm();
	
	// opserr << "Inno3DPnPJoint::formTransfMat xNorm: " << xNorm << endln;
	// opserr << "Inno3DPnPJoint::formTransfMat yNorm: " << yNorm << endln;
	// opserr << "Inno3DPnPJoint::formTransfMat zNorm: " << zNorm << endln;
	
    // check valid x and y vectors, i.e. not parallel and of zero length
    if (xNorm == 0 || yNorm == 0 || zNorm == 0)
	{
      opserr << "ERROR: Inno3DPnPJoint::formTransfMat -- invalid vectors. Cannot compute Transformation Matrix. Check node coordinate definition." << endln;
	  return;
    }
    
    // create transformation matrix of direction cosines
    for (int i=0; i<3; i++ )
	{
		Tran(0,i) = Node2(i)/xNorm;
		Tran(1,i) = crossProd_12(i)/yNorm;
		Tran(2,i) = Node1(i)/zNorm;
	}
	
	// opserr << "Inno3DPnPJoint::formTransfMat -- Tran(3x3): " << Tran << endln;

	Transf.Assemble(Tran,0,0,1.0);  	//  1
	Transf.Assemble(Tran,3,3,1.0); 		//  2
	Transf.Assemble(Tran,6,6,1.0);  	//  3
	Transf.Assemble(Tran,9,9,1.0);  	//  4
	Transf.Assemble(Tran,12,12,1.0);	//  5
	Transf.Assemble(Tran,15,15,1.0);	//  6
	Transf.Assemble(Tran,18,18,1.0);  	//  7
	Transf.Assemble(Tran,21,21,1.0);	//  8
	Transf.Assemble(Tran,24,24,1.0);	//  9
	Transf.Assemble(Tran,27,27,1.0);	// 10
	
	// opserr << "Inno3DPnPJoint::formTransfMat -- Transf(30x30): " << Transf << endln;
	
	// opserr << "formTransfMat: END" << endln;
}


void Inno3DPnPJoint::getdDef_du()
{   
	// opserr << "getdDef_du: START" << endln;
	
	dDef_du.Zero();
	// extract last 4 columns of compatibility matrix (=internal DOFs)
	dDef_du.Extract(matA,0,30,1.0);

	// opserr << "dDef_du " <<  dDef_du << endln;
	
	// opserr << "getdDef_du: END" << endln;
}


void Inno3DPnPJoint::matDiag(Vector k,Matrix &dfd)
{
	// opserr << "matDiag: START" << endln;
	
	// THIS FUNC RETURS THE DIAG MATRIX (2 input)
    dfd.Zero();
	// the values are from the tcl file!
    // takes in a vector and converts it to a diagonal matrix
	for (int ja=0; ja<32; ja++)
	{
		dfd(ja,ja) = k(ja);
	}
	// opserr << "matDiag: " << dfd << endln;	

	// opserr << "matDiag: END" << endln;
}


void Inno3DPnPJoint::formR(Vector f)
{
	// opserr << "formR: START" << endln;
	
    // develops the element residual force vector
	Vector rForceTemp(34);
	rForceTemp.Zero();
	// opserr << "formR: rForceTemp: " << rForceTemp << endln;
	
	// rForceTemp = rForceTemp + (matA')*f
	rForceTemp.addMatrixTransposeVector(0.0,matA,f,1.0);
	
	Vector R_ext(30);
	R_ext.Zero();
	
	Vector R_int(4);
	R_int.Zero();
	
	R_ext.Extract(rForceTemp,0,1.0);
	R_int.Extract(rForceTemp,30,1.0);

	Vector R_ext_TR(30);
	R_ext_TR.Zero();
	
	// R_ext_TR = Transf'*R_ext
	R_ext_TR.addMatrixTransposeVector(0.0, Transf, R_ext, 1.0);
	
	// Vector R_new(34);
	// R_new.Zero();
	
	// R_new.Assemble(R_ext_TR,0,1.0);  	//  1
	// R_new.Assemble(R_int,30,1.0);  	//  1
	
	R.Zero();
	R.Assemble(R_ext_TR,0,1.0);  	//  1
	R.Assemble(R_int,30,1.0);  	//  1
	
	// opserr << "formR: rForceTemp: " 		<< rForceTemp 	<< endln;
	// opserr << "Reactions R_ext: "		<< R_ext 		<< endln;
	// opserr << "Reactions R_int: "		<< R_int 		<< endln;
	// opserr << "Reactions R_ext_TR: "		<< R_ext_TR 	<< endln;
	// opserr << "Reactions R: "			<< R 			<< endln;

	// opserr << "formR: END" << endln;
}


void Inno3DPnPJoint::formK(Vector k)
{
	// opserr << "formK: START" << endln;
	
    // develops the element stiffness matrix
    Matrix kSprDiag(32,32);
    kSprDiag.Zero();
    
	Matrix kRForce(34,34);
    kRForce.Zero();
    
	Matrix kRFT1(4,30);
    kRFT1.Zero();
    
	Matrix kRFT2(4,4);
    kRFT2.Zero();
    
	Matrix kRFT3(30,4);
    kRFT3.Zero();
    
	Matrix I(4,4);
    I.Zero();
    
	Matrix kRSTinv(4,4);
    kRSTinv.Zero();
    
	Matrix kRF(30,30);
    kRF.Zero();
    
	Matrix K2Temp(30,4);
    K2Temp.Zero();
    
	Matrix K2(30,30);
    K2.Zero();

    matDiag(k,kSprDiag);

	// kRForce = matA'*kSprDiag*matA
    kRForce.addMatrixTripleProduct(0.0,matA,kSprDiag,1.0);   

	// opserr << "formK: matA " << matA << endln;	
    
	kRFT2.Extract(kRForce,30,30,1.0);
    kRFT1.Extract(kRForce,30,0,1.0);
    kRFT3.Extract(kRForce,0,30,1.0);
    kRF.Extract(kRForce,0,0,1.0);
	
	// these matrices don't change!!!
	// opserr << " formK: kRForce! "	<< kRForce << endln;
	// opserr << " formK: kRFT1! "	<< kRFT1 << endln;
	// opserr << " formK: kRFT2! "	<< kRFT2 << endln;
	// opserr << " formK: kRFT3! "	<< kRFT3 << endln;
	// opserr << " formK: kRF!"	<< kRF	 << endln;
    
    for (int ic=0; ic<4; ic++)
    {
        I(ic,ic) = 1.0;
    }
	
    kRFT2.Solve(I,kRSTinv);
	// opserr << " formK 1: kRSTinv!"	<< kRSTinv	 << endln;
    
    K2Temp.addMatrixProduct(0.0,kRFT3,kRSTinv,1.0);
	// opserr << " formK 2: K2Temp!"	<< K2Temp	 << endln;

    for(int i = 0; i <30; ++i)
    {   
        for(int j = 0; j < 4; ++j)
        {
            if (fabs(K2Temp(i,j)) < 1e-15)
                K2Temp(i,j) = 0.;
        }
    }
	
	// // opserr << " formK 3: K2Temp!"	<< K2Temp	 << endln;

    K2.addMatrixProduct(0.0,K2Temp,kRFT1,1.0);
	// opserr << " formK 4: K2!"	<< K2	 << endln;

    for(int i1 = 0; i1 < 30; ++i1)
    {   
        for(int j1 = 0; j1 < 30; ++j1)
        {
            if(fabs(K2(i1,j1)) < 1e-15)
			{
                K2(i1,j1) = 0.;
			}
        }
    }
	
	// opserr << " formK 5: K2!"	<< K2	 << endln;
	
	kRF.addMatrix(1.0,K2,-1.0);
	// opserr << " formK 6: kRF!"	<< kRF	 << endln;
	
	Matrix Knew(30,30);
    Knew.Zero();
	
	// K = Transf'*kRF*Transf
	K.addMatrixTripleProduct(0.0,Transf,kRF,1.0);   

	
	// opserr << "formK: K "     << K << endln;
	// opserr << "formK: K_new " << Knew << endln;
	// opserr << "formK: K "     << K << endln;
	// opserr << "formK: K!" << K << endln;
	
	// opserr << "formK: END" << endln;
}


void Inno3DPnPJoint::getdg_df()
{
	// opserr << "getdg_df: START" << endln;
	
	dg_df.Zero();

	dg_df.AssembleTranspose(dDef_du, 0, 0, 1.0);

	// opserr << "getdg_df:" << dg_df << endln;

	// opserr << "getdg_df: END" << endln;
}


void Inno3DPnPJoint::getmatA()
{
	// opserr << "getmatA: START" << endln;

	// the compatibility matrix of the element
    matA.Zero(); 

	// // PLANE X-Z
	matA(0,0) = 1.00;
	matA(0,24) = -1.00;
	matA(0,28) = -dcZ/2.00;
	matA(1,1) = 1.00;
	matA(1,25) = -1.00;
	matA(1,27) = dcZ/2.00;
	matA(2,2) = 1.00;
	matA(2,30) = -1.00;
	matA(3,3) = -1.00;
	matA(3,27) = -1.00;
	matA(4,4) = 1.00;
	matA(4,28) = -1.00;
	matA(5,5) = 1.00;
	matA(5,29) = -1.00;
	matA(6,6) = 1.00;
	matA(6,31) = -1.00;
	matA(7,7) = -1.00;
	matA(7,25) = 1.00;
	matA(7,29) = dcX/2.00;
	matA(8,8) = 1.00;
	matA(8,26) = -1.00;
	matA(8,28) = dcX/2.00;
	matA(9,9) = 1.00;
	matA(9,27) = -1.00;
	matA(10,10) = 1.00;
	matA(10,28) = -1.00;
	matA(11,11) = 1.00;
	matA(11,29) = -1.00;
	matA(12,12) = -1.00;
	matA(12,24) = 1.00;
	matA(12,28) = -dcZ/2.00;
	matA(13,13) = -1.00;
	matA(13,25) = 1.00;
	matA(13,27) = dcZ/2.00;
	matA(14,14) = -1.00;
	matA(14,32) = 1.00;
	matA(15,15) = 1.00;
	matA(15,27) = 1.00;
	matA(16,16) = -1.00;
	matA(16,28) = 1.00;
	matA(17,17) = -1.00;
	matA(17,29) = 1.00;
	matA(18,18) = -1.00;
	matA(18,33) = 1.00;
	matA(19,19) = 1.00;
	matA(19,25) = -1.00;
	matA(19,29) = dcX/2.00;
	matA(20,20) = -1.00;
	matA(20,26) = 1.00;
	matA(20,28) = dcX/2.00;
	matA(21,21) = -1.00;
	matA(21,27) = 1.00;
	matA(22,22) = -1.00;
	matA(22,28) = 1.00;
	matA(23,23) = -1.00;
	matA(23,29) = 1.00;
	matA(24,26) = -1.00;
	matA(24,30) = 1.00;
	matA(25,24) = -1.00;
	matA(25,31) = 1.00;
	matA(26,26) = 1.00;
	matA(26,32) = -1.00;
	matA(27,24) = 1.00;
	matA(27,33) = -1.00;
	matA(28,24) = 1.00;
	matA(28,26) = -1.00;
	matA(28,30) = 1.00;
	matA(28,33) = -1.00;
	matA(29,24) = -1.00;
	matA(29,26) = -1.00;
	matA(29,30) = 1.00;
	matA(29,31) = 1.00;
	matA(30,24) = 1.00;
	matA(30,26) = 1.00;
	matA(30,32) = -1.00;
	matA(30,33) = -1.00;
	matA(31,24) = -1.00;
	matA(31,26) = 1.00;
	matA(31,31) = 1.00;
	matA(31,32) = -1.00;
	
	// // // PLANE X-Y
	// matA(0,0) = 1.00;
	// matA(0,24) = -1.00;
	// matA(0,29) = -dcY/2.00;
	// matA(1,1) = -1.00;
	// matA(1,30) = 1.00;
	// matA(2,2) = -1.00;
	// matA(2,26) = 1.00;
	// matA(2,27) = -dcY/2.00;
	// matA(3,3) = 1.00;
	// matA(3,27) = -1.00;
	// matA(4,4) = -1.00;
	// matA(4,28) = 1.00;
	// matA(5,5) = 1.00;
	// matA(5,29) = -1.00;
	// matA(6,6) = 1.00;
	// matA(6,31) = -1.00;
	// matA(7,7) = -1.00;
	// matA(7,25) = 1.00;
	// matA(7,29) = dcX/2.00;
	// matA(8,8) = 1.00;
	// matA(8,26) = -1.00;
	// matA(8,28) = dcX/2.00;
	// matA(9,9) = 1.00;
	// matA(9,27) = -1.00;
	// matA(10,10) = 1.00;
	// matA(10,28) = -1.00;
	// matA(11,11) = 1.00;
	// matA(11,29) = -1.00;
	// matA(12,12) = -1.00;
	// matA(12,24) = 1.00;
	// matA(12,29) = -dcY/2.00;
	// matA(13,13) = 1.00;
	// matA(13,32) = -1.00;
	// matA(14,14) = 1.00;
	// matA(14,26) = -1.00;
	// matA(14,27) = -dcY/2.00;
	// matA(15,15) = -1.00;
	// matA(15,27) = 1.00;
	// matA(16,16) = 1.00;
	// matA(16,28) = -1.00;
	// matA(17,17) = -1.00;
	// matA(17,29) = 1.00;
	// matA(18,18) = -1.00;
	// matA(18,33) = 1.00;
	// matA(19,19) = 1.00;
	// matA(19,25) = -1.00;
	// matA(19,29) = dcX/2.00;
	// matA(20,20) = -1.00;
	// matA(20,26) = 1.00;
	// matA(20,28) = dcX/2.00;
	// matA(21,21) = -1.00;
	// matA(21,27) = 1.00;
	// matA(22,22) = -1.00;
	// matA(22,28) = 1.00;
	// matA(23,23) = -1.00;
	// matA(23,29) = 1.00;
	// matA(24,25) = 1.00;
	// matA(24,30) = -1.00;
	// matA(25,24) = -1.00;
	// matA(25,31) = 1.00;
	// matA(26,25) = -1.00;
	// matA(26,32) = 1.00;
	// matA(27,24) = 1.00;
	// matA(27,33) = -1.00;
	// matA(28,24) = 1.00;
	// matA(28,25) = 1.00;
	// matA(28,30) = -1.00;
	// matA(28,33) = -1.00;
	// matA(29,24) = -1.00;
	// matA(29,25) = 1.00;
	// matA(29,30) = -1.00;
	// matA(29,31) = 1.00;
	// matA(30,24) = 1.00;
	// matA(30,25) = -1.00;
	// matA(30,32) = 1.00;
	// matA(30,33) = -1.00;
	// matA(31,24) = -1.00;
	// matA(31,25) = -1.00;
	// matA(31,31) = 1.00;
	// matA(31,32) = 1.00;

	// opserr << "getmatA: matA" << matA << endln;
	
	// opserr << "getmatA: END" << endln;
}


double Inno3DPnPJoint::getStepSize(double s0,double s1,Vector uExt,Vector duExt,Vector uInt,Vector duInt,double tol)
{
	// opserr << "getStepSize: START" << endln;
	
    Vector u(34);    u.Zero();
    Vector fSpr(32); fSpr.Zero();
    Vector kSpr(32); kSpr.Zero();
    Vector intEq(4); intEq.Zero();

	// tolerance check for line-search
    double r0 = 0.0;
	
	// slack region tolerance set for line-search	
    double tolerance = 0.8;     
    
    if (s0 != 0.0)
	{
        r0 = fabs(s1/s0);
	}

    if (r0 <= tolerance)
	{
		// Linsearch Not required residual decrease less than tolerance
		return 1.0;   
	}

    if (s1 == s0)
	{
		// Bisection would have divide by zero error if continued
        return 1.0;   
	}

    // set some variables
    double etaR;
    double eta = 1.0;
    double s = s1;
    double etaU = 1.0;
    double etaL = 0.0;
    double sU = s1;
    double sL = s0;
    double r = r0;
    double etaJ = 1.0;

    double minEta = 0.1;
    double maxEta = 10.0;
    int maxIter = 1000;

	// opserr << "enterring bracket region " << endln;
    // bracket the region
    int count = 0;
    while ((sU*sL > 0.0) && (etaU < maxEta))
	{
        count = count + 1;
        etaU = etaU * 2.0;
        etaR = etaU - etaJ;

        for (int i = 0; i < 30; i++)
		{
            u(i) = uExt(i) + duExt(i);
        }
        u(30) = uInt(0) - etaR*duInt(0);
		u(31) = uInt(1) - etaR*duInt(1);
		u(32) = uInt(2) - etaR*duInt(2);
		u(33) = uInt(3) - etaR*duInt(3);
		
		// opserr << "getStepSize(1) - u: " << u << endln;

        getMatResponse(u,fSpr,kSpr);
		
		// // internal equilibrium
		// // plane Z-X
		intEq(0) =  fSpr(2)  - fSpr(24) - fSpr(28) - fSpr(29);

		intEq(1) = -fSpr(6)  + fSpr(25) + fSpr(29) + fSpr(31);

		intEq(2) = -fSpr(14) + fSpr(26) + fSpr(30) + fSpr(31);

		intEq(3) =  fSpr(18) - fSpr(27) - fSpr(28) - fSpr(30);

		// // plane Y-X
		// intEq(0) =  fSpr(1)  - fSpr(24) - fSpr(28) - fSpr(29);

		// intEq(1) = -fSpr(6)  + fSpr(25) + fSpr(29) + fSpr(31);

		// intEq(2) = -fSpr(13) + fSpr(26) + fSpr(30) + fSpr(31);

		// intEq(3) =  fSpr(18) - fSpr(27) - fSpr(28) - fSpr(30);

		// opserr << "getStepSize(1) - fSpr: " << fSpr << endln;
		// opserr << "getStepSize(1) - intEq: " << intEq << endln;
		
		// opserr << "getStepSize(1) - duInt: " << duInt << endln;
		// opserr << "getStepSize(1) - intEq: " << intEq << endln;
		
        sU = duInt^intEq;
		
		// opserr << "getStepSize(1) - sU: " << sU << endln;

        // check if the solution is ok
        r = fabs(sU/s0);

        if (r < tolerance)
		{
            return etaU;
		}

        etaJ = etaU;
    }

	// no bracketing could be done
    if (sU*sL > 0.0) 
	{		
        return 1.0;
	}

    count = 0;
    while (r > tolerance && count < maxIter) {
        count = count + 1;
        eta = (etaU + etaL)/2.0;

        if (r > r0)
		{
			eta = 1.0;
		}

        etaR = eta - etaJ;

        for (int i = 0; i < 30; i++) {
            u(i) = uExt(i) + duExt(i);
        }
        u(30) = uInt(0) - etaR*duInt(0);
		u(31) = uInt(1) - etaR*duInt(1);
		u(32) = uInt(2) - etaR*duInt(2);
		u(33) = uInt(3) - etaR*duInt(3);
		
		// opserr << "getStepSize(2) - u: " << u << endln;

        getMatResponse(u,fSpr,kSpr);

		// // internal equilibrium
		// // plane Z-X
		intEq(0) =  fSpr(2)  - fSpr(24) - fSpr(28) - fSpr(29);

		intEq(1) = -fSpr(6)  + fSpr(25) + fSpr(29) + fSpr(31);

		intEq(2) = -fSpr(14) + fSpr(26) + fSpr(30) + fSpr(31);

		intEq(3) =  fSpr(18) - fSpr(27) - fSpr(28) - fSpr(30);

		// // // plane Y-X
		// intEq(0) =  fSpr(1)  - fSpr(24) - fSpr(28) - fSpr(29);

		// intEq(1) = -fSpr(6)  + fSpr(25) + fSpr(29) + fSpr(31);

		// intEq(2) = -fSpr(13) + fSpr(26) + fSpr(30) + fSpr(31);

		// intEq(3) =  fSpr(18) - fSpr(27) - fSpr(28) - fSpr(30);

		// opserr << "getStepSize(2) - fSpr: " << fSpr << endln;
		// opserr << "getStepSize(2) - intEq: " << intEq << endln;
		
        s = duInt^intEq;

        // check if the solution is ok
        r = fabs(s/s0);

        // set variables for next iteration 
        etaJ = eta;
        
        if (s*sU < 0.0)
		{
            etaL = eta;
            sL = s;
		}
		else if (s*sU == 0.0)
		{
			count = maxIter;
		}
		else
		{
            etaU = eta;
            sU = s;
        }

        if (sL == sU)
		{
            count = maxIter;
		}
    }
	
	// opserr << "getStepSize - u: " << u << endln;
	
	// opserr << "getStepSize: END" << endln;
	
    return eta;
}

  
const Matrix& Inno3DPnPJoint::getDamp(void)
{
	// opserr << "getDamp: START" << endln;
    //not applicable (stiffness being returned)
	// opserr << " getDamp!" << endln;	
	K.Zero();
	
	// opserr << "getDamp: END" << endln;
	return K;
}


const Matrix& Inno3DPnPJoint::getMass(void)
{
	// opserr << "getMass: START" << endln;

    //not applicable  (stiffness being returned)
	K.Zero();
	
	// opserr << "getMass: END" << endln;
    return K;
}


void Inno3DPnPJoint::zeroLoad(void)
{
    // not applicable  
	// opserr << "zeroLoad: START/END" << endln;
    return;
}


int Inno3DPnPJoint::addLoad(ElementalLoad *theLoad, double loadFactor)
{
    // not applicable
	// opserr << "addLoad: START/END" << endln;
    return 0;
}


int Inno3DPnPJoint::addInertiaLoadToUnbalance(const Vector &accel)
{
    // not applicable
	// opserr << "addInertiaLoadToUnbalance: START/END" << endln;
	return 0;
}


const Vector& Inno3DPnPJoint::getResistingForceIncInertia()
{
	// not applicable (residual being returned)
    // opserr << "getResistingForceIncInertia: START/END" << endln;
    return R;
}


int Inno3DPnPJoint::sendSelf(int commitTag, Channel &theChannel)
{
	// opserr << "sendSelf: START" << endln;
	int res;
	
	int dataTag = this->getDbTag();
	
	static ID idData(70); // to be confirmed !!!
	
	idData(0) = this->getTag();
	
	idData(1) = ExternalNodes(0);
	idData(2) = ExternalNodes(1);
	idData(3) = ExternalNodes(2);
	idData(4) = ExternalNodes(3);
	idData(5) = ExternalNodes(4);
	
	for (int i=0; i<32; i++)
	{
		// opserr << "sendSelf: FOR i " << i << endln;
		if (MaterialPtr[i] != NULL)
		{
			// opserr << "sendSelf: IF i " << i << endln;
			idData(i+6) = MaterialPtr[i] -> getClassTag();
			// opserr << "sendSelf: idData " << idData(i+6) << endln;
			
			int SpringDbTag = MaterialPtr[i]->getDbTag();
			// opserr << "sendSelf: SpringDbTag " << SpringDbTag << endln;
			
			if (SpringDbTag == 0)
			{
				SpringDbTag = theChannel.getDbTag();
				
				if (SpringDbTag != 0)
				{
					MaterialPtr[i] -> setDbTag(SpringDbTag);
				}
			}
			
			idData(i+38) = SpringDbTag;
			
		}
		else
		{
			// opserr << "sendSelf: ELSE i " << i << endln;
			idData(i+ 5) = 0;
			idData(i+38) = 0;
		}
	}
	
	// opserr << "sendSelf: idData " << idData << endln;
	
	// send the ID vector
	res += theChannel.sendID(dataTag, commitTag, idData);

	if (res < 0)
	{
		opserr << "Inno3DPnPJoint::sendSelf -- failed to send ID data." << endln;
		return -1;
	}
	
	
	for (int i=0 ; i<32 ; i++ )
	{
		if (MaterialPtr[i] != NULL )
		{
			res = MaterialPtr[i]->sendSelf(commitTag, theChannel);
			if (res < 0)
			{
				opserr << "Inno3DPnPJoint::sendSelf() - "<< this->getTag() << " failed to send its Spring " << (i+1) << " material\n";
				return -3;
			}
		}
	}
	
	// opserr << "sendSelf: END" << endln;
	// opserr << "sendSelf: res " << res << endln;
	
    return res;
	
	// return -1;
}


int Inno3DPnPJoint::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
	int res;

	int dataTag = this->getDbTag();

	static ID idData(70); // to be confirmed !!!
	
	res = theChannel.recvID(dataTag, commitTag, idData);
	
	// opserr << "recvID: res " << res << endln;
	
	if (res < 0)
	{
		opserr << "ERROR: Inno3DPnPJoint::recvSelf() - failed to receive ID data." << endln;
		return -1;
	}
	
	this->setTag((int)idData(0));
	
	ExternalNodes(0) = idData(1);
	ExternalNodes(1) = idData(2);
	ExternalNodes(2) = idData(3);
	ExternalNodes(3) = idData(4);
	ExternalNodes(4) = idData(5);

	for (int i=0; i<32; i++)
	{
		int SpringClass = idData(i+ 6);
		int SpringDb    = idData(i+38);
		
		if (SpringClass != 0 && SpringDb != 0)
		{
			// check if we have a material object already & if we do, if of right type
			if ((MaterialPtr[i] == 0) || (MaterialPtr[i]->getClassTag() != SpringClass))
			{
				// if old one .. delete it
				if (MaterialPtr[i] != 0)
				{
					delete MaterialPtr[i];
				}

				// create a new material object

				MaterialPtr[i] = theBroker.getNewUniaxialMaterial(SpringClass);
				
				if (MaterialPtr[i] == 0)
				{
					opserr << "WARNING Inno3DPnPJoint::recvSelf() - " << (i+6) << " failed to get a blank Material of type " << this->getTag() << " for Spring " << SpringClass << endln;
					
					return -3;
				}
			}
			
			MaterialPtr[i]->setDbTag(SpringDb); // note: we set the dbTag before we receive the material
			
			res = MaterialPtr[i]->recvSelf(commitTag, theChannel, theBroker);
			
			if (res < 0)
			{
				opserr << "WARNING Inno3DPnPJoint::recvSelf() - " << this->getTag() << " failed to receive its Material for Spring " << (i+6) << endln;
				
				return -3;
			}
		}		
		else
		{
			MaterialPtr[i] = NULL;
		}
	}
	
	return res;
}

int Inno3DPnPJoint::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode)
{
	// opserr << "displaySelf: START" << endln;
	
	// get node coordinates
    const Vector &node1Crd = nodePtr[0]->getCrds();
    const Vector &node2Crd = nodePtr[1]->getCrds(); 
    const Vector &node3Crd = nodePtr[2]->getCrds();
    const Vector &node4Crd = nodePtr[3]->getCrds();
	const Vector &node5Crd = nodePtr[4]->getCrds();

	// get node disp
    const Vector &node1Disp = nodePtr[0]->getDisp();
    const Vector &node2Disp = nodePtr[1]->getDisp();    
    const Vector &node3Disp = nodePtr[2]->getDisp();
    const Vector &node4Disp = nodePtr[3]->getDisp();
	const Vector &node5Disp = nodePtr[4]->getDisp();
	
	// define vector variables for node coord+disp.
    static Vector v1(3);
    static Vector v2(3);
    static Vector v3(3);
    static Vector v4(3);
	static Vector v5(3);
    
    // calculate the current coordinates of external nodes
    for (int i=0; i<3; i++) 
    {
        v1(i) = node1Crd(i)+node1Disp(i)*fact;
        v2(i) = node2Crd(i)+node2Disp(i)*fact;
        v3(i) = node3Crd(i)+node3Disp(i)*fact;
        v4(i) = node4Crd(i)+node4Disp(i)*fact;
		v5(i) = node5Crd(i)+node5Disp(i)*fact;
    }
	
	// draw the center lines 5-1, 5-2, 5-3, 5-4
	int dummy;
	dummy = theViewer.drawLine(v5, v1, 1.0, 1.0);
	dummy = theViewer.drawLine(v5, v2, 1.0, 1.0);
    dummy = theViewer.drawLine(v5, v3, 1.0, 1.0);
	dummy = theViewer.drawLine(v5, v4, 1.0, 1.0);
	
	
	// // draw diagonals 1-2, 1-4, 3-2, 3-4
	// dummy = theViewer.drawLine(v1, v4, 0.0, 1.0);
	// dummy = theViewer.drawLine(v1, v2, 0.0, 1.0);
	
	// dummy = theViewer.drawLine(v3, v2, 0.0, 1.0);
	// dummy = theViewer.drawLine(v3, v4, 0.0, 1.0);
	
	
	// identify corners of cross-section
	static Vector v12(3);
    static Vector v23(3);
	static Vector v34(3);
	static Vector v41(3);
	
	v12 = v1 + v2;
	v23 = v2 + v3;
	v34 = v3 + v4;
	v41 = v4 + v1;
	
	// draw lines to corners 1-12, 2-12, 2-23, 3-23, 3-34, 4-34, 4-41, 1-41
	dummy = theViewer.drawLine(v1, v12, 1.0, 1.0);
	dummy = theViewer.drawLine(v2, v12, 1.0, 1.0);
	
	dummy = theViewer.drawLine(v2, v23, 1.0, 1.0);
	dummy = theViewer.drawLine(v3, v23, 1.0, 1.0);
	
	dummy = theViewer.drawLine(v3, v34, 1.0, 1.0);
	dummy = theViewer.drawLine(v4, v34, 1.0, 1.0);
	
	dummy = theViewer.drawLine(v4, v41, 1.0, 1.0);
	dummy = theViewer.drawLine(v1, v41, 1.0, 1.0);
	
	// opserr << "displaySelf: END" << endln;
	return 0;
}


void Inno3DPnPJoint::Print(OPS_Stream &s, int flag)
{
	// opserr << "Print: START" << endln;
	
	opserr << "Inno3DJointND element - Written by Cristian V. Miculas @ISISE - University of Coimbra, Portugal." << endln;
	
	opserr << "Element Tag: " << this->getTag() << "; Type: Beam Column Inno3DPnPJoint." << endln;
	
	for (int i = 0; i<5; i++)
	{
		opserr << "Node: " << ExternalNodes(i)  << "; ";
		opserr << "DOFs: " << nodePtr[i]->getNumberDOF() << "; ";
		opserr << "Coordinates: " << nodePtr[i]->getCrds();
	}
	
	opserr << "Cross-section dcX x dcZ (width X height Local X-Z): " << dcX << " x "<< dcZ << endln;
	
	// opserr << "Print: END" << endln;	
	return;
}


Response* Inno3DPnPJoint::setResponse(const char **argv, int argc, OPS_Stream &output)
{
	// opserr << "setResponse: START" << endln;
	
    // // // OUTPUT AT COMPONENT LEVEL
    /////////////////////////////////////////////////////////////////////
	// material response for the 32 springs/components
	if ( (strcmp(argv[0],"spring")==0) || (strcmp(argv[0],"-spring") == 0) ||
	    (strcmp(argv[0],"material")==0) || (strcmp(argv[0],"-material") == 0) )
		{	
			// // prints number of args. should be 3
            // opserr << "setResponse: argc1: " << argc << endln;
			
            // // prints each argument (the 3 from above)
            // for (int ii = 0; ii < argc; ii ++)
				// {
					// opserr << "setResponse: argv[" << ii << "]: " << argv[ii] << endln;
				// }
			
            // parses the "str", interprets its content "int"
            int springNo = atoi(argv[1]);
			// opserr << "setResponse: springNo: " << springNo << endln;
				
			if (springNo > 0 && springNo < 33)
            {
				if (MaterialPtr[springNo-1] != 0)
                {
					// opserr << "setResponse: output for spring " << springNo << " done!" << endln;
                    
					return MaterialPtr[springNo-1]->setResponse(&argv[2], argc-1, output);
                }
            }
            else
            {
                opserr << "ERROR: Inno3DPnPJoint::setResponse number of springs out of range: "<< springNo << endln;
                opserr << "Spring numbers go from 1 to 32." << endln;
                return 0;
            }
		}


    // // // OUTPUT AT ELEMENT LEVEL
    /////////////////////////////////////////////////////////////////////
    
    // compare argv[0] for known response types for the Truss
	// this goes with the cases from getResponse
	
    // displacement for external DOFs (30)
	else if (strcmp(argv[0],"extDisp") == 0 || strcmp(argv[0],"extdisp") == 0)
    {
        return new ElementResponse(this,1,Vector(30));
	}
	
    // displacement for internal DOFs (4)
	else if (strcmp(argv[0],"intDisp") == 0 || strcmp(argv[0],"intdisp") == 0)
    {
        return new ElementResponse(this,2,Vector(4));
	}
	
    // displacement for external and internal DOFs (30 + 4)
	else if (strcmp(argv[0],"disp") == 0 || strcmp(argv[0],"Disp") == 0)
    {
		return new ElementResponse(this,3,Vector(34));
	}
	
    // residual force for external and internal DOFs (30 + 4)
	else if (strcmp(argv[0],"reaction") == 0 || strcmp(argv[0],"Reaction") == 0)
    {
        return new ElementResponse(this,4,Vector(34));
	}
	
	// material output for springs/components at element level (prints out all)
    
	// force/bmom || stress
	else if (strcmp(argv[0],"matStress") == 0 || strcmp(argv[0],"matstress") == 0 || strcmp(argv[0],"stress") == 0 || strcmp(argv[0],"Stress") == 0)
    {
		return new ElementResponse(this,5,Vector(32));
	}
	
	// disp/rot || strain
	else if (strcmp(argv[0],"matStrain") == 0 || strcmp(argv[0],"matstrain") == 0 || strcmp(argv[0],"strain") == 0 || strcmp(argv[0],"Strain") == 0)
    {
		return new ElementResponse(this,6,Vector(32));
	}
	
	// force/bmom + disp/rot || stress + strain
	else if (strcmp(argv[0],"matStressStrain") == 0 || strcmp(argv[0],"matstressstrain") == 0 || strcmp(argv[0],"stressStrain") == 0 || strcmp(argv[0],"StressStrain") == 0 )
    {
		return new ElementResponse(this,7,Vector(64));
	}
	
	else
    {
		return 0;
	}
	
	// opserr << "setResponse: END" << endln;
}


int Inno3DPnPJoint::getResponse(int responseID, Information &eleInfo)
{
	// opserr << "getResponse: START" << endln;
	
	switch (responseID) {
	case -1:
		return -1;

	case 1:       
		if(eleInfo.theVector!=0)
		{
			// transformation from LOCAL to GLOBAL (EXT disp)
			Vector UeprCommit_G(30);
			UeprCommit_G.Zero();
			
			UeprCommit_G.addMatrixTransposeVector(0.0, Transf, UeprCommit, 1.0);
			
			// opserr << "getResponse -- case 1: " << UeprCommit_G << endln;
			// opserr << "getResponse -- case 1: " << UeprCommit << endln;
			
			for (int j=0; j<30; j++)
				{
					(*(eleInfo.theVector))(j) =  UeprCommit_G(j);
				}
		}
		// opserr << "getResponse: END " << "case 1" <<endln;
		return 0;
		
	case 2:       
		if(eleInfo.theVector!=0)
		{
			for (int j=0; j<4; j++)
				{
					(*(eleInfo.theVector))(j) =  UeprIntCommit(j);

				}
		}
		// opserr << "getResponse: END " << "case 2" <<endln;
		return 0;
		
	case 3:       
		if(eleInfo.theVector!=0)
		{
			// transformation from LOCAL to GLOBAL (EXT disp)
			Vector UeprCommit_G(30);
			UeprCommit_G.Zero();
			
			UeprCommit_G.addMatrixTransposeVector(0.0, Transf, UeprCommit, 1.0);
			
			// opserr << "getResponse -- case 1: " << UeprCommit_G <<endln;
			// opserr << "getResponse -- case 1: " << UeprCommit <<endln;
			
			for (int j=0; j<30; j++)
				{
					(*(eleInfo.theVector))(j) =  UeprCommit_G(j);
				}

			// no transformation is needed (INT disp)
			for (int k=0; k<4; k++)
				{
					(*(eleInfo.theVector))(30+k) =  UeprIntCommit(k);

				}
				
		}
		// opserr << "getResponse: END " << "case 3" <<endln;
		return 0;

	case 4:       
		if(eleInfo.theVector!=0)
		{
			for (int j=0; j<34; j++)
				{
					(*(eleInfo.theVector))(j) =  R(j);
				}
		}
		// opserr << "getResponse: END " << "case 4" <<endln;
		return 0;

	case 5:
		if(eleInfo.theVector!=0)
		{
			for ( int i =0 ; i<32 ; i++ )
			{
				(*(eleInfo.theVector))(i) = 0.0;
				if ( MaterialPtr[i] != NULL ) 
					(*(eleInfo.theVector))(i) = MaterialPtr[i]->getStress();
			}
		}
		// opserr << "getResponse: END " << "case 5" <<endln;
		return 0;
		
	case 6:
		if(eleInfo.theVector!=0)
		{
			for ( int i =0 ; i<32 ; i++ )
			{
				(*(eleInfo.theVector))(i) = 0.0;
				if ( MaterialPtr[i] != NULL ) 
					(*(eleInfo.theVector))(i) = MaterialPtr[i]->getStrain();
			}
		}
		// opserr << "getResponse: END " << "case 6" <<endln;
		return 0;

	case 7:
		if(eleInfo.theVector!=0)
		{
			for ( int i =0 ; i<32 ; i++ )
			{
				(*(eleInfo.theVector))(i) = 0.0;
				(*(eleInfo.theVector))(i+32) = 0.0;
				if ( MaterialPtr[i] != NULL )
				{
					(*(eleInfo.theVector))(i) 	 = MaterialPtr[i]->getStrain();
					(*(eleInfo.theVector))(i+32) = MaterialPtr[i]->getStress();
				}
			}
		}
		// opserr << "getResponse: END " << "case 7" <<endln;
		return 0;

	// case 6:
		// return eleInfo.setMatrix(this->getTangentStiff());
	
	// case 7:
		// if(eleInfo.theVector!=0)
		// {
			// for ( int i=0 ; i<38 ; i++ )
			// {
				// (*(eleInfo.theVector))(i) = 0.0;
				// if ( MaterialPtr[i] != NULL && MaterialPtr[i]->getInitialTangent() != 0.0 )
				// {
					// (*(eleInfo.theVector))(i) = 
						// MaterialPtr[i]->getStrain() - MaterialPtr[i]->getStress()/MaterialPtr[i]->getInitialTangent();
				// }
				
			// }			
		// }
		// return 0;

	default:
		// opserr << "getResponse: END " << "default return -1" <<endln;
		return -1;
	}
	// opserr << "getResponse: END " << "return -1" <<endln;
	return -1;
	
	// opserr << "getResponse: END" << endln;
}


int Inno3DPnPJoint::setParameter (char **argv, int argc, Information &info)
{
	// opserr << "setParameter: START/END" << endln;
	return -1;
}
 
 
int Inno3DPnPJoint::updateParameter (int parameterID, Information &info)
{
	// opserr << "updateParameter: START/END" << endln;
	return -1;
}