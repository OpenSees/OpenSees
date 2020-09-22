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

// Written: csasj
// $Revision: 1.11 $
// $Date: 21/09/2020 $

// Description: This file contains the implementation of the
// TclModelBuilder_addPISAMaterial() function. 

#include <TclModelBuilder.h>
#include <PyClayPISA.h> // csasj
#include <MtClayPISA.h> // csasj
#include <HbClayPISA.h> // csasj
#include <MbClayPISA.h> // csasj

#include <PySandPISA.h> // csasj
#include <MtSandPISA.h> // csasj
#include <HbSandPISA.h> // csasj
#include <MbSandPISA.h> // csasj

#include <tcl.h>

#include <Vector.h>
#include <string.h>
#include <fstream>
using std::ifstream;

#include <iomanip>
using std::ios;



static void printCommand(int argc, TCL_Char** argv)
{
	opserr << "Input command: ";
	for (int i = 0; i < argc; i++)
		opserr << argv[i] << " ";
	opserr << endln;
}

UniaxialMaterial*
TclModelBuilder_addPISAMaterial(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char** argv)
{
	// Make sure there is a minimum number of arguments
	if (argc < 3) {
		opserr << "WARNING insufficient number of arguments\n";
		printCommand(argc, argv);
		return 0;
	}

	int tag;
	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
		opserr << "WARNING invalid uniaxialMaterial tag\n";
		printCommand(argc, argv);
		return 0;
	}

	// Pointer to a uniaxial material that will be added to the model builder
	UniaxialMaterial* theMaterial = 0;

	// Check argv[1] for soil curve component
	if (strcmp(argv[1], "PyClayPISA") == 0) {

		if (argc < 8) {
			opserr << "WARNING insufficient arguments\n";
			printCommand(argc, argv);
			opserr << "Want: uniaxialMaterial PyClayPISA tag? depth? diameter? undrained shear strength? G0? grid spacing? A0? Au? " << endln;
			return 0;
		}

		int tag;
		double depth, diameter, undrained_shear_strength, gmax, grid, A0, Au;

		if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
			opserr << "WARNING invalid uniaxialMaterial PyClayPISA tag" << endln;
			return 0;
		}
		if (Tcl_GetDouble(interp, argv[3], &depth) != TCL_OK) {
			opserr << "WARNING invalid depth (m)\n";
			opserr << "uniaxialMaterial PyClayPISA: " << tag << endln;
			return 0;
		}
		if (Tcl_GetDouble(interp, argv[4], &diameter) != TCL_OK) {
			opserr << "WARNING invalid diameter\n";
			opserr << "uniaxialMaterial PyClayPISA: " << tag << endln;
			return 0;
		}
		if (Tcl_GetDouble(interp, argv[5], &undrained_shear_strength) != TCL_OK) {
			opserr << "WARNING invalid undrained shear strength (kPa)\n";
			opserr << "uniaxialMaterial PyClayPISA: " << tag << endln;
			return 0;
		}
		if (Tcl_GetDouble(interp, argv[6], &gmax) != TCL_OK) {
			opserr << "WARNING invalid small-strain shear modulus (kPa)\n";
			opserr << "uniaxialMaterial PyClayPISA: " << tag << endln;
			return 0;
		}
		if (Tcl_GetDouble(interp, argv[7], &grid) != TCL_OK) {
			opserr << "WARNING invalid grid (m)\n";
			opserr << "uniaxialMaterial PyClayPISA: " << tag << endln;
			return 0;
		}

		if (argc == 8)
		{
			A0 = 1.0;
			Au = 1.0;
		}

		if (argc > 8)
		{ 
			if (Tcl_GetDouble(interp, argv[8], &A0) != TCL_OK) {
				opserr << "WARNING invalid scaling factor A0\n";
				opserr << "uniaxialMaterial PyClayPISA: " << tag << endln;
				return 0;
			}
			else if  (Tcl_GetDouble(interp, argv[9], &Au) != TCL_OK) {
				opserr << "WARNING invalid scaling factor Au\n";
				opserr << "uniaxialMaterial PyClayPISA: " << tag << endln;
				return 0;
			}
		}
	
		theMaterial = new PyClayPISA (tag, depth, diameter, undrained_shear_strength, gmax, grid, A0, Au);
	}
	
	else if (strcmp(argv[1], "MtClayPISA") == 0) {

		if (argc < 8) {
			opserr << "WARNING insufficient arguments\n";
			printCommand(argc, argv);
			opserr << "Want: uniaxialMaterial MtClayPISA tag? depth? diameter? undrained shear strength? G0? grid spacing? A0? Au? " << endln;
			return 0;
		}

		int tag;
		double depth, diameter, undrained_shear_strength, gmax, grid, A0, Au;

		if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
			opserr << "WARNING invalid uniaxialMaterial PySimple1 tag" << endln;
			return 0;
		}
		if (Tcl_GetDouble(interp, argv[3], &depth) != TCL_OK) {
			opserr << "WARNING invalid depth (m)\n";
			opserr << "uniaxialMaterial MtClayPISA: " << tag << endln;
			return 0;
		}
		if (Tcl_GetDouble(interp, argv[4], &diameter) != TCL_OK) {
			opserr << "WARNING invalid diameter\n";
			opserr << "uniaxialMaterial MtClayPISA: " << tag << endln;
			return 0;
		}
		if (Tcl_GetDouble(interp, argv[5], &undrained_shear_strength) != TCL_OK) {
			opserr << "WARNING invalid undrained shear strength (kPa)\n";
			opserr << "uniaxialMaterial MtClayPISA: " << tag << endln;
			return 0;
		}
		if (Tcl_GetDouble(interp, argv[6], &gmax) != TCL_OK) {
			opserr << "WARNING invalid small-strain shear modulus (kPa)\n";
			opserr << "uniaxialMaterial MtClayPISA: " << tag << endln;
			return 0;
		}
		if (Tcl_GetDouble(interp, argv[7], &grid) != TCL_OK) {
			opserr << "WARNING invalid grid (m)\n";
			opserr << "uniaxialMaterial MtClayPISA: " << tag << endln;
			return 0;
		}

		if (argc == 8)
		{
			A0 = 1.0;
			Au = 1.0;
		}

		if (argc > 8)
		{
			if (Tcl_GetDouble(interp, argv[8], &A0) != TCL_OK) {
				opserr << "WARNING invalid scaling factor A0\n";
				opserr << "uniaxialMaterial MtClayPISA: " << tag << endln;
				return 0;
			}
			else if (Tcl_GetDouble(interp, argv[9], &Au) != TCL_OK) {
				opserr << "WARNING invalid scaling factor Au\n";
				opserr << "uniaxialMaterial MtClayPISA: " << tag << endln;
				return 0;
			}
		}

		theMaterial = new MtClayPISA(tag, depth, diameter, undrained_shear_strength, gmax, grid, A0, Au);
	}

	else if (strcmp(argv[1], "HbClayPISA") == 0)  {

	if (argc < 7) {
		opserr << "WARNING insufficient arguments\n";
		printCommand(argc, argv);
		opserr << "Want: uniaxialMaterial HbClayPISA tag? diameter? embedded length? undrained shear strength? G0? A0? Au? " << endln;
		return 0;
	}

	int tag;
	double  diameter, embedded_length, undrained_shear_strength, gmax, A0, Au;

	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
		opserr << "WARNING invalid uniaxialMaterial HbClayPISA tag" << endln;
		return 0;
	}
	if (Tcl_GetDouble(interp, argv[3], &diameter) != TCL_OK) {
		opserr << "WARNING invalid diameter (m)\n";
		opserr << "uniaxialMaterial HbClayPISA: " << tag << endln;
		return 0;
	}
	if (Tcl_GetDouble(interp, argv[4], &embedded_length) != TCL_OK) {
		opserr << "WARNING invalid embedded_length (m) \n";
		opserr << "uniaxialMaterial HbClayPISA: " << tag << endln;
		return 0;
	}
	if (Tcl_GetDouble(interp, argv[5], &undrained_shear_strength) != TCL_OK) {
		opserr << "WARNING invalid undrained shear strength (kPa)\n";
		opserr << "uniaxialMaterial HbClayPISA: " << tag << endln;
		return 0;
	}
	if (Tcl_GetDouble(interp, argv[6], &gmax) != TCL_OK) {
		opserr << "WARNING invalid small-strain shear modulus (kPa)\n";
		opserr << "uniaxialMaterial HbClayPISA: " << tag << endln;
		return 0;
	}

	if (argc == 7)
	{
		A0 = 1.0;
		Au = 1.0;
	}

	if (argc > 7)
	{
		if (Tcl_GetDouble(interp, argv[7], &A0) != TCL_OK) {
			opserr << "WARNING invalid scaling factor A0\n";
			opserr << "uniaxialMaterial HbClayPISA: " << tag << endln;
			return 0;
		}
		else if (Tcl_GetDouble(interp, argv[8], &Au) != TCL_OK) {
			opserr << "WARNING invalid scaling factor Au\n";
			opserr << "uniaxialMaterial HbClayPISA: " << tag << endln;
			return 0;
		}
	}

	theMaterial = new HbClayPISA (tag, diameter, embedded_length, undrained_shear_strength, gmax, A0, Au);

	}
	
	else if (strcmp(argv[1], "MbClayPISA") == 0)  {

		if (argc < 7) {
			opserr << "WARNING insufficient arguments\n";
			printCommand(argc, argv);
			opserr << "Want: uniaxialMaterial MbClayPISA tag? diameter? embedded length? undrained shear strength? G0? A0? Au? " << endln;
			return 0;
		}

		int tag;
		double  diameter, embedded_length, undrained_shear_strength, gmax, A0, Au;

		if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
			opserr << "WARNING invalid uniaxialMaterial MbClayPISA tag" << endln;
			return 0;
		}
		if (Tcl_GetDouble(interp, argv[3], &diameter) != TCL_OK) {
			opserr << "WARNING invalid diameter (m)\n";
			opserr << "uniaxialMaterial MbClayPISA: " << tag << endln;
			return 0;
		}
		if (Tcl_GetDouble(interp, argv[4], &embedded_length) != TCL_OK) {
			opserr << "WARNING invalid embedded_length (m) \n";
			opserr << "uniaxialMaterial MbClayPISA: " << tag << endln;
			return 0;
		}
		if (Tcl_GetDouble(interp, argv[5], &undrained_shear_strength) != TCL_OK) {
			opserr << "WARNING invalid undrained shear strength (kPa)\n";
			opserr << "uniaxialMaterial MbClayPISA: " << tag << endln;
			return 0;
		}
		if (Tcl_GetDouble(interp, argv[6], &gmax) != TCL_OK) {
			opserr << "WARNING invalid small-strain shear modulus (kPa)\n";
			opserr << "uniaxialMaterial MbClayPISA: " << tag << endln;
			return 0;
		}

		if (argc == 7)
		{
			A0 = 1.0;
			Au = 1.0;
		}

		if (argc > 7)
		{
			if (Tcl_GetDouble(interp, argv[7], &A0) != TCL_OK) {
				opserr << "WARNING invalid scaling factor A0\n";
				opserr << "uniaxialMaterial MbClayPISA: " << tag << endln;
				return 0;
			}
			else if (Tcl_GetDouble(interp, argv[8], &Au) != TCL_OK) {
				opserr << "WARNING invalid scaling factor Au\n";
				opserr << "uniaxialMaterial MbClayPISA: " << tag << endln;
				return 0;
			}
		}

		theMaterial = new MbClayPISA(tag, diameter, embedded_length, undrained_shear_strength, gmax, A0, Au);
	}

	else if (strcmp(argv[1], "PySandPISA") == 0) {

		if (argc < 9) {
			opserr << "WARNING insufficient arguments\n";
			printCommand(argc, argv);
			opserr << " Want: uniaxialMaterial PySandPISA tag? depth? diameter? pile embedded lenght? sigma_vo_eff? G0?  grid spacing A0? Au? ? " << endln;
			return 0;
		}

		int tag;
		double depth, diameter, embedded_length, sigma_vo_eff, gmax, grid, A0, Au;

		if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
			opserr << "WARNING invalid uniaxialMaterial PySandPISA tag" << endln;
			return 0;
		}
		if (Tcl_GetDouble(interp, argv[3], &depth) != TCL_OK) {
			opserr << "WARNING invalid depth (m)\n";
			opserr << "uniaxialMaterial PySandPISA: " << tag << endln;
			return 0;
		}
		if (Tcl_GetDouble(interp, argv[4], &diameter) != TCL_OK) {
			opserr << "WARNING invalid diameter\n";
			opserr << "uniaxialMaterial PySandPISA: " << tag << endln;
			return 0;
		}
		if (Tcl_GetDouble(interp, argv[5], &embedded_length) != TCL_OK) {
			opserr << "WARNING invalid pile embedded length (m)\n";
			opserr << "uniaxialMaterial PySandPISA: " << tag << endln;
			return 0;
		}
		if (Tcl_GetDouble(interp, argv[6], &sigma_vo_eff) != TCL_OK) {
			opserr << "WARNING invalid vertical effective soil stress (kPa)\n";
			opserr << "uniaxialMaterial PySandPISA: " << tag << endln;
			return 0;
		}
		if (Tcl_GetDouble(interp, argv[7], &gmax) != TCL_OK) {
			opserr << "WARNING invalid small-strain shear modulus (kPa)\n";
			opserr << "uniaxialMaterial PySandPISA: " << tag << endln;
			return 0;
		}
		if (Tcl_GetDouble(interp, argv[8], &grid) != TCL_OK) {
			opserr << "WARNING invalid grid (m)\n";
			opserr << "uniaxialMaterial PySandPISA: " << tag << endln;
			return 0;
		}

		if (argc == 9)
		{
			A0 = 1.0;
			Au = 1.0;
		}

		if (argc > 9)
		{
			if (Tcl_GetDouble(interp, argv[9], &A0) != TCL_OK) {
				opserr << "WARNING invalid scaling factor A0\n";
				opserr << "uniaxialMaterial PySandPISA: " << tag << endln;
				return 0;
			}
			else if (Tcl_GetDouble(interp, argv[10], &Au) != TCL_OK) {
				opserr << "WARNING invalid scaling factor Au\n";
				opserr << "uniaxialMaterial PySandPISA: " << tag << endln;
				return 0;
			}
		}

		theMaterial = new PySandPISA(tag, depth, diameter, embedded_length, sigma_vo_eff, gmax, grid, A0, Au);
	}
	
	else if (strcmp(argv[1], "MtSandPISA") == 0) {
		
	if (argc < 10) {
		opserr << "WARNING insufficient arguments\n";
		printCommand(argc, argv);
		opserr << " Want: uniaxialMaterial MtSandPISA tag? depth? diameter? pile embedded lenght? sigma_vo_eff? G0? Local distributed lateral load? grid spacing? A0? Au? ? " << endln;
		return 0;
	}

	int tag;
	double depth, diameter, embedded_length, sigma_vo_eff, gmax, Pvalue, grid, A0, Au;

	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
		opserr << "WARNING invalid uniaxialMaterial MtSandPISA tag" << endln;
		return 0;
	}
	if (Tcl_GetDouble(interp, argv[3], &depth) != TCL_OK) {
		opserr << "WARNING invalid depth (m)\n";
		opserr << "uniaxialMaterial MtSandPISA: " << tag << endln;
		return 0;
	}
	if (Tcl_GetDouble(interp, argv[4], &diameter) != TCL_OK) {
		opserr << "WARNING invalid diameter\n";
		opserr << "uniaxialMaterial MtSandPISA: " << tag << endln;
		return 0;
	}
	if (Tcl_GetDouble(interp, argv[5], &embedded_length) != TCL_OK) {
		opserr << "WARNING invalid pile embedded length (m)\n";
		opserr << "uniaxialMaterial MtSandPISA: " << tag << endln;
		return 0;
	}
	if (Tcl_GetDouble(interp, argv[6], &sigma_vo_eff) != TCL_OK) {
		opserr << "WARNING invalid vertical effective soil stress (kPa)\n";
		opserr << "uniaxialMaterial MtSandPISA: " << tag << endln;
		return 0;
	}
	if (Tcl_GetDouble(interp, argv[7], &gmax) != TCL_OK) {
		opserr << "WARNING invalid small-strain shear modulus (kPa)\n";
		opserr << "uniaxialMaterial MtSandPISA: " << tag << endln;
		return 0;
	}
	if (Tcl_GetDouble(interp, argv[8], &Pvalue) != TCL_OK) {
		opserr << "WARNING invalid current value of distributed lateral load (kN/m)\n";
		opserr << "uniaxialMaterial MtSandPISA: " << tag << endln;
		return 0;
	}
	if (Tcl_GetDouble(interp, argv[9], &grid) != TCL_OK) {
		opserr << "WARNING invalid grid (m)\n";
		opserr << "uniaxialMaterial MtSandPISA: " << tag << endln;
		return 0;
	}

	if (argc == 10)
	{
		A0 = 1.0;
		Au = 1.0;
	}

	if (argc > 10)
	{
		if (Tcl_GetDouble(interp, argv[10], &A0) != TCL_OK) {
			opserr << "WARNING invalid scaling factor A0\n";
			opserr << "uniaxialMaterial MtSandPISA: " << tag << endln;
			return 0;
		}
		else if (Tcl_GetDouble(interp, argv[11], &Au) != TCL_OK) {
			opserr << "WARNING invalid scaling factor Au\n";
			opserr << "uniaxialMaterial MtSandPISA: " << tag << endln;
			return 0;
		}
	}

	theMaterial = new MtSandPISA(tag, depth, diameter, embedded_length, sigma_vo_eff, gmax, Pvalue, grid, A0, Au);

	}

	else if (strcmp(argv[1], "HbSandPISA") == 0) {

	if (argc < 7) {
		opserr << "WARNING insufficient arguments\n";
		printCommand(argc, argv);
		opserr << "Want: uniaxialMaterial HbSandPISA tag? diameter? pile embedded lenght? sigma_vo_eff? G0? A0? Au? " << endln;
		return 0;
	}

	int tag;
	double  diameter, embedded_length, sigma_vo_eff, gmax, A0, Au;

	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
		opserr << "WARNING invalid uniaxialMaterial HbSandPISA tag" << endln;
		return 0;
	}
	if (Tcl_GetDouble(interp, argv[3], &diameter) != TCL_OK) {
		opserr << "WARNING invalid diameter (m)\n";
		opserr << "uniaxialMaterial HbSandPISA: " << tag << endln;
		return 0;
	}
	if (Tcl_GetDouble(interp, argv[4], &embedded_length) != TCL_OK) {
		opserr << "WARNING invalid embedded_length (m) \n";
		opserr << "uniaxialMaterial HbSandPISA: " << tag << endln;
		return 0;
	}
	if (Tcl_GetDouble(interp, argv[5], &sigma_vo_eff) != TCL_OK) {
		opserr << "WARNING invalid vertical effective soil stress (kPa)\n";
		opserr << "uniaxialMaterial HbSandPISA: " << tag << endln;
		return 0;
	}
	if (Tcl_GetDouble(interp, argv[6], &gmax) != TCL_OK) {
		opserr << "WARNING invalid small-strain shear modulus (kPa)\n";
		opserr << "uniaxialMaterial HbSandPISA: " << tag << endln;
		return 0;
	}

	if (argc == 7)
	{
		A0 = 1.0;
		Au = 1.0;
	}

	if (argc > 7)
	{
		if (Tcl_GetDouble(interp, argv[7], &A0) != TCL_OK) {
			opserr << "WARNING invalid scaling factor A0\n";
			opserr << "uniaxialMaterial HbSandPISA: " << tag << endln;
			return 0;
		}
		else if (Tcl_GetDouble(interp, argv[8], &Au) != TCL_OK) {
			opserr << "WARNING invalid scaling factor Au\n";
			opserr << "uniaxialMaterial HbSandPISA: " << tag << endln;
			return 0;
		}
	}

	theMaterial = new HbSandPISA(tag, diameter, embedded_length, sigma_vo_eff, gmax, A0, Au);

	}

	else if (strcmp(argv[1], "MbSandPISA") == 0) {

	if (argc < 7) {
		opserr << "WARNING insufficient arguments\n";
		printCommand(argc, argv);
		opserr << "Want: uniaxialMaterial MbSandPISA tag? diameter? pile embedded lenght? sigma_vo_eff? G0? A0? Au? " << endln;
		return 0;
	}

	int tag;
	double  diameter, embedded_length, sigma_vo_eff, gmax, A0, Au;

	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
		opserr << "WARNING invalid uniaxialMaterial MbSandPISA tag" << endln;
		return 0;
	}
	if (Tcl_GetDouble(interp, argv[3], &diameter) != TCL_OK) {
		opserr << "WARNING invalid diameter (m)\n";
		opserr << "uniaxialMaterial MbSandPISA: " << tag << endln;
		return 0;
	}
	if (Tcl_GetDouble(interp, argv[4], &embedded_length) != TCL_OK) {
		opserr << "WARNING invalid embedded_length (m) \n";
		opserr << "uniaxialMaterial MbSandPISA: " << tag << endln;
		return 0;
	}
	if (Tcl_GetDouble(interp, argv[5], &sigma_vo_eff) != TCL_OK) {
		opserr << "WARNING invalid vertical effective soil stress (kPa)\n";
		opserr << "uniaxialMaterial MbSandPISA: " << tag << endln;
		return 0;
	}
	if (Tcl_GetDouble(interp, argv[6], &gmax) != TCL_OK) {
		opserr << "WARNING invalid small-strain shear modulus (kPa)\n";
		opserr << "uniaxialMaterial MbSandPISA: " << tag << endln;
		return 0;
	}

	if (argc == 7)
	{
		A0 = 1.0;
		Au = 1.0;
	}

	if (argc > 7){
		if (Tcl_GetDouble(interp, argv[7], &A0) != TCL_OK) {
			opserr << "WARNING invalid scaling factor A0\n";
			opserr << "uniaxialMaterial MbSandPISA: " << tag << endln;
			return 0;
		}
		else if (Tcl_GetDouble(interp, argv[8], &Au) != TCL_OK) {
			opserr << "WARNING invalid scaling factor Au\n";
			opserr << "uniaxialMaterial MbSandPISA: " << tag << endln;
			return 0;
		}
	}

	theMaterial = new MbSandPISA(tag, diameter, embedded_length, sigma_vo_eff, gmax, A0, Au);

	}

	return theMaterial;
}