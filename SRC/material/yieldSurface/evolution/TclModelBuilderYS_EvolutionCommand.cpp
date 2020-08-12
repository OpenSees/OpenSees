#include "YieldSurface_BC.h"
#include <TclModelBuilder.h>
#include <string.h>
#include <Vector.h>

#include "NullEvolution.h"
#include "Kinematic2D01.h"
#include "PeakOriented2D01.h"
#include "Isotropic2D01.h"
#include "CombinedIsoKin2D01.h"

#include "Kinematic2D02.h"
#include "PeakOriented2D02.h"
#include "CombinedIsoKin2D02.h"

#include "PlasticHardeningMaterial.h"
#include "YieldSurface_BC.h"



int addTclYS_Evolution(TclModelBuilder *theBuilder, YS_Evolution *theModel)
{
	if(theModel ==0)
		return TCL_ERROR;

	if(!theModel)
	{
	 	opserr << "Model not created\n";
	    return TCL_ERROR;
	}

	if (theBuilder->addYS_EvolutionModel(*theModel) < 0)
	{
		opserr << "WARNING could not add hardening model to the domain\n";
		opserr << *theModel << endln;
		delete theModel; // invoke the material objects destructor, otherwise mem leak
		return TCL_ERROR;
	}

	return TCL_OK;

}

PlasticHardeningMaterial * getTclPlasticMaterial(Tcl_Interp *interp, TCL_Char *arg, TclModelBuilder *theBuilder)
{
int id;
	if (Tcl_GetInt(interp, arg, &id) != TCL_OK)
	{
		opserr << "WARNING: TclModelYS_EvolutionCommand - Invalid plastic material tag \n";
		return 0;
	}

	PlasticHardeningMaterial *theMat = theBuilder->getPlasticMaterial(id);
	if(theMat == 0)
	{
		opserr << "WARNING: TclModelYS_EvolutionCommand - no PlasticHardeningMaterial with id = "
			 << id << " exists\n";
		return 0;
	}
	else
		return theMat;
}

YieldSurface_BC * getTclYieldSurface_BC(Tcl_Interp *interp, TCL_Char *arg, TclModelBuilder *theBuilder)
{
int id;
	if (Tcl_GetInt(interp, arg, &id) != TCL_OK)
	{
		opserr << "WARNING: TclModelYS_EvolutionCommand - Invalid YieldSurface_BC tag \n";
		return 0;
	}

	YieldSurface_BC *theYS = theBuilder->getYieldSurface_BC(id);
	if(theYS == 0)
	{
		opserr << "WARNING: TclModelYS_EvolutionCommand - no YieldSurface_BC with id = "
			 << id << " exists\n";
		return 0;
	}
	else
		return theYS;
}

//////////////////////////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////////////////////////

int TclNullEvolutionCommand(ClientData clienData, Tcl_Interp *interp, int argc,
				 TCL_Char **argv, TclModelBuilder *theBuilder)
{
YS_Evolution *theModel = 0;
		int tag;
		double isox;
		double isoy;
		double isoz;
		int dim=0;

		if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK)
			return TCL_ERROR;

		if(argc > 3)
		{
			if (Tcl_GetDouble(interp, argv[3], &isox) != TCL_OK)
				return TCL_ERROR;
			dim++;
		}
		if(argc > 4)
		{
			if (Tcl_GetDouble(interp, argv[4], &isoy) != TCL_OK)
				return TCL_ERROR;
			dim++;
		}
		if(argc > 5)
		{
			if (Tcl_GetDouble(interp, argv[5], &isoz) != TCL_OK)
				return TCL_ERROR;
			dim++;
		}

//		opserr << "Dim = " << dim << endln;
//		opserr << "\a";
		
		// Parsing was successful, allocate the material
		if(dim==1)
			theModel = new NullEvolution(tag, isox);
		else if(dim == 2)
			theModel = new NullEvolution(tag, isox, isoy);
		else if(dim == 3)
			theModel = new NullEvolution(tag, isox, isoy, isoz);
		else
			theModel = 0;

return addTclYS_Evolution(theBuilder, theModel);
}

int TclKinematic2D01Command(ClientData clienData, Tcl_Interp *interp, int argc,
				 TCL_Char **argv, TclModelBuilder *theBuilder)
{
YS_Evolution *theModel = 0;
		int tag;
		double minIsoFactor;

		if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK)
			return TCL_ERROR;

		if (Tcl_GetDouble(interp, argv[3], &minIsoFactor) != TCL_OK)
			return TCL_ERROR;

		PlasticHardeningMaterial *theMatX = getTclPlasticMaterial(interp, argv[4], theBuilder);
		if(theMatX == 0)
			return TCL_ERROR;

		PlasticHardeningMaterial *theMatY = getTclPlasticMaterial(interp, argv[5], theBuilder);
		if(theMatY == 0)
			return TCL_ERROR;

double dir;
		if (Tcl_GetDouble(interp, argv[6], &dir) != TCL_OK)
			return TCL_ERROR;

		// Parsing was successful, allocate the material
		theModel = new Kinematic2D01(tag, minIsoFactor, *theMatX, *theMatY, dir);


return addTclYS_Evolution(theBuilder, theModel);
}


int TclIsotropic2D01Command(ClientData clienData, Tcl_Interp *interp, int argc,
				 TCL_Char **argv, TclModelBuilder *theBuilder)
{
YS_Evolution *theModel = 0;

		int tag;
		double minIsoFactor;

		if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK)
			return TCL_ERROR;

		if (Tcl_GetDouble(interp, argv[3], &minIsoFactor) != TCL_OK)
			return TCL_ERROR;

		PlasticHardeningMaterial *theMatX = getTclPlasticMaterial(interp, argv[4], theBuilder);
		if(theMatX == 0)
			return TCL_ERROR;

		PlasticHardeningMaterial *theMatY = getTclPlasticMaterial(interp, argv[5], theBuilder);
		if(theMatY == 0)
			return TCL_ERROR;

		// Parsing was successful, allocate the material
		theModel = new Isotropic2D01(tag, minIsoFactor, *theMatX, *theMatY);

return addTclYS_Evolution(theBuilder, theModel);
}


int TclPeakOriented2D01Command(ClientData clienData, Tcl_Interp *interp, int argc,
				 TCL_Char **argv, TclModelBuilder *theBuilder)
{
YS_Evolution *theModel = 0;
		int tag;
		double minIsoFactor;


		if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK)
			return TCL_ERROR;

		if (Tcl_GetDouble(interp, argv[3], &minIsoFactor) != TCL_OK)
			return TCL_ERROR;

		PlasticHardeningMaterial *theMatX = getTclPlasticMaterial(interp, argv[4], theBuilder);
		if(theMatX == 0)
			return TCL_ERROR;

		PlasticHardeningMaterial *theMatY = getTclPlasticMaterial(interp, argv[5], theBuilder);
		if(theMatY == 0)
			return TCL_ERROR;

		// Parsing was successful, allocate the material
		theModel = new PeakOriented2D01(tag, minIsoFactor, *theMatX, *theMatY);


return addTclYS_Evolution(theBuilder, theModel);
}



int TclCombinedIsoKin2D01Command(ClientData clienData, Tcl_Interp *interp, int argc,
				 TCL_Char **argv, TclModelBuilder *theBuilder)
{
YS_Evolution *theModel = 0;

		int tag;
		double minIsoFactor, iso_ratio, kin_ratio, shr_iso_ratio, shr_kin_ratio;
		int deformable;
		bool deform = false;


		if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK)
			return TCL_ERROR;
		if (Tcl_GetDouble(interp, argv[3], &iso_ratio) != TCL_OK)
			return TCL_ERROR;
		if (Tcl_GetDouble(interp, argv[4], &kin_ratio) != TCL_OK)
			return TCL_ERROR;
		if (Tcl_GetDouble(interp, argv[5], &shr_iso_ratio) != TCL_OK)
			return TCL_ERROR;
		if (Tcl_GetDouble(interp, argv[6], &shr_kin_ratio) != TCL_OK)
			return TCL_ERROR;
		if (Tcl_GetDouble(interp, argv[7], &minIsoFactor) != TCL_OK)
			return TCL_ERROR;

		PlasticHardeningMaterial *kpx_pos = getTclPlasticMaterial(interp, argv[8], theBuilder);
		if(kpx_pos == 0)
			return TCL_ERROR;

		PlasticHardeningMaterial *kpx_neg = getTclPlasticMaterial(interp, argv[9], theBuilder);
		if(kpx_neg == 0)
			return TCL_ERROR;

		PlasticHardeningMaterial *kpy_pos = getTclPlasticMaterial(interp, argv[10], theBuilder);
		if(kpx_pos == 0)
			return TCL_ERROR;

		PlasticHardeningMaterial *kpy_neg = getTclPlasticMaterial(interp, argv[11], theBuilder);
		if(kpx_neg == 0)
			return TCL_ERROR;

		if (Tcl_GetInt(interp, argv[12], &deformable) != TCL_OK)
			return TCL_ERROR;

double dir;
		if (Tcl_GetDouble(interp, argv[13], &dir) != TCL_OK)
			return TCL_ERROR;

		if (deformable == 1)
			deform = true;

		// Parsing was successful, allocate the material
		theModel = new CombinedIsoKin2D01(tag, iso_ratio, kin_ratio,
							shr_iso_ratio, shr_kin_ratio, minIsoFactor,
							*kpx_pos, *kpx_neg, *kpy_pos, *kpy_neg, deform, dir);

return addTclYS_Evolution(theBuilder, theModel);
}


////////////////////////////////////////////////////////////////////////////////////////
int TclKinematic2D02Command(ClientData clienData, Tcl_Interp *interp, int argc,
				 TCL_Char **argv, TclModelBuilder *theBuilder)
{
YS_Evolution *theModel = 0;
int tag;
double minIsoFactor;

	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK)
		return TCL_ERROR;

	if (Tcl_GetDouble(interp, argv[3], &minIsoFactor) != TCL_OK)
		return TCL_ERROR;

	YieldSurface_BC *ys = getTclYieldSurface_BC(interp, argv[4], theBuilder);
	if(ys==0)
		return TCL_ERROR;

	PlasticHardeningMaterial *theMatX = getTclPlasticMaterial(interp, argv[5], theBuilder);
	if(theMatX == 0)
		return TCL_ERROR;

	PlasticHardeningMaterial *theMatY = getTclPlasticMaterial(interp, argv[6], theBuilder);
	if(theMatY == 0)
		return TCL_ERROR;

int algo;
double resfact, appfact, dir;

	if (Tcl_GetInt(interp, argv[7], &algo) != TCL_OK)
		return TCL_ERROR;

	if (Tcl_GetDouble(interp, argv[8], &resfact) != TCL_OK)
		return TCL_ERROR;

	if (Tcl_GetDouble(interp, argv[9], &appfact) != TCL_OK)
		return TCL_ERROR;

	if (Tcl_GetDouble(interp, argv[10], &dir) != TCL_OK)
		return TCL_ERROR;

	// Parsing was successful, allocate the material
	theModel = new Kinematic2D02(tag, minIsoFactor, *ys,
	                               *theMatX, *theMatY, algo, resfact, appfact, dir);


return addTclYS_Evolution(theBuilder, theModel);
}


int TclPeakOriented2D02Command(ClientData clienData, Tcl_Interp *interp, int argc,
				 TCL_Char **argv, TclModelBuilder *theBuilder)
{
YS_Evolution *theModel = 0;

int tag;
double minIsoFactor;

	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK)
		return TCL_ERROR;

	if (Tcl_GetDouble(interp, argv[3], &minIsoFactor) != TCL_OK)
		return TCL_ERROR;

	YieldSurface_BC *ys = getTclYieldSurface_BC(interp, argv[4], theBuilder);
	if(ys==0)
		return TCL_ERROR;

	PlasticHardeningMaterial *kinX = getTclPlasticMaterial(interp, argv[5], theBuilder);
	if(kinX == 0)
		return TCL_ERROR;

	PlasticHardeningMaterial *kinY = getTclPlasticMaterial(interp, argv[6], theBuilder);
	if(kinY == 0)
		return TCL_ERROR;

	PlasticHardeningMaterial *isoX = getTclPlasticMaterial(interp, argv[7], theBuilder);
	if(isoX == 0)
		return TCL_ERROR;

	PlasticHardeningMaterial *isoY = getTclPlasticMaterial(interp, argv[8], theBuilder);
	if(isoY == 0)
		return TCL_ERROR;
int algo;
	if (Tcl_GetInt(interp, argv[9], &algo) != TCL_OK)
		return TCL_ERROR;

	// Parsing was successful, allocate the material
	theModel = new PeakOriented2D02(tag, minIsoFactor, *ys, *kinX, *kinY, *isoX, *isoY, algo);

return addTclYS_Evolution(theBuilder, theModel);
}


int TclCombinedIsoKin2D02Command(ClientData clienData, Tcl_Interp *interp, int argc,
				 TCL_Char **argv, TclModelBuilder *theBuilder)
{
YS_Evolution *theModel = 0;
int tag, deformable;
bool deform = false;
double minIsoFactor, isoRatio, kinRatio;

	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK)
		return TCL_ERROR;

	if (Tcl_GetDouble(interp, argv[3], &minIsoFactor) != TCL_OK)
		return TCL_ERROR;

	if (Tcl_GetDouble(interp, argv[4], &isoRatio) != TCL_OK)
		return TCL_ERROR;

	if (Tcl_GetDouble(interp, argv[5], &kinRatio) != TCL_OK)
		return TCL_ERROR;

	YieldSurface_BC *ys = getTclYieldSurface_BC(interp, argv[6], theBuilder);
	if(ys==0)
		return TCL_ERROR;

	PlasticHardeningMaterial *kinX = getTclPlasticMaterial(interp, argv[7], theBuilder);
	if(kinX == 0)
		return TCL_ERROR;

	PlasticHardeningMaterial *kinY = getTclPlasticMaterial(interp, argv[8], theBuilder);
	if(kinY == 0)
		return TCL_ERROR;

	PlasticHardeningMaterial *isoXPos = getTclPlasticMaterial(interp, argv[9], theBuilder);
	if(isoXPos == 0)
		return TCL_ERROR;

	PlasticHardeningMaterial *isoXNeg = getTclPlasticMaterial(interp, argv[10], theBuilder);
	if(isoXNeg == 0)
		return TCL_ERROR;

	PlasticHardeningMaterial *isoYPos = getTclPlasticMaterial(interp, argv[11], theBuilder);
	if(isoYPos == 0)
		return TCL_ERROR;

	PlasticHardeningMaterial *isoYNeg = getTclPlasticMaterial(interp, argv[12], theBuilder);
	if(isoYNeg == 0)
		return TCL_ERROR;

	if (Tcl_GetInt(interp, argv[13], &deformable) != TCL_OK)
		return TCL_ERROR;

	if(deformable == 1)
		deform = true;
int algo;
	if (Tcl_GetInt(interp, argv[14], &algo) != TCL_OK)
		return TCL_ERROR;

double resfact, appfact, dir;

	if (Tcl_GetDouble(interp, argv[15], &resfact) != TCL_OK)
		return TCL_ERROR;

	if (Tcl_GetDouble(interp, argv[16], &appfact) != TCL_OK)
		return TCL_ERROR;

	if (Tcl_GetDouble(interp, argv[17], &dir) != TCL_OK)
		return TCL_ERROR;

		
	// Parsing was successful, allocate the material
	theModel = new CombinedIsoKin2D02(tag, minIsoFactor, isoRatio, kinRatio, *ys, *kinX, *kinY,
										*isoXPos, *isoXNeg, *isoYPos, *isoYNeg,
										deform, algo, resfact, appfact, dir);

return addTclYS_Evolution(theBuilder, theModel);
}




int
TclModelBuilderYS_EvolutionModelCommand (ClientData clientData, Tcl_Interp *interp, int argc,
				 TCL_Char **argv, TclModelBuilder *theBuilder)
{
    if (strcmp(argv[1],"null") == 0)
	{
		return TclNullEvolutionCommand(clientData, interp, argc, argv, theBuilder);
	}
	else if (strcmp(argv[1],"kinematic2D01") == 0)
	{
		return TclKinematic2D01Command(clientData, interp, argc, argv, theBuilder);
	}
	else if (strcmp(argv[1],"isotropic2D01") == 0)
	{
		return TclIsotropic2D01Command(clientData, interp, argc, argv, theBuilder);
	}
	else if (strcmp(argv[1],"peakOriented2D01") == 0)
	{
		return TclPeakOriented2D01Command(clientData, interp, argc, argv, theBuilder);
	}
	else if (strcmp(argv[1],"combinedIsoKin2D01") == 0)
	{
		return TclCombinedIsoKin2D01Command(clientData, interp, argc, argv, theBuilder);
	}

	else if (strcmp(argv[1],"kinematic2D02") == 0)
	{
		return TclKinematic2D02Command(clientData, interp, argc, argv, theBuilder);
	}
	else if (strcmp(argv[1],"peakOriented2D02") == 0)
	{
		return TclPeakOriented2D02Command(clientData, interp, argc, argv, theBuilder);
	}
	else if (strcmp(argv[1],"combinedIsoKin2D02") == 0)
	{
		return TclCombinedIsoKin2D02Command(clientData, interp, argc, argv, theBuilder);
	}
	else
	{
		opserr << "Unknown YS_Evolution type: " << argv[1] << endln;
		return TCL_ERROR;
	}


}

///







