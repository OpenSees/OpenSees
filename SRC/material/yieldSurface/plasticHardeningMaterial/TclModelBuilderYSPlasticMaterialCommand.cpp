#include <TclModelBuilder.h>
#include <string.h>
#include <Vector.h>

#include "MultiLinearKp.h"
#include "ExponReducing.h"
#include "NullPlasticMaterial.h"

int TclMultiLinearCommand (ClientData clienData, Tcl_Interp *interp, int argc,
					TCL_Char **argv, TclModelBuilder *theTclBuilder)
{
    // Pointer to a uniaxial material that will be added to the model builder
    PlasticHardeningMaterial *theMaterial = 0;

	if (strcmp(argv[1],"multiLinearKp") == 0)
	{
		int tag;
		int numPoints = (argc - 3)/2;

		if(numPoints < 2)
		{
		 	opserr << "WARNING invalid uniaxialMaterial MultilinearUniaxial tag" << endln;
		 	opserr << "Minimum of 2 points are required\n";
			return TCL_ERROR;
		}
		
		if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK)
		{
	    	opserr << "WARNING invalid uniaxialMaterial MultilinearUniaxial tag" << endln;
	    	return TCL_ERROR;
		}

		Vector defo(numPoints);
		Vector kp(numPoints);
		
		double temp;
        int indx=3, j1, j2;

        for (j1 = 0; j1 < numPoints; j1++)
        {
			if (Tcl_GetDouble(interp, argv[indx], &temp) != TCL_OK)
			{
				opserr << "WARNING invalid input, data " << temp << '\n';
				opserr << "MultilinearUniaxial material: " << tag << endln;
				return TCL_ERROR;
			}
			
			defo(j1) = temp;
			indx++;
		}

		for (j2 = 0; j2 < numPoints; j2++)
        {
			if (Tcl_GetDouble(interp, argv[indx], &temp) != TCL_OK)
			{
				opserr << "WARNING invalid input, data " << temp << '\n';
				opserr << "MultilinearUniaxial material: " << tag << endln;
				return TCL_ERROR;
			}
			
			kp(j2) = temp;
			indx++;
		}
		
		// Parsing was successful, allocate the material
		theMaterial = new MultiLinearKp(tag, defo, kp);
    }


    // Ensure we have created the Material, out of memory if got here and no material
    if (theMaterial == 0) {
	opserr << "WARNING: ran out of memory creating uniaxialMaterial\n";
	opserr << argv[1] << endln;
	return TCL_ERROR;
    }

    // Now add the material to the modelBuilder
    if (theTclBuilder->addPlasticMaterial(*theMaterial) < 0)
	{
		opserr << "WARNING could not add uniaxialMaterial to the domain\n";
		opserr << *theMaterial << endln;
		delete theMaterial; // invoke the material objects destructor, otherwise mem leak
		return TCL_ERROR;
    }


	return TCL_OK;
}


// QuadrReducing(int tag, double kp0, double kp_half);
/*int TclQuadrReducingCommand(ClientData clienData, Tcl_Interp *interp, int argc,
					char **argv, TclModelBuilder *theTclBuilder)
{
    // Pointer to a uniaxial material that will be added to the model builder
    PlasticHardeningMaterial *theMaterial = 0;

	int tag;
	double kp_0;
	double kp_half;

	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK)
	{
	    	opserr << "WARNING invalid  PlaticHardening quadrReducing tag" << endln;
	    	return TCL_ERROR;
	}
	if (Tcl_GetDouble(interp, argv[3], &kp_0) != TCL_OK)
	{
	    	opserr << "WARNING invalid  PlaticHardening quadrReducing kp_0" << endln;
	    	return TCL_ERROR;
	}
	if (Tcl_GetDouble(interp, argv[4], &kp_half) != TCL_OK)
	{
	    	opserr << "WARNING invalid  PlaticHardening quadrReducing kp_half" << endln;
	    	return TCL_ERROR;
	}

	theMaterial = new QuadrReducing(tag, kp_0, kp_half);
    if (theTclBuilder->addPlasticMaterial(*theMaterial) < 0)
	{
		opserr << "WARNING could not add uniaxialMaterial to the domain\n";
		opserr << *theMaterial << endln;
		delete theMaterial; // invoke the material objects destructor, otherwise mem leak
		return TCL_ERROR;
    }

	return TCL_OK;
}
*/

int TclExponReducingCommand(ClientData clienData, Tcl_Interp *interp, int argc,
					TCL_Char **argv, TclModelBuilder *theTclBuilder)
{
	if(argc < 5)
	{
		opserr << "TclExponReducingCommand - argc != 5 \n";
		return TCL_ERROR;
	}
	
	PlasticHardeningMaterial  *theMaterial = 0;

int tag;
double arg1, arg2, arg3;


//plasticMaterial exponReducing (int tag, double kp0, double alfa); //5
//plasticMaterial exponReducing (int tag, double kp0, double x0, double tol); //6
	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK)
	{
	    	opserr << "WARNING invalid  PlaticHardening exponReducing tag" << endln;
	    	return TCL_ERROR;
	}
	
	if (Tcl_GetDouble(interp, argv[3], &arg1) != TCL_OK)
	{
	    	opserr << "WARNING invalid double PlaticHardening exponReducing" << endln;
	    	return TCL_ERROR;
	}

	if (Tcl_GetDouble(interp, argv[4], &arg2) != TCL_OK)
	{
	    	opserr << "WARNING invalid double PlaticHardening exponReducing" << endln;
	    	return TCL_ERROR;
	}

	if(argc == 6)
	{
		if (Tcl_GetDouble(interp, argv[5], &arg3) != TCL_OK)
		{
				opserr << "WARNING invalid double PlaticHardening exponReducing" << endln;
				return TCL_ERROR;
		}

		theMaterial = new ExponReducing(tag, arg1, arg2, arg3);
//		opserr << "factor = " << arg3 << endln;
	}
	else
		theMaterial = new ExponReducing(tag, arg1, arg2);
	
	if (theTclBuilder->addPlasticMaterial(*theMaterial) < 0)
	{
		opserr << "WARNING could not add uniaxialMaterial to the domain\n";
		opserr << *theMaterial << endln;
		delete theMaterial; // invoke the material objects destructor, otherwise mem leak
		return TCL_ERROR;
    }


	return TCL_OK;
}


int TclNullPlasticMaterialCommand(ClientData clienData, Tcl_Interp *interp, int argc,
					TCL_Char **argv, TclModelBuilder *theTclBuilder)
{
    PlasticHardeningMaterial *theMaterial = 0;

	int tag;

	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK)
	{
	    	opserr << "WARNING invalid  PlaticHardening quadrReducing tag" << endln;
	    	return TCL_ERROR;
	}

	theMaterial = new NullPlasticMaterial(tag);
    if (theTclBuilder->addPlasticMaterial(*theMaterial) < 0)
	{
		opserr << "WARNING could not add uniaxialMaterial to the domain\n";
		opserr << *theMaterial << endln;
		delete theMaterial; // invoke the material objects destructor, otherwise mem leak
		return TCL_ERROR;
    }

	return TCL_OK;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////





int
TclModelBuilderPlasticMaterialCommand (ClientData clientData, Tcl_Interp *interp, int argc,
				 TCL_Char **argv, TclModelBuilder *theTclBuilder)
{

	if (strcmp(argv[1],"multiLinearKp") == 0)
	{
		return TclMultiLinearCommand(clientData, interp, argc, argv, theTclBuilder);
	}
	/*else if(strcmp(argv[1],"quadrReducing") == 0)
	{
		return TclQuadrReducingCommand(clientData, interp, argc, argv, theTclBuilder);
	}*/
	else if(strcmp(argv[1],"exponReducing") == 0)
	{
		return TclExponReducingCommand(clientData, interp, argc, argv, theTclBuilder);
	}
	else if(strcmp(argv[1],"null")==0)
	{
		return TclNullPlasticMaterialCommand(clientData, interp, argc, argv, theTclBuilder);
	}
	else
	{
		opserr << "Unknown PlasticMaterial: \nValid types: null, multiLinearKp, "
		<< "quadrReducing, exponReducing \n";

		return TCL_ERROR;
	}

}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////





