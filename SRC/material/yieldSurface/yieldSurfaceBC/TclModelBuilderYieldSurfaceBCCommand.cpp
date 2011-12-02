#include "YieldSurface_BC.h"
#include <TclModelBuilder.h>
#include <string.h>
#include <Vector.h>

#include "NullYS2D.h"
#include "Attalla2D.h"
#include "Orbison2D.h"
#include "Hajjar2D.h"
#include "ElTawil2D.h"
#include "ElTawil2DUnSym.h"


static void printCommand(int argc, TCL_Char **argv)
{
    opserr << "Input command: ";
    for (int i=0; i<argc; i++)
	opserr << argv[i] << " ";
    opserr << endln;
}

int
TclModelBuilderYieldSurface_BCCommand (ClientData clienData, Tcl_Interp *interp, int argc,
				       TCL_Char **argv, TclModelBuilder *theBuilder)
{
    // Make sure there is a minimum number of arguments
    if (argc < 3) {
	opserr << "WARNING insufficient number of uniaxial material arguments\n";
	opserr << "Want: yieldSurfaceBC type? tag? <specific material args>" << endln;
	return TCL_ERROR;
    }

    // Pointer to a ys that will be added to the model builder
    YieldSurface_BC *theYS = 0;

    if(strcmp(argv[1],"null") == 0)
    {
		if(argc < 4)
		{
			opserr << "WARNING invalid number of arguments\n";
			printCommand(argc,argv);
			opserr << "Want: yieldSurfaceBC null tag? dimensions?" << endln;
			return TCL_ERROR;
		}
		int tag, dim;
		if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK)
		{
			opserr << "WARNING invalid yieldSurfaceBC null tag" << endln;
			return TCL_ERROR;
		}
		if (Tcl_GetInt(interp, argv[3], &dim) != TCL_OK)
		{
			opserr << "WARNING invalid yieldSurfaceBC null dimensions" << endln;
			return TCL_ERROR;
		}

		switch (dim)
		{
			// case 1: 1D YS
			case 2: theYS = new NullYS2D(tag); break;
			// case 3: 3D YS
            default:
				opserr << "incorrect dimension for null ys\n";
				return TCL_ERROR;
		}

	}
    else if (strcmp(argv[1],"Orbison2D") == 0)
	{
		if (argc < 6)
		{
			opserr << "WARNING invalid number of arguments\n";
			printCommand(argc,argv);
			// Orbison2D(int tag, double xmax, double ymax, YS_HardeningModel &model);
			opserr << "Want: yieldSurfaceBC Orbison2D tag? xCap? yCap? ys_model_tag?" << endln;
			return TCL_ERROR;
		}

		int tag;
		double xCap, yCap;
		// int matID1, matID2;
		int modelID;
//		double isoRatio;

		if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK)
		{
			opserr << "WARNING invalid yieldSurfaceBC Orbison2D tag" << endln;
			return TCL_ERROR;
		}

		if (Tcl_GetDouble(interp, argv[3], &xCap) != TCL_OK)
		{
			opserr << "WARNING invalid xCap\n";
			opserr << "yieldSurfaceBC Orbison2D tag: " << tag << endln;
			return TCL_ERROR;
		}

		if (Tcl_GetDouble(interp, argv[4], &yCap) != TCL_OK)
		{
			opserr << "WARNING invalid yCap\n";
			opserr << "yieldSurfaceBC Orbison2D tag: " << tag << endln;
			return TCL_ERROR;
		}

		if (Tcl_GetInt(interp, argv[5], &modelID) != TCL_OK)
		{
			opserr << "WARNING invalid yieldSurfaceBC Orbison2D matID1" << modelID << endln;
			return TCL_ERROR;
		}

		YS_Evolution *theModel = theBuilder->getYS_EvolutionModel(modelID);
		if(theModel == 0)
		{
			opserr << "WARNING yieldSurfaceBC Orbison2D no ys_model exixts with tag: " << modelID << endln;
			return TCL_ERROR;
		}

		// Parsing was successful, allocate the material
		theYS = new Orbison2D(tag, xCap, yCap, *theModel);
    }

    else if (strcmp(argv[1],"ElTawil2D") == 0)
	{
		if (argc < 7)
		{
			opserr << "WARNING invalid number of arguments\n";
			printCommand(argc,argv);
			// Orbison2D(int tag, double xmax, double ymax, YS_HardeningModel &model);
			opserr << "Want: yieldSurfaceBC ElTawil2D tag? xCap? yCap? ys_model_tag?" << endln;
			return TCL_ERROR;
		}

		int tag;
		double xBal, yBal;
		double yPos, yNeg;

		int modelID;

		if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK)
		{
			opserr << "WARNING invalid yieldSurfaceBC ElTawil2D tag" << endln;
			return TCL_ERROR;
		}

		if (Tcl_GetDouble(interp, argv[3], &xBal) != TCL_OK)
		{
			opserr << "WARNING invalid xBal\n";
			opserr << "yieldSurfaceBC ElTawil2D tag: " << tag << endln;
			return TCL_ERROR;
		}

		if (Tcl_GetDouble(interp, argv[4], &yBal) != TCL_OK)
		{
			opserr << "WARNING invalid yBal\n";
			opserr << "yieldSurfaceBC ElTawil2D tag: " << tag << endln;
			return TCL_ERROR;
		}
		
		if (Tcl_GetDouble(interp, argv[5], &yPos) != TCL_OK)
		{
			opserr << "WARNING invalid xPos\n";
			opserr << "yieldSurfaceBC ElTawil2D tag: " << tag << endln;
			return TCL_ERROR;
		}

		if (Tcl_GetDouble(interp, argv[6], &yNeg) != TCL_OK)
		{
			opserr << "WARNING invalid yNeg\n";
			opserr << "yieldSurfaceBC ElTawil2D tag: " << tag << endln;
			return TCL_ERROR;
		}
		

		if (Tcl_GetInt(interp, argv[7], &modelID) != TCL_OK)
		{
			opserr << "WARNING invalid yieldSurfaceBC ElTawil2D matID1" << modelID << endln;
			return TCL_ERROR;
		}

		YS_Evolution *theModel = theBuilder->getYS_EvolutionModel(modelID);
		if(theModel == 0)
		{
			opserr << "WARNING yieldSurfaceBC ElTawil2D no ys_model exixts with tag: " << modelID << endln;
			return TCL_ERROR;
		}

		// Parsing was successful, allocate the material
		theYS = new ElTawil2D(tag, xBal, yBal, yPos, yNeg, *theModel);
    }
	
	    else if (strcmp(argv[1],"ElTawil2DUnSym") == 0)
	{
		if (argc < 9)
		{
			opserr << "WARNING invalid number of arguments\n";
			printCommand(argc,argv);
			// Orbison2D(int tag, double xmax, double ymax, YS_HardeningModel &model);
			opserr << "Want: yieldSurfaceBC ElTawil2DUnSym tag? xPosBal? yPosBal? "
			     << "xNegBal? yPos? yNeg? ys_model_tag?" << endln;
			return TCL_ERROR;
		}

		int tag;
		double xPosBal, yPosBal;
		double xNegBal, yNegBal;
		double yPos, yNeg;
		
		int modelID;

		if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK)
		{
			opserr << "WARNING invalid yieldSurfaceBC ElTawil2DUnSym tag" << endln;
			return TCL_ERROR;
		}

		if (Tcl_GetDouble(interp, argv[3], &xPosBal) != TCL_OK)
		{
			opserr << "WARNING invalid xPosBal\n";
			opserr << "yieldSurfaceBC ElTawil2DUnSym tag: " << tag << endln;
			return TCL_ERROR;
		}

		if (Tcl_GetDouble(interp, argv[4], &yPosBal) != TCL_OK)
		{
			opserr << "WARNING invalid yPosBal\n";
			opserr << "yieldSurfaceBC ElTawil2DUnSym tag: " << tag << endln;
			return TCL_ERROR;
		}
		if (Tcl_GetDouble(interp, argv[5], &xNegBal) != TCL_OK)
		{
			opserr << "WARNING invalid xNegBal\n";
			opserr << "yieldSurfaceBC ElTawil2DUnSym tag: " << tag << endln;
			return TCL_ERROR;
		}

		if (Tcl_GetDouble(interp, argv[6], &yNegBal) != TCL_OK)
		{
			opserr << "WARNING invalid yNegBal\n";
			opserr << "yieldSurfaceBC ElTawil2DUnSym tag: " << tag << endln;
			return TCL_ERROR;
		}
		
		if (Tcl_GetDouble(interp, argv[7], &yPos) != TCL_OK)
		{
			opserr << "WARNING invalid xPos\n";
			opserr << "yieldSurfaceBC ElTawil2DUnSym tag: " << tag << endln;
			return TCL_ERROR;
		}

		if (Tcl_GetDouble(interp, argv[8], &yNeg) != TCL_OK)
		{
			opserr << "WARNING invalid yNeg\n";
			opserr << "yieldSurfaceBC ElTawil2DUnSym tag: " << tag << endln;
			return TCL_ERROR;
		}
		

		if (Tcl_GetInt(interp, argv[9], &modelID) != TCL_OK)
		{
			opserr << "WARNING invalid yieldSurfaceBC ElTawil2DUnSym matID1" 
			     << modelID << endln;
			return TCL_ERROR;
		}

		YS_Evolution *theModel = theBuilder->getYS_EvolutionModel(modelID);
		if(theModel == 0)
		{
			opserr << "WARNING yieldSurfaceBC ElTawil2D no ys_model exixts with tag: " << modelID << endln;
			return TCL_ERROR;
		}

		// Parsing was successful, allocate the material
		theYS = new ElTawil2DUnSym(tag, xPosBal, yPosBal, xNegBal, yNegBal, yPos, yNeg, *theModel);
    }

	else if (strcmp(argv[1],"Attalla2D") == 0)
	{
// 	     Attalla2D( int tag, double xmax, double ymax, YS_HardeningModel &model,
// 					double x_offset=0, double y_offset=0,
// 					double a01=0.19,  double a02=0.54, double a03=-1.4,
// 					double a04=-1.64, double a05=2.21, double a06=2.10);

		if (argc < 6 || argc > 14)
		{
			opserr << "WARNING invalid number of arguments\n";
			printCommand(argc,argv);
			opserr << "Want: yieldSurfaceBC Attalla2D tag? xCap? yCap? matXTag? maxYTag? isoRatio? <..>" << endln;
			return TCL_ERROR;
		}

		int tag;
		double xCap, yCap;
		//int matID1, matID2;
		int modelID;
		//double isoRatio;
		//double x_offset = 0, y_offset = 0;
		Vector param(6);

		param[0] = 0.19;
		param[1] = 0.54;
		param[2] =-1.40;
		param[3] =-1.64;
		param[4] = 2.21;
		param[5] = 2.10;


		if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK)
		{
			opserr << "WARNING invalid yieldSurfaceBC Attalla2D tag" << endln;
			return TCL_ERROR;
		}

		if (Tcl_GetDouble(interp, argv[3], &xCap) != TCL_OK)
		{
			opserr << "WARNING invalid xCap\n";
			opserr << "yieldSurfaceBC Attalla2D tag: " << tag << endln;
			return TCL_ERROR;
		}

		if (Tcl_GetDouble(interp, argv[4], &yCap) != TCL_OK)
		{
			opserr << "WARNING invalid yCap\n";
			opserr << "yieldSurfaceBC Attalla2D tag: " << tag << endln;
			return TCL_ERROR;
		}

		if (Tcl_GetInt(interp, argv[5], &modelID) != TCL_OK)
		{
			opserr << "WARNING invalid yieldSurfaceBC Attalla2D modelID" << modelID << endln;
			return TCL_ERROR;
		}
		
		YS_Evolution *theModel = theBuilder->getYS_EvolutionModel(modelID);
		if(theModel == 0)
		{
			opserr << "WARNING yieldSurfaceBC Orbison2D no ys_model exixts with tag: " << modelID << endln;
			return TCL_ERROR;
		}



		if(argc > 6)
		{
			int count = 6;
			double temp;

			for(int i=0; i < 6; i++)
			{
				if (Tcl_GetDouble(interp, argv[count], &temp) != TCL_OK)
				{
					opserr << "WARNING invalid parameter " << i+1 << "\n";
					opserr << "yieldSurfaceBC Attalla2D tag: " << tag << endln;
					return TCL_ERROR;
				}
				param(i) = temp;
				count++;
			}

		}

		// Parsing was successful, allocate the material
		theYS = new Attalla2D(tag, xCap, yCap, *theModel,
							  param(0),param(1),param(2),param(3),param(4),param(5));

	}

	else if (strcmp(argv[1],"Hajjar2D") == 0)
	{
		if (argc < 9)
		{
			opserr << "WARNING invalid number of arguments\n";
			printCommand(argc,argv);
			opserr << "Want: yieldSurfaceBC Hajjar2D tag? ysModelTag? D? b? t? fc? fy?" << endln;
			return TCL_ERROR;
		}

		int tag;
		//int matID1, matID2;
		int modelID;
//		double isoRatio;
		double D, b, t, fc, fy;

		if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK)
		{
			opserr << "WARNING invalid yieldSurfaceBC Hajjar2D  tag" << endln;
			return TCL_ERROR;
		}

		if (Tcl_GetInt(interp, argv[3], &modelID) != TCL_OK)
		{
			opserr << "WARNING invalid yieldSurfaceBC Hajjar2D  matID1" << modelID << endln;
			return TCL_ERROR;
		}

		if (Tcl_GetDouble(interp, argv[4], &D) != TCL_OK)
		{
			opserr << "WARNING invalid D \n";
			opserr << "yieldSurfaceBC Hajjar2D  tag: " << tag << endln;
			return TCL_ERROR;
		}

		if (Tcl_GetDouble(interp, argv[5], &b) != TCL_OK)
		{
			opserr << "WARNING invalid b \n";
			opserr << "yieldSurfaceBC Hajjar2D  tag: " << tag << endln;
			return TCL_ERROR;
		}

		if (Tcl_GetDouble(interp, argv[6], &t) != TCL_OK)
		{
			opserr << "WARNING invalid t \n";
			opserr << "yieldSurfaceBC Hajjar2D  tag: " << tag << endln;
			return TCL_ERROR;
		}

		if (Tcl_GetDouble(interp, argv[7], &fc) != TCL_OK)
		{
			opserr << "WARNING invalid fc \n";
			opserr << "yieldSurfaceBC Hajjar2D  tag: " << tag << endln;
			return TCL_ERROR;
		}

		if (Tcl_GetDouble(interp, argv[8], &fy) != TCL_OK)
		{
			opserr << "WARNING invalid fy \n";
			opserr << "yieldSurfaceBC Hajjar2D  tag: " << tag << endln;
			return TCL_ERROR;
		}
		
		YS_Evolution *theModel = theBuilder->getYS_EvolutionModel(modelID);
		if(theModel == 0)
		{
			opserr << "WARNING yieldSurfaceBC Orbison2D no ys_model exixts with tag: " << modelID << endln;
			return TCL_ERROR;
		}

		// Parsing was successful, allocate the material
		theYS = new Hajjar2D(tag, *theModel, D, b, t, fc, fy);
    }

	else
	{
		opserr << "Warning - unknown yield surface type \n";
		printCommand(argc,argv);
	}

	///////////////////////////////////////////////////////////////
	// Now add the ys to the modelBuilder
	///////////////////////////////////////////////////////////////

	if (theBuilder->addYieldSurface_BC(*theYS) < 0)
	{
		opserr << "WARNING could not add YieldSurfaceBC to the domain\n";
		opserr << *theYS << endln;
		delete theYS; // invoke the material objects destructor, otherwise mem leak
		return TCL_ERROR;
	}


    return TCL_OK;
}







