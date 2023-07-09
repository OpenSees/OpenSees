#include <stdlib.h>
#include <string.h>
#include <OPS_Stream.h>

#include <Domain.h>
#include <Node.h>
#include <Matrix.h>

#include <CyclicModel.h>
#include <LinearCyclic.h>
#include <BilinearCyclic.h>
#include <QuadraticCyclic.h>
#include <TclModelBuilder.h>


int TclModelBuilder_addLinearCylic(ClientData clientData, Tcl_Interp *interp,
				 int argc, TCL_Char **argv,
				 TclModelBuilder *theBuilder)
{
int tag;
	
	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK)
	{
		opserr << "WARNING invalid CyclicModel tag" << endln;
		return TCL_ERROR;
	}

	CyclicModel *cModel = new LinearCyclic(tag);
	if(!cModel)
	{
		opserr << "TclModelBuilder_addLinearCycylic - could not allocate memory\n";
		return TCL_ERROR;
	}
	if (theBuilder->addCyclicModel(*cModel) < 0)
	{
		opserr << "WARNING TclElmtBuilder - could not add cycModel to domain ";
		opserr << tag << endln;
		opserr << "\a";
		return TCL_ERROR;
	}	
	
	return TCL_OK;
}

int TclModelBuilder_addBilinearCyclic(ClientData clientData, Tcl_Interp *interp,
				 int argc, TCL_Char **argv,
				 TclModelBuilder *theBuilder)
{
int tag;
double wt;
	
	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK)
	{
		opserr << "WARNING invalid CyclicModel tag" << endln;
		return TCL_ERROR;
	}
	if (Tcl_GetDouble(interp, argv[3], &wt) != TCL_OK)
	{
		opserr << "WARNING invalid arg[3]" << endln;
		return TCL_ERROR;
	}
	
	CyclicModel *cModel = new BilinearCyclic(tag, wt);
	if (theBuilder->addCyclicModel(*cModel) <0)
	{
		opserr << "WARNING TclElmtBuilder - could not add cycModel to domain ";
		opserr << tag << endln;
		opserr << "\a";
		return TCL_ERROR;
	}	
	
	return TCL_OK;
}


int TclModelBuilder_addQuadraticCyclic(ClientData clientData, Tcl_Interp *interp,
				 int argc, TCL_Char **argv,
				 TclModelBuilder *theBuilder)
{
int tag;
double wt, qy;
	
	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK)
	{
		opserr << "WARNING invalid CyclicModel tag" << endln;
		return TCL_ERROR;
	}
	if (Tcl_GetDouble(interp, argv[3], &wt) != TCL_OK)
	{
		opserr << "WARNING invalid arg[3]" << endln;
		return TCL_ERROR;
	}
	if (Tcl_GetDouble(interp, argv[4], &qy) != TCL_OK)
	{
		opserr << "WARNING invalid arg[4]" << endln;
		return TCL_ERROR;
	}

	
	CyclicModel *cModel = new QuadraticCyclic(tag, wt, qy);
	if (theBuilder->addCyclicModel(*cModel) < 0)
	{
		opserr << "WARNING TclElmtBuilder - could not add cycModel to domain ";
		opserr << tag << endln;
		opserr << "\a";
		return TCL_ERROR;
	}	
	return TCL_OK;
}


/*******************************************************************************************/
int
TclModelBuilderCyclicModelCommand (ClientData clientData, Tcl_Interp *interp, int argc,
					               TCL_Char **argv, TclModelBuilder *theTclBuilder)
{

  if (strcmp(argv[1],"linear") == 0) {
	  int result = TclModelBuilder_addLinearCylic(clientData, interp, argc, argv,
						 theTclBuilder);
    return result;
  }
  else if (strcmp(argv[1],"bilinear") == 0) {
	  int result = TclModelBuilder_addBilinearCyclic(clientData, interp, argc, argv,
						 theTclBuilder);
    return result;
  }
  
   else if (strcmp(argv[1],"quadratic") == 0) {
	  int result = TclModelBuilder_addQuadraticCyclic(clientData, interp, argc, argv,
						  theTclBuilder);
    return result;
  }

  else

	return TCL_ERROR;

}
