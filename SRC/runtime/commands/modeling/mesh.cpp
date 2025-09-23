//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation    
//
//===----------------------------------------------------------------------===//
//
//
extern int OPS_LineMesh(Domain& domain, int ndm);
extern int OPS_TriMesh(Domain& domain);
extern int OPS_TriReMesh(Domain& domain, int ndf);

int
TclCommand_mesh(ClientData clientData, Tcl_Interp *interp, int argc,
                TCL_Char ** const argv)
{
  // ensure the destructor has not been called -
  assert(clientData != nullptr);

  // make sure corect number of arguments on command line
  if (argc < 2) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: mesh type? ...>\n";
    return TCL_ERROR;
  }

  OPS_ResetInput(clientData, interp, 2, argc, argv, theTclDomain,
                 theTclBuilder);

  // mesh type
  int res = 0;
  // if (strcmp(argv[1], "line") == 0) {
  // 	res = OPS_LineMesh(*theTclDomain,ndm);
  // } else if (strcmp(argv[1], "tri") == 0) {
  // 	res = OPS_TriMesh(*theTclDomain);
  // } else {
  // 	opserr<<"WARNING: mesh type "<<argv[1]<<" is unknown\n";
  // 	return TCL_ERROR;
  // }

  if (res < 0) {
    return TCL_ERROR;
  }

  return 0;
}

int
TclCommand_remesh(ClientData clientData, Tcl_Interp *interp, int argc,
                  TCL_Char ** const argv)
{
  // ensure the destructor has not been called -
  if (theTclBuilder == nullptr) {
    opserr << "WARNING builder has been destroyed" << endln;
    return TCL_ERROR;
  }


  // make sure corect number of arguments on command line
  if (argc < 2) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: mesh type? ...>\n";
    return TCL_ERROR;
  }

  OPS_ResetInput(clientData, interp, 2, argc, argv, theTclDomain,
                 theTclBuilder);

  // mesh type
  int res = 0;
  if (strcmp(argv[1], "line") == 0) {
  	//res = OPS_LineMesh(*theTclDomain,ndm);
  } else if (strcmp(argv[1], "tri") == 0) {
  	res = OPS_TriReMesh(*theTclDomain,ndf);
  } else {
  	opserr<<"WARNING: remesh type "<<argv[1]<<" is unknown\n";
  	return TCL_ERROR;
  }

  if (res < 0) {
    return TCL_ERROR;
  }

  return 0;
}
