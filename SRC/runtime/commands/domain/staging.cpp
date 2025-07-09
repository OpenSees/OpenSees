//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation
//
//===----------------------------------------------------------------------===//
// 
int
elementActivate(ClientData clientData, Tcl_Interp *interp, int argc,
                TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  Domain* the_domain = (Domain*)clientData; 

  int eleTag;
  int argLoc = 1;
  int Nelements = argc;
  ID activate_us(0, Nelements);

  while (argLoc < argc && Tcl_GetInt(interp, argv[argLoc], &eleTag) == TCL_OK) {
    activate_us.insert(eleTag);
    ++argLoc;
  }

  the_domain->activateElements(activate_us);

  return TCL_OK;
}

int
elementDeactivate(ClientData clientData, Tcl_Interp *interp, int argc,
                  TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  Domain* the_domain = (Domain*)clientData; 

  int eleTag;
  int argLoc = 1;
  int Nelements = argc;
  ID deactivate_us(0, Nelements);

  while (argLoc < argc && Tcl_GetInt(interp, argv[argLoc], &eleTag) == TCL_OK) {
    deactivate_us.insert(eleTag);
    ++argLoc;
  }

  the_domain->deactivateElements(deactivate_us);
  return TCL_OK;
}
