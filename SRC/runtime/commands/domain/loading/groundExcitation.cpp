//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation
//
//===----------------------------------------------------------------------===//
//
//
#if 0
int
groundExcitation(ClientData clientData, Tcl_Interp *interp, int argc,
                  TCL_Char ** const argv)
{
  G3_Runtime *rt = G3_getRuntime(interp);
  Domain* the_domain = G3_getDomain(rt);

  // make sure at least one other argument to contain integrator
  if (argc < 2) {
      opserr << "WARNING need to specify the commitTag \n";
      return TCL_ERROR;
  }

  if (strcmp(argv[1],"Single") == 0) {
      if (argc < 4) {
        opserr << "WARNING quake single dof motion\n";
        return TCL_ERROR;
      }

      int dof;
      if (Tcl_GetInt(interp, argv[2], &dof) != TCL_OK)
          return TCL_ERROR;

      // read in the ground motion
      GroundMotion *theMotion;
      if (strcmp(argv[3],"ElCentro") == 0) {
          double fact = 1.0;
          if (argc == 5) {
              if (Tcl_GetDouble(interp, argv[4], &fact) != TCL_OK)
                  return TCL_ERROR;
          }
          theMotion = new ElCentroGroundMotion(fact);
      } else {
          opserr << "WARNING quake Single motion - no motion type exists \n";
          return TCL_ERROR;
      }

      Load *theLoad = new SingleExcitation(*theMotion, dof, nextTag++);
      the_domain->addOtherLoad(theLoad);
      return TCL_OK;
  }

  else {
    opserr << "WARNING No quake type exists \n";
    return TCL_ERROR;
  }
}
#endif

