/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** With a lot of additions by                                         **
**   Boris Jeremic    (jeremic@ucdavis.edu)                           **
**   Zaohui Yang      (zhyang@ucdavis.edu)                            **
**   Zhao Cheng       (zcheng@ucdavis.edu)                            **
**                                                                    **
** ****************************************************************** */
//
#include <tcl.h>
#include <string.h>
#include <Logging.h>
#include <Parsing.h>
#include <BasicModelBuilder.h>

#include <PressureDependentElastic3D.h>

#include <PlaneStressMaterial.h>
#include <PlateFiberMaterial.h>
#include <BeamFiberMaterial.h>

#include <PressureIndependMultiYield.h>
#include <PressureDependMultiYield.h>
#include <PressureDependMultiYield02.h>
#include <FluidSolidPorousMaterial.h>


Tcl_CmdProc TclCommand_newPlasticMaterial;
Tcl_CmdProc TclCommand_newElasticMaterial;



int
TclCommand_addMaterial(ClientData clientData, Tcl_Interp* interp, 
                        Tcl_Size argc, TCL_Char** const argv)
{
  static
  std::unordered_map<std::string, Tcl_CmdProc*> MaterialLibrary = {
    {"ElasticIsotropic",          TclCommand_newElasticMaterial},
    {"Elastic",                   TclCommand_newElasticMaterial},
    {"Isotropic",                 TclCommand_newElasticMaterial},
    {"J2",                        TclCommand_newPlasticMaterial},
    {"J2Simplified",              TclCommand_newPlasticMaterial},
    {"J2BeamFiber",               TclCommand_newPlasticMaterial},
  };


  if (argc < 2) {
    opserr << OpenSees::PromptValueError
           << "missing argument type"
           << "\n";
    return TCL_ERROR;
  }


  auto cmd = MaterialLibrary.find(std::string(argv[1]));
  if (cmd != MaterialLibrary.end())
    return (*cmd->second)(clientData, interp, argc, &argv[0]);

  return TCL_ERROR;

}
