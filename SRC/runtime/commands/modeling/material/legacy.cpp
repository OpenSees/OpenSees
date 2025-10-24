//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation    
//
//===----------------------------------------------------------------------===//
//
#include <tcl.h>
#include <string.h>
#include <assert.h>
#include <Logging.h>
#include <Parsing.h>
#include <BasicModelBuilder.h>

#include <ENTMaterial.h>             // MHS
#include <Elastic2Material.h>        // ZHY
#include <Concrete01WithSITC.h>      // Won Lee
#include <BarSlipMaterial.h>
#include <ShearPanelMaterial.h>      // NM


int
TclDispatch_LegacyUniaxials(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char**const argv)
{
  assert(clientData != nullptr);
  BasicModelBuilder *builder = static_cast<BasicModelBuilder*>(clientData);
  UniaxialMaterial  *theMaterial = nullptr;

  if (strcmp(argv[1], "Elastic2") == 0) {
    if (argc < 4 || argc > 5) {
      opserr << OpenSees::PromptValueError << "invalid number of arguments\n";
      opserr << "Want: uniaxialMaterial Elastic tag? E? <eta?>" << "\n";
      return TCL_ERROR;
    }

    int tag;
    double E;
    double eta = 0.0;

    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid uniaxialMaterial Elastic tag" << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[3], &E) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid E\n";
      opserr << "uniaxiaMaterial Elastic: " << tag << "\n";
      return TCL_ERROR;
    }

    if (argc == 5) {
      if (Tcl_GetDouble(interp, argv[4], &eta) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "invalid eta\n";
        opserr << "uniaxialMaterial Elastic: " << tag << "\n";
        return TCL_ERROR;
      }
    }

    // Parsing was successful, allocate the material
    theMaterial = new Elastic2Material(tag, E, eta);

  } else if (strcmp(argv[1], "ENT") == 0) {
    if (argc < 4) {
      opserr << OpenSees::PromptValueError << "invalid number of arguments\n";
      opserr << "Want: uniaxialMaterial ENT tag? E?" << "\n";
      return TCL_ERROR;
    }

    int tag;
    double E;

    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid uniaxialMaterial ENT tag" << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[3], &E) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid E\n";
      opserr << "uniaxiaMaterial ENT: " << tag << "\n";
      return TCL_ERROR;
    }

    // Parsing was successful, allocate the material
    theMaterial = new ENTMaterial(tag, E);

  }

  else if (strcmp(argv[1], "BarSlip") == 0) {
    if (argc != 17 && argc != 15) {
      opserr << OpenSees::PromptValueError << "insufficient arguments\n";
      opserr << "Want: uniaxialMaterial BarSlip tag? fc? fy? Es? fu? Eh? db? "
                "ld? nb? width? depth? bsflag? type? <damage? unit?>"
             << "\n";
      return TCL_ERROR;
    }

    int tag, nb, bsf, typ, dmg=1, unt=1;
    double fc, fy, Es, fu, Eh, ld, width, depth, db;

    int argStart = 2;

    if (Tcl_GetInt(interp, argv[argStart++], &tag) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid tag\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[argStart++], &fc) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid fc\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[argStart++], &fy) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid fy\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[argStart++], &Es) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid Es\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[argStart++], &fu) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid fu\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[argStart++], &Eh) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid Eh\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[argStart++], &db) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid db\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[argStart++], &ld) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid ld\n";
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[argStart++], &nb) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid nbars\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[argStart++], &width) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid width\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[argStart++], &depth) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid depth\n";
      return TCL_ERROR;
    }

    int y;
    y = argStart;

    if ((strcmp(argv[y], "strong") == 0) ||
        (strcmp(argv[y], "Strong") == 0) || (strcmp(argv[y], "weak") == 0) ||
        (strcmp(argv[y], "Weak") == 0)) {
      if ((strcmp(argv[y], "strong") == 0) ||
          (strcmp(argv[y], "Strong") == 0)) {
        bsf = 0;
      }

      if ((strcmp(argv[y], "weak") == 0) || (strcmp(argv[y], "Weak") == 0)) {
        bsf = 1;
      }
    } else {
      opserr << OpenSees::PromptValueError << "invalid bond strength specified\n";
      return TCL_ERROR;
    }
    y++;

    if ((strcmp(argv[y], "beamtop") == 0) ||
        (strcmp(argv[y], "beamTop") == 0) ||
        (strcmp(argv[y], "beambot") == 0) ||
        (strcmp(argv[y], "beamBot") == 0) ||
        (strcmp(argv[y], "beambottom") == 0) ||
        (strcmp(argv[y], "beamBottom") == 0) ||
        (strcmp(argv[y], "beam") == 0) || (strcmp(argv[y], "Beam") == 0) ||
        (strcmp(argv[y], "Column") == 0) ||
        (strcmp(argv[y], "column") == 0)) {
      if ((strcmp(argv[y], "beamtop") == 0) ||
          (strcmp(argv[y], "beamTop") == 0) ||
          (strcmp(argv[y], "beam") == 0) || (strcmp(argv[y], "Beam") == 0)) {
        typ = 0;
      }

      if ((strcmp(argv[y], "beambot") == 0) ||
          (strcmp(argv[y], "beamBot") == 0) ||
          (strcmp(argv[y], "beambottom") == 0) ||
          (strcmp(argv[y], "beamBottom") == 0)) {
        typ = 1;
      }

      if ((strcmp(argv[y], "column") == 0) ||
          (strcmp(argv[y], "Column") == 0)) {
        typ = 2;
      }
    } else {
      opserr << OpenSees::PromptValueError << "invalid location of bar specified\n";
      return TCL_ERROR;
    }
    if (argc == 17) {
      y++;

      if ((strcmp(argv[y], "damage1") == 0) ||
          (strcmp(argv[y], "Damage1") == 0) ||
          (strcmp(argv[y], "damage2") == 0) ||
          (strcmp(argv[y], "Damage2") == 0) ||
          (strcmp(argv[y], "nodamage") == 0) ||
          (strcmp(argv[y], "Nodamage") == 0) ||
          (strcmp(argv[y], "NoDamage") == 0) ||
          (strcmp(argv[y], "noDamage") == 0)) {
        if ((strcmp(argv[y], "damage1") == 0) ||
            (strcmp(argv[y], "Damage1") == 0)) {
          dmg = 1;
        } else if ((strcmp(argv[y], "damage2") == 0) ||
                   (strcmp(argv[y], "Damage2") == 0)) {
          dmg = 2;
        } else if ((strcmp(argv[y], "nodamage") == 0) ||
                   (strcmp(argv[y], "Nodamage") == 0) ||
                   (strcmp(argv[y], "NoDamage") == 0) ||
                   (strcmp(argv[y], "noDamage") == 0)) {
          dmg = 0;
        }

      } else {
        opserr << OpenSees::PromptValueError << "invalid damage specified\n";
        return TCL_ERROR;
      }

      y++;

      if ((strcmp(argv[y], "mpa") == 0) || (strcmp(argv[y], "MPa") == 0) ||
          (strcmp(argv[y], "mPa") == 0) || (strcmp(argv[y], "Mpa") == 0) ||
          (strcmp(argv[y], "psi") == 0) || (strcmp(argv[y], "Psi") == 0) ||
          (strcmp(argv[y], "PSI") == 0) || (strcmp(argv[y], "Pa") == 0) ||
          (strcmp(argv[y], "pa") == 0) || (strcmp(argv[y], "psf") == 0) ||
          (strcmp(argv[y], "Psf") == 0) || (strcmp(argv[y], "PSF") == 0) ||
          (strcmp(argv[y], "ksi") == 0) || (strcmp(argv[y], "Ksi") == 0) ||
          (strcmp(argv[y], "KSI") == 0) || (strcmp(argv[y], "ksf") == 0) ||
          (strcmp(argv[y], "Ksf") == 0) || (strcmp(argv[y], "KSF") == 0)) {
        if ((strcmp(argv[y], "mpa") == 0) || (strcmp(argv[y], "MPa") == 0) ||
            (strcmp(argv[y], "mPa") == 0) || (strcmp(argv[y], "Mpa") == 0)) {
          unt = 1;
        } else if ((strcmp(argv[y], "psi") == 0) ||
                   (strcmp(argv[y], "Psi") == 0) ||
                   (strcmp(argv[y], "PSI") == 0)) {
          unt = 2;
        } else if ((strcmp(argv[y], "Pa") == 0) ||
                   (strcmp(argv[y], "pa") == 0)) {
          unt = 3;
        } else if ((strcmp(argv[y], "psf") == 0) ||
                   (strcmp(argv[y], "Psf") == 0) ||
                   (strcmp(argv[y], "PSF") == 0)) {
          unt = 4;
        } else if ((strcmp(argv[y], "ksi") == 0) ||
                   (strcmp(argv[y], "Ksi") == 0) ||
                   (strcmp(argv[y], "KSI") == 0)) {
          unt = 5;
        } else if ((strcmp(argv[y], "ksf") == 0) ||
                   (strcmp(argv[y], "Ksf") == 0) ||
                   (strcmp(argv[y], "KSF") == 0)) {
          unt = 6;
        }
      } else {
        opserr << OpenSees::PromptValueError << "invalid unit specified\n";
        return TCL_ERROR;
      }
    }

    // allocate the material
    if (argc == 15) {
      theMaterial = new BarSlipMaterial(tag, fc, fy, Es, fu, Eh, db, ld, nb,
                                        width, depth, bsf, typ);
    }

    if (argc == 17) {
      theMaterial = new BarSlipMaterial(tag, fc, fy, Es, fu, Eh, db, ld, nb,
                                        width, depth, bsf, typ, dmg, unt);
    }

  }

  else if (strcmp(argv[1], "ShearPanel") == 0) {
    if (argc != 42 && argc != 31) {
      opserr << OpenSees::PromptValueError << "insufficient arguments\n";
      opserr << "Want: uniaxialMaterial ShearPanel tag? stress1p? strain1p? "
                "stress2p? strain2p? stress3p? strain3p? stress4p? strain4p? "
             << "\n<stress1n? strain1n? stress2n? strain2n? stress3n? "
                "strain3n? stress4n? strain4n?> rDispP? rForceP? uForceP? "
             << "\n<rDispN? rForceN? uForceN?> gammaK1? gammaK2? gammaK3? "
                "gammaK4? gammaKLimit? gammaD1? gammaD2? gammaD3? gammaD4? "
             << "\ngammaDLimit? gammaF1? gammaF2? gammaF3? gammaF4? "
                "gammaFLimit? gammaE? YieldStress? ";
      return TCL_ERROR;
    }

    int tag;
    double stress1p, stress2p, stress3p, stress4p;
    double strain1p, strain2p, strain3p, strain4p;
    double stress1n, stress2n, stress3n, stress4n;
    double strain1n, strain2n, strain3n, strain4n;
    double rDispP, rForceP, uForceP, rDispN, rForceN, uForceN;
    double gammaK1, gammaK2, gammaK3, gammaK4, gammaKLimit;
    double gammaD1, gammaD2, gammaD3, gammaD4, gammaDLimit;
    double gammaF1, gammaF2, gammaF3, gammaF4, gammaFLimit;
    double gammaE, yStr;

    int i = 2;

    if (Tcl_GetInt(interp, argv[i++], &tag) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid uniaxialMaterial ShearPanel tag" << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[i++], &stress1p) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid stress1p\n";
      opserr << "ShearPanel material: " << tag << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[i++], &strain1p) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid strain1p\n";
      opserr << "ShearPanel material: " << tag << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[i++], &stress2p) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid stress2p\n";
      opserr << "ShearPanel material: " << tag << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[i++], &strain2p) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid strain2p\n";
      opserr << "ShearPanel material: " << tag << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[i++], &stress3p) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid stress3p\n";
      opserr << "ShearPanel material: " << tag << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[i++], &strain3p) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid strain3p\n";
      opserr << "ShearPanel material: " << tag << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[i++], &stress4p) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid stress4p\n";
      opserr << "ShearPanel material: " << tag << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[i++], &strain4p) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid strain4p\n";
      opserr << "ShearPanel material: " << tag << "\n";
      return TCL_ERROR;
    }

    if (argc == 42) {
      if (Tcl_GetDouble(interp, argv[i++], &stress1n) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "invalid stress1n\n";
        opserr << "ShearPanel material: " << tag << "\n";
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[i++], &strain1n) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "invalid strain1n\n";
        opserr << "ShearPanel material: " << tag << "\n";
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[i++], &stress2n) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "invalid stress2n\n";
        opserr << "ShearPanel material: " << tag << "\n";
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[i++], &strain2n) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "invalid strain2n\n";
        opserr << "ShearPanel material: " << tag << "\n";
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[i++], &stress3n) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "invalid stress3n\n";
        opserr << "ShearPanel material: " << tag << "\n";
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[i++], &strain3n) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "invalid strain3n\n";
        opserr << "ShearPanel material: " << tag << "\n";
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[i++], &stress4n) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "invalid stress4n\n";
        opserr << "ShearPanel material: " << tag << "\n";
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[i++], &strain4n) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "invalid strain4n\n";
        opserr << "ShearPanel material: " << tag << "\n";
        return TCL_ERROR;
      }
    }

    if (Tcl_GetDouble(interp, argv[i++], &rDispP) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid rDispP\n";
      opserr << "ShearPanel material: " << tag << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[i++], &rForceP) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid rForceP\n";
      opserr << "ShearPanel material: " << tag << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[i++], &uForceP) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid uForceP\n";
      opserr << "ShearPanel material: " << tag << "\n";
      return TCL_ERROR;
    }

    if (argc == 42) {
      if (Tcl_GetDouble(interp, argv[i++], &rDispN) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "invalid rDispN\n";
        opserr << "ShearPanel material: " << tag << "\n";
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[i++], &rForceN) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "invalid rForceN\n";
        opserr << "ShearPanel material: " << tag << "\n";
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[i++], &uForceN) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "invalid uForceN\n";
        opserr << "ShearPanel material: " << tag << "\n";
        return TCL_ERROR;
      }
    }

    if (Tcl_GetDouble(interp, argv[i++], &gammaK1) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid gammaK1\n";
      opserr << "ShearPanel material: " << tag << "\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[i++], &gammaK2) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid gammaK2\n";
      opserr << "ShearPanel material: " << tag << "\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[i++], &gammaK3) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid gammaK3\n";
      opserr << "ShearPanel material: " << tag << "\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[i++], &gammaK4) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid gammaK4\n";
      opserr << "ShearPanel material: " << tag << "\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[i++], &gammaKLimit) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid gammaKLimit\n";
      opserr << "ShearPanel material: " << tag << "\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[i++], &gammaD1) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid gammaD1\n";
      opserr << "ShearPanel material: " << tag << "\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[i++], &gammaD2) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid gammaD2\n";
      opserr << "ShearPanel material: " << tag << "\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[i++], &gammaD3) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid gammaD3\n";
      opserr << "ShearPanel material: " << tag << "\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[i++], &gammaD4) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid gammaD4\n";
      opserr << "ShearPanel material: " << tag << "\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[i++], &gammaDLimit) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid gammaDLimit\n";
      opserr << "ShearPanel material: " << tag << "\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[i++], &gammaF1) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid gammaF1\n";
      opserr << "ShearPanel material: " << tag << "\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[i++], &gammaF2) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid gammaF2\n";
      opserr << "ShearPanel material: " << tag << "\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[i++], &gammaF3) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid gammaF3\n";
      opserr << "ShearPanel material: " << tag << "\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[i++], &gammaF4) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid gammaF4\n";
      opserr << "ShearPanel material: " << tag << "\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[i++], &gammaFLimit) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid gammaFLimit\n";
      opserr << "ShearPanel material: " << tag << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[i++], &gammaE) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid gammaE\n";
      opserr << "ShearPanel material: " << tag << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[i++], &yStr) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid yield stress\n";
      opserr << "ShearPanel material: " << tag << "\n";
      return TCL_ERROR;
    }

    // allocate the pinching material
    if (argc == 42) {
      theMaterial = new ShearPanelMaterial(
          tag, stress1p, strain1p, stress2p, strain2p, stress3p, strain3p,
          stress4p, strain4p, stress1n, strain1n, stress2n, strain2n,
          stress3n, strain3n, stress4n, strain4n, rDispP, rForceP, uForceP,
          rDispN, rForceN, uForceN, gammaK1, gammaK2, gammaK3, gammaK4,
          gammaKLimit, gammaD1, gammaD2, gammaD3, gammaD4, gammaDLimit,
          gammaF1, gammaF2, gammaF3, gammaF4, gammaFLimit, gammaE, yStr);
    }
    if (argc == 31) {
      theMaterial = new ShearPanelMaterial(
          tag, stress1p, strain1p, stress2p, strain2p, stress3p, strain3p,
          stress4p, strain4p, rDispP, rForceP, uForceP, gammaK1, gammaK2,
          gammaK3, gammaK4, gammaKLimit, gammaD1, gammaD2, gammaD3, gammaD4,
          gammaDLimit, gammaF1, gammaF2, gammaF3, gammaF4, gammaFLimit,
          gammaE, yStr);
    }
  }

  else if (strcmp(argv[1], "Concrete01WithSITC") == 0) {
    if (argc < 7) {
      opserr << OpenSees::PromptValueError << "insufficient arguments\n";
      opserr << "Want: uniaxialMaterial Concrete01 tag? fpc? epsc0? fpcu? "
                "epscu? <endStrainSITC?>"
             << "\n";
      return TCL_ERROR;
    }

    int tag;

    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid uniaxialMaterial Concrete01 tag" << "\n";
      return TCL_ERROR;
    }

    // Read required Concrete01 material parameters
    double fpc, epsc0, fpcu, epscu;

    if (Tcl_GetDouble(interp, argv[3], &fpc) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid fpc\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[4], &epsc0) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid epsc0\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[5], &fpcu) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid fpcu\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[6], &epscu) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid epscu\n";
      return TCL_ERROR;
    }

    if (argc == 7)
      // Parsing was successful, allocate the material
      theMaterial = new Concrete01WithSITC(tag, fpc, epsc0, fpcu, epscu);
    else {
      double endStrainSITC;
      if (Tcl_GetDouble(interp, argv[7], &endStrainSITC) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "invalid epscu\n";
        return TCL_ERROR;
      }
      theMaterial =
          new Concrete01WithSITC(tag, fpc, epsc0, fpcu, epscu, endStrainSITC);
    }
  }

  if (theMaterial == nullptr)
    return TCL_ERROR;

  if (builder->addTaggedObject<UniaxialMaterial>(*theMaterial) != TCL_OK) {
    delete theMaterial;
    return TCL_ERROR;
  }

  return TCL_OK;
}

