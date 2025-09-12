//===----------------------------------------------------------------------===//
//
//                                   xara
//                              https://xara.so
//
//===----------------------------------------------------------------------===//
//
// Copyright (c) 2025, OpenSees/Xara Developers
// All rights reserved.  No warranty, explicit or implicit, is provided.
//
// This source code is licensed under the BSD 2-Clause License.
// See LICENSE file or https://opensource.org/licenses/BSD-2-Clause
//
//===----------------------------------------------------------------------===//
//
#include <string>
#include <Logging.h>
#include <Parsing.h>

#include <Domain.h>
#include "DomainModalProperties.h"
#include "ResponseSpectrumAnalysis.h"


namespace OpenSees {

Tcl_CmdProc responseSpectrumAnalysis;

namespace DomainCommands {

int
modalProperties(ClientData clientData, 
                Tcl_Interp *interp, Tcl_Size argc, TCL_Char ** const argv)
{
  // modalProperties <-print> <-file $fileName> <-unorm>

  // some kudos
  static bool first_done = false;
  if (!first_done) {
      opserr << "Using DomainModalProperties - Developed by: Massimo Petracca, Guido Camata, ASDEA Software Technology\n";
      first_done = true;
  }

  // init default values
  bool unorm = false; // by default do not displacement-normalize eigenvectors 
  bool print_on_console = false; // by default do not print on console
  bool print_on_file = false;    // by default do no print on file
  std::string fname;

  // check options
  int loc = 1;
  while (loc < argc) {
    if (strcmp(argv[loc], "-unorm") == 0) {
      unorm = true;
    }
    else if (strcmp(argv[loc], "-print") == 0) {
      print_on_console = true;
    }
    else if (strcmp(argv[loc], "-file") == 0) {
      print_on_file = true;
      if (loc < argc - 1) {
        ++loc;
        fname = argv[loc];
      }
      else {
          opserr << "Error in modalProperties; "
                "After the keyword -file you should specify the file name.\n";
          return TCL_ERROR;
      }
    }
    ++loc;
  }


  Domain* domain = static_cast<Domain*>(clientData);

  // create the modal properties, compute them,
  // and add them to the domain
  DomainModalProperties *modal_props = new DomainModalProperties(domain, unorm);
  if (!modal_props->compute())
    return TCL_ERROR;
  // domain->setModalProperties(modal_props);

  // report
  if (print_on_console)
    modal_props->print();

  if (print_on_file)
    modal_props->print(fname);


  Tcl_CreateCommand(interp, "responseSpectrumAnalysis",
                    responseSpectrumAnalysis, modal_props, nullptr);

  return TCL_OK;
}

} // namespace DomainCommands


int
responseSpectrumAnalysis(ClientData clientData, 
                  Tcl_Interp* interp,
                  Tcl_Size argc, 
                  const char** const argv)
{
  // responseSpectrum $tsTag $dir <-scale $scale>
  //

  // some kudos
  static bool first_done = false;
  if (!first_done) {
    opslog << "Using ResponseSpectrumAnalysis - Developed by: Massimo Petracca, Guido Camata, ASDEA Software Technology\n";
    first_done = true;
  }

  if (clientData == nullptr) {
    opserr << "ResponseSpectrumAnalysis - eigen and modalProperties have not been called" << "\n";
    return TCL_ERROR;
  }


  // make sure eigenvalue and modal properties have been called before
  DomainModalProperties& modal_props = *static_cast<DomainModalProperties*>(clientData);

  int ndf = modal_props.totalMass().Size();


  //
  // parse
  //
  // ResponseSpectrumAnalysis $tsTag $dir <-scale $scale> <-damp $damp>
  // or
  // ResponseSpectrumAnalysis $dir -Tn $TnValues -fn $fnValues -Sa $SaValues <-scale $scale> <-damp $damp>
  //

  // init default arguments
  TimeSeries* ts = nullptr;
  int dir = 1;
  double scale = 1.0;
  std::vector<double> Tn;
  std::vector<double> Sa;
  int mode_id = 0;
  bool single_mode = false;


  int nargs = argc - 1; // skip the command name
  if (nargs < 2) {
    opserr << "ResponseSpectrumAnalysis $tsTag $dir <-scale $scale> <-damp $damp>\n"
      << "or\n"
      << "ResponseSpectrumAnalysis $dir -Tn $TnValues -fn $fnValues -Sa $SaValues <-scale $scale> <-damp $damp>\n"
      "Error: at least 2 arguments should be provided.\n";
    return TCL_ERROR;
  }

  // search for -Tn -Sa, if both found we use them as lists
  // otherwise we fallback to the old implementation of the timeSeries
  bool found_Tn = false;
  bool found_Sa = false;
  for (int i=1; i<argc; i++) {
    if (!found_Tn && (strcmp(argv[i], "-Tn") == 0 || strcmp(argv[i], "-fn") == 0))
      found_Tn = true;
    if (!found_Sa && strcmp(argv[i], "-Sa") == 0)
      found_Sa = true;
  }

  bool use_lists = found_Tn && found_Sa;
  if (found_Tn && !found_Sa) {
    opserr << "ResponseSpectrumAnalysis Error: found -Tn or -fn but not -Sa, please specify both of them or use a timeSeries\n";
    return TCL_ERROR;
  }
  if (found_Sa && !found_Tn) {
    opserr << "ResponseSpectrumAnalysis Error: found -Sa but not -Tn or -fn, please specify both of them or use a timeSeries\n";
    return TCL_ERROR;
  }

  int argi = 1;

  // get time series;
  if (!use_lists) {
    // command is 
    //   ResponseSpectrumAnalysis $tsTag $dir <-scale $scale> <-damp $damp>

    int tstag;
    if (Tcl_GetInt(interp, argv[argi], &tstag) != TCL_OK) {
      opserr << "ResponseSpectrumAnalysis Error: Failed to get timeSeries tag.\n";
      return TCL_ERROR;
    }

    // ts = OPS_getTimeSeries(tstag);
    if (ts == nullptr) {
      opserr << "ResponseSpectrumAnalysis Error: Failed to get timeSeries with tag = " << tstag << ".\n";
      return TCL_ERROR;
    }

    argi++; // move to the next argument
  }

  // get direction
  if (Tcl_GetInt(interp, argv[argi], &dir) != TCL_OK) {
    opserr << "ResponseSpectrumAnalysis Error: Failed to get direction.\n";
    return TCL_ERROR;
  }
  if (dir < 1 || dir > ndf) {
    opserr << "ResponseSpectrumAnalysis Error: provided direction (" << dir << ") should be in the range 1-" << ndf << ".\n";
    return TCL_ERROR;
  }

  // parse optional data
  for (int i=argi; i<argc; i++) {
    const char* value = argv[i];

    if (strcmp(value, "-scale") == 0) {
      if (i+1 >= argc) {
        opserr << "ResponseSpectrumAnalysis Error: scale factor requested but not provided.\n";
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[i+1], &scale) < 0) {
        opserr << "ResponseSpectrumAnalysis Error: Failed to get scale factor.\n";
        return TCL_ERROR;
      }
      i++;
    }

    else if (strcmp(value, "-mode") == 0) {
      if (i + 1 >= argc) {
        opserr << "ResponseSpectrumAnalysis Error: mode_id requested but not provided.\n";
        return TCL_ERROR;
      }

      if (Tcl_GetInt(interp, argv[i + 1], &mode_id) < 0) {
        opserr << "ResponseSpectrumAnalysis Error: Failed to get mode_id.\n";
        return TCL_ERROR;
      }
      --mode_id; // make it 0-based
      single_mode = true;
      i++;
    }

    else if (strcmp(value, "-Tn") == 0 || strcmp(value, "-fn") == 0) {
      // first try expanded list like {*}$the_list,
      // also used in python like *the_list
      Tn.clear();
      if (i + 1 >= argc) {
        opserr << "ResponseSpectrumAnalysis Error: Tn values requested but not provided.\n";
        return TCL_ERROR;
      }

      for (int j=i+1; j<argc; j++) {
        double item;
        if (Tcl_GetDouble(interp, argv[j], &item) != TCL_OK) {
          break;
        }

        if (strcmp(value, "-fn") == 0 && item != 0) {
          item = 1.0 / item;
        }
        Tn.push_back(item);
      }

      // try proper list
      if (Tn.size() == 0) {
        int tn_argc;
        const char** tn_argv = nullptr;
        if (Tcl_SplitList(interp, argv[i + 1], &tn_argc, &tn_argv) != TCL_OK) {
          opserr << "ResponseSpectrumAnalysis Error: cannot parse the Tn list.\n";
          return TCL_ERROR;
        }
        for (int j = 0; j < tn_argc; ++j) {
          double item;
          if (Tcl_GetDouble(interp, tn_argv[j], &item) != TCL_OK) {
            opserr << "ResponseSpectrumAnalysis Error: cannot parse the Tn list.\n";
            Tcl_Free((char*)tn_argv);
            return TCL_ERROR;
          }
          if (strcmp(value, "-fn") == 0 && item != 0) {
            item = 1.0 / item;
          }
          Tn.push_back(item);
        }
        Tcl_Free((char*)tn_argv);
      }
    }

    else if (strcmp(value, "-Sa") == 0) {
      // first try expanded list like {*}$the_list,
      // also used in python like *the_list
      Sa.clear();
      if (i + 1 >= argc) {
        opserr << "ResponseSpectrumAnalysis Error: Sa values requested but not provided.\n";
        return TCL_ERROR;
      }

      for (int j = i + 1; j < argc; j++) {
        double item;
        if (Tcl_GetDouble(interp, argv[j], &item) != TCL_OK) {
          break;
        }
        Sa.push_back(item);
      }

      if (Sa.size() == 0) {
        // try Tcl list
        int sa_argc;
        const char** sa_argv = nullptr;
        if (Tcl_SplitList(interp, argv[i + 1], &sa_argc, &sa_argv) != TCL_OK) {
          opserr << "ResponseSpectrumAnalysis Error: cannot parse the Sa list.\n";
          return TCL_ERROR;
        }
        for (int j = 0; j < sa_argc; ++j) {
          double item;
          if (Tcl_GetDouble(interp, sa_argv[j], &item) != TCL_OK) {
            opserr << "ResponseSpectrumAnalysis Error: cannot parse the Sa list.\n";
            Tcl_Free((char*)sa_argv);
            return TCL_ERROR;
          }
          Sa.push_back(item);
        }
        Tcl_Free((char*)sa_argv);
      }
    }
  }

  // check Tn and Sa vectors
  if (use_lists) {
    if (Tn.size() != Sa.size()) {
      opserr << "ResponseSpectrumAnalysis Error: Sa and Tn lists must have the same length\n";
      opserr << (int)Tn.size() << " != " << (int)Sa.size() << "\n";
      return TCL_ERROR;
    }
    if (Tn.size() == 0) {
      opserr << "ResponseSpectrumAnalysis Error: Sa and Tn lists cannot be empty\n";
      return TCL_ERROR;
    }
    for (std::size_t i = 0; i < Tn.size(); ++i) {
      if (Tn[i] < 0.0) {
        opserr << "ResponseSpectrumAnalysis Error: Tn values must be positive (found " << Tn[i] << ")\n";
        return TCL_ERROR;
      }
      if (i > 0 && Tn[i] <= Tn[i - 1]) {
        opserr << "ResponseSpectrumAnalysis Error: Tn values must be monotonically increasing (found " << Tn[i] << " after " << Tn[i - 1] << ")\n";
        return TCL_ERROR;
      }
    }
    for (double item : Sa) {
      if (item < 0.0) {
        opserr << "ResponseSpectrumAnalysis Error: Sa values must be positive (found " << item << ")\n";
        return TCL_ERROR;
      }
    }
  }

  // Create the response spectrum analysis and run it
  ResponseSpectrumAnalysis rsa(modal_props, ts, Tn, Sa, dir, scale);
  int result;
  if (single_mode)
    result = rsa.analyze(mode_id);
  else
    result = rsa.analyze();

  return result == 0? TCL_OK : TCL_ERROR;
}


} // namespace OpenSees