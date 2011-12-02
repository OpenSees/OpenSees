///////////////////////////////////////////////////////////////////////////////
//   COPYLEFT (C): Woody's viral GPL-like license (by BJ):
//                 ``This    source  code is Copyrighted in
//                 U.S.,  for  an  indefinite  period,  and anybody
//                 caught  using it without our permission, will be
//                 mighty good friends of ourn, cause we don't give
//                 a  darn.  Hack it. Compile it. Debug it. Run it.
//                 Yodel  it.  Enjoy it. We wrote it, that's all we
//                 wanted to do.''
//
//
// COPYRIGHT (C):     :-))
// PROJECT:           Object Oriented Finite Element Program
// FILE:              
// CLASS:             
// MEMBER FUNCTIONS:
//
// MEMBER VARIABLES
//
// PURPOSE:           
//
// RETURN:
// VERSION:
// LANGUAGE:          C++
// TARGET OS:         
// DESIGNER:          Zhao Cheng, Boris Jeremic
// PROGRAMMER:        Zhao Cheng, 
// DATE:              Fall 2005
// UPDATE HISTORY:    
//
///////////////////////////////////////////////////////////////////////////////
//

#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream.h>

#include <Domain.h>

#include <tcl.h>
#include <ErrorHandler.h>
#include <TclModelBuilder.h>
#include <OPS_Globals.h>
#include <ConsoleErrorHandler.h>

#include <stresst.h>

#include "NewTemplate3Dep.h"

#include "MaterialParameter.h"

#include "ElasticState.h"
#include "Isotropic_Elastic.h"
#include "elnp_Elastic.h"
#include "PressureDependent_Elastic.h"
#include "DM04_Elastic.h"

#include "YieldFunction.h"
#include "DP_YF.h"
#include "VM_YF.h"
#include "CC_YF.h"
#include "DM04_YF.h"

#include "PlasticFlow.h"
#include "DP_PF.h"
#include "VM_PF.h"
#include "CC_PF.h"
#include "DM04_PF.h"

#include "ScalarEvolution.h"
#include "Linear_Eeq.h"
#include "CC_Ev.h"


#include "TensorEvolution.h"
#include "Linear_Eij.h"
#include "AF_Eij.h"
#include "DM04_alpha_Eij.h"
#include "DM04_z_Eij.h"

int n_mc = 0;
int n_is = 0;
int n_it = 0;

MaterialParameter*  EvaluateMaterialParameter(ClientData, Tcl_Interp*, TCL_Char* tclString);
ElasticState*       EvaluateElasticState(ClientData, Tcl_Interp*, TCL_Char* tclString);
YieldFunction*      EvaluateYieldFunction(ClientData, Tcl_Interp*, TCL_Char* tclString);
PlasticFlow*        EvaluatePlasticFlow(ClientData, Tcl_Interp*, TCL_Char* tclString);
ScalarEvolution**   EvaluateSE(ClientData, Tcl_Interp*, TCL_Char* tclString);
TensorEvolution**   EvaluateTE(ClientData, Tcl_Interp*, TCL_Char* tclString);


static void cleanup(TCL_Char **argv) {
    Tcl_Free((char *) argv);
}


//************************************************************************************************
//************************************************************************************************
NewTemplate3Dep*
TclModelBuilder_addNewTemplate3Dep(ClientData clientData, Tcl_Interp *interp,  int argc,
          TCL_Char **argv, TclModelBuilder *theTclBuilder, int eleArgStart)
{
  int tag = 0;
  MaterialParameter*  MatP= NULL;
  ElasticState*       ES = NULL;
  YieldFunction*      YF = NULL;
  PlasticFlow*        PF = NULL;
  ScalarEvolution**   SSE = NULL;
  TensorEvolution**   TTE = NULL;
  int CI = 0;

  int loc = eleArgStart;

  if (Tcl_GetInt(interp, argv[loc++], &tag) != TCL_OK) {
    cout << "Warning:  NewTemplate3Dep - invalid tag " << argv[loc] << endl;
    exit (1);
  }

  while (loc < argc) {
    
    if ( strcmp(argv[loc],"-MaterialParameter") == 0 ) {
      MatP = EvaluateMaterialParameter(clientData, interp, argv[loc+1]);
      if (MatP == NULL) {
        cout << "Warning:  TclNewTemplate3Dep - could not create elastic state from " << argv[loc+1] << endl;
        exit (1);
      }
      loc += 2;
    }
    
    else if ( strcmp(argv[loc],"-ElasticState") == 0 ) {
      ES = EvaluateElasticState(clientData, interp, argv[loc+1]);
      if (ES == NULL) {
        cout << "Warning:  TclNewTemplate3Dep - could not create elastic state from " << argv[loc+1] << endl;
        exit (1);
      }
      loc += 2;
    }

    else if ( strcmp(argv[loc],"-YieldFunction") == 0 ) {
      YF = EvaluateYieldFunction(clientData, interp, argv[loc+1]);
      if (YF == NULL) {
        cout << "Warning:  TclNewTemplate3Dep - could not create a yield function from " << argv[loc+1] << endl;
        exit (1);
      }
      loc += 2;
    }

    else if ( strcmp(argv[loc],"-PlasticFlow") == 0 ) {
      PF = EvaluatePlasticFlow(clientData, interp, argv[loc+1]);
      if (PF == NULL) {
        cout << "Warning:  TclNewTemplate3Dep - could not create a yield function from " << argv[loc+1] << endl;
        exit (1);
      }
      loc += 2;
    }

    else if ( strcmp(argv[loc],"-ScalarEvolution") == 0 ) {
      SSE = EvaluateSE(clientData, interp, argv[loc+1]);
      if (SSE == NULL) {
        cout << "Warning:  TclNewTemplate3Dep - could not create a scalar evolution from " << argv[loc+1] << endl;
        exit (1);
      }
      loc += 2;
    }

    else if ( strcmp(argv[loc],"-TensorEvolution") == 0 ) {
      TTE = EvaluateTE(clientData, interp, argv[loc+1]);
      if (TTE == NULL) {
        cout << "Warning:  TclNewTemplate3Dep - could not create a tensor evolution from " << argv[loc+1] << endl;
        exit (1);
      }
      loc += 2;
    }

    else if ( strcmp(argv[loc],"-Algorithm") == 0 ) {
      if (Tcl_GetInt(interp, argv[loc+1], &CI) != TCL_OK) {
        cout << "Warning:  TclNewTemplate3Dep - invalid  algorithm index " << argv[loc] << endl;
        exit (1);
      }
      loc += 2;
    }
    
    else {
        cout << "Warning:  TclNewTemplate3Dep - unknown keyword/command : " << argv[loc] << endl;
        exit (1);
      }

  }

  NewTemplate3Dep *theMaterial = 0;

  if ( (MatP != NULL) && (ES != NULL) && (YF != NULL) && (PF != NULL) && (SSE != NULL) )
    theMaterial = new NewTemplate3Dep(tag, MatP, ES, YF, PF, SSE, TTE, CI);
  else if ( (MatP != NULL) && (ES != NULL) && (YF != NULL) && (PF != NULL) && (SSE == NULL) )
    theMaterial = new NewTemplate3Dep(tag, MatP, ES, YF, PF, TTE, CI);
  else
    cout << "Warning: invalid args used to create a NewTemplate3Dep material." << endl;

  return theMaterial;
}


//**************************************************************************************
// Function - to create an MaterialParameter  object
MaterialParameter *EvaluateMaterialParameter(ClientData clientData, Tcl_Interp *interp, TCL_Char *tclString)
{
  int argc;
  TCL_Char **argv;

  // split the list
  if (Tcl_SplitList(interp, tclString, &argc, &argv) != TCL_OK)
    exit (1);

  if (argc == 0)
    exit (1);

  int loc = 0;
  double* mc = NULL;
  double* is = NULL;
  stresstensor* it = NULL;
  MaterialParameter *inp = NULL;

  while (loc < argc) {
    
    if ( strcmp(argv[loc],"MaterialConstant") == 0  ) {      

      loc++;  
      
      if (Tcl_GetInt(interp, argv[loc], &n_mc) != TCL_OK) {
        cout << "Warning:  TclNewTemplate3Dep - invalid input " << argv[loc] << endl;
        exit (1);
      } 
 
      //cout << "n_mc = " << n_mc << endl;
      
      if (n_mc > 0)
        mc = new double[n_mc];

      if (mc == 0) {
        cout << "Warning:  TclNewTemplate3Dep - could not create material constants " << endl;
        exit (1);
      }

      for (int i = 0; i < n_mc; i++) {
        loc++;
        if (Tcl_GetDouble(interp, argv[loc], &mc[i]) != TCL_OK) {
          cout << "Warning:  TclNewTemplate3Dep - invalid input " << argv[loc] << endl;
          exit (1);
        }
      }            
    }
    
    else if ( strcmp(argv[loc],"InternalScalar") == 0 ) {
 
      loc++;

      if (Tcl_GetInt(interp, argv[loc], &n_is) != TCL_OK) {
        cout << "Warning:  TclNewTemplate3Dep - invalid input " << argv[loc] << endl;
        exit (1);
      } 
      
      //cout << "n_is = " << n_is << endl;
      
      if (n_is > 0) {
        is = new double[n_is];
        if (is == 0) {
          cout << "Warning:  TclNewTemplate3Dep - could not create scalar variables " << endl;
          exit (1);
        }

        for (int i = 0; i < n_is; i++) {
          loc++;
          if (Tcl_GetDouble(interp, argv[loc], &is[i]) != TCL_OK) {
            cout << "Warning:  TclNewTemplate3Dep - invalid input " << argv[loc] << endl;
            exit (1);
          }
        }
      }
    }
    
    else if ( strcmp(argv[loc],"InternalTensor") == 0 ) {

      loc++;

      if (Tcl_GetInt(interp, argv[loc], &n_it) != TCL_OK) {
        cout << "Warning:  TclNewTemplate3Dep - invalid input " << argv[loc] << endl;
        exit (1);
      } 

      if (n_it > 0) {
        it = new stresstensor[n_it];
        double a[9];
        if (it == 0) {
          cout << "Warning:  TclNewTemplate3Dep - could not create tensor variables " << endl;
          exit (1);
        }
 
        //cout << "n_it = " << n_it << endl;
        
        for (int i = 0; i < n_it; i++) {
          for (int j = 0; j < 9; j++) {
            loc++;
            if (Tcl_GetDouble(interp, argv[loc], &a[j]) != TCL_OK) {
              cout << "Warning:  TclNewTemplate3Dep - invalid input " << argv[loc] << endl;
              exit (1);
            }
          }
          it[i].val(1,1) = a[0]; it[i].val(1,2) = a[1]; it[i].val(1,3) = a[2];
          it[i].val(2,1) = a[3]; it[i].val(2,2) = a[4]; it[i].val(2,3) = a[5];
          it[i].val(3,1) = a[6]; it[i].val(3,2) = a[7]; it[i].val(3,3) = a[8];
        }
      }         
    }

    loc++;       
  }

  // now parse the list & construct the required object
  if ( (mc != NULL) && (is != NULL) && (it != NULL) )
    inp = new MaterialParameter(mc, n_mc, is, n_is, it, n_it);
  else if ( (mc != NULL) && (is != NULL) && (it == NULL) )
    inp = new MaterialParameter(mc, n_mc, is, n_is);
  else if ( (mc != NULL) && (is == NULL) && (it != NULL) )
    inp = new MaterialParameter(mc, n_mc, it, n_it);
  else if ( (mc != NULL) && (is == NULL) && (it == NULL) )
    inp = new MaterialParameter(mc, n_mc);
  else {
    cout << "Warning:  TclNewTemplate3Dep - could not create input parameter " << endl;
    exit (1);
  }

  cleanup(argv);
  return inp;
}


//**************************************************************************************
// Function - to create a Elastic State object
ElasticState *EvaluateElasticState(ClientData clientData, Tcl_Interp *interp, TCL_Char *tclString)
{
  int argc;
  TCL_Char **argv;

  // split the list
  if (Tcl_SplitList(interp, tclString, &argc, &argv) != TCL_OK) {
    exit (1);
  }

  if (argc == 0)
    exit (1);

  // now parse the list & construct the required object
  ElasticState *ElS = NULL;

  // 1. Linear elastic model
  if ( strcmp(argv[0],"Isotropic") == 0 ) {
    int a1 = 0;
    int a2 = 0;
    stresstensor initStre;
    straintensor initStra;

    if (argc > 2 ) {
      if (Tcl_GetInt(interp, argv[1], &a1) != TCL_OK) {
        cout << "Warning:  TclNewTemplate3Dep - invalid input " << argv[1] << endl;
        exit (1);
      }
      if (Tcl_GetInt(interp, argv[2], &a2) != TCL_OK) {
        cout << "Warning:  TclNewTemplate3Dep - invalid input " << argv[2] << endl;
        exit (1);
      }
    }
    if (argc > (2 + 9) ) {
      double a[9];
        for (int j = 0; j < 9; j++) {
          if (Tcl_GetDouble(interp, argv[4+1+j], &a[j]) != TCL_OK) {
            cout << "Warning:  TclNewTemplate3Dep - invalid input " << argv[2+1+j] << endl;
            exit (1);
          }
        }
        initStre.val(1,1) = a[0]; initStre.val(1,2) = a[1]; initStre.val(1,3) = a[2];
        initStre.val(2,1) = a[3]; initStre.val(2,2) = a[4]; initStre.val(2,3) = a[5];
        initStre.val(3,1) = a[6]; initStre.val(3,2) = a[7]; initStre.val(3,3) = a[8];
    }
    if (argc > (2 + 18) ) {
      double a[9];
        for (int j = 0; j < 9; j++) {
          if (Tcl_GetDouble(interp, argv[2+10+j], &a[j]) != TCL_OK) {
            cout << "Warning:  TclNewTemplate3Dep - invalid input " << argv[2+10+j] << endl;
            exit (1);
          }
        }
        initStra.val(1,1) = a[0]; initStra.val(1,2) = a[1]; initStra.val(1,3) = a[2];
        initStra.val(2,1) = a[3]; initStra.val(2,2) = a[4]; initStra.val(2,3) = a[5];
        initStra.val(3,1) = a[6]; initStra.val(3,2) = a[7]; initStra.val(3,3) = a[8];
    }

    ElS = new Isotropic_Elastic(a1, a2, initStre, initStra);
  }

  // 2. elnp elastic model
  else if ( strcmp(argv[0],"elnp") == 0  ||  strcmp(argv[0],"Cam-Clay") == 0 ) {
    int a1 = 0;
    int a2 = 0;
    int a3 = 0;
    int a4 = 0;
    stresstensor initStre;
    straintensor initStra;

    if (argc > 4 ) {
      if (Tcl_GetInt(interp, argv[1], &a1) != TCL_OK) {
        cout << "Warning:  TclNewTemplate3Dep - invalid input " << argv[1] << endl;
        exit (1);
      }
      if (Tcl_GetInt(interp, argv[2], &a2) != TCL_OK) {
        cout << "Warning:  TclNewTemplate3Dep - invalid input " << argv[2] << endl;
        exit (1);
      }
      if (Tcl_GetInt(interp, argv[3], &a3) != TCL_OK) {
        cout << "Warning:  TclNewTemplate3Dep - invalid input " << argv[3] << endl;
        exit (1);
      }
      if (Tcl_GetInt(interp, argv[4], &a4) != TCL_OK) {
        cout << "Warning:  TclNewTemplate3Dep - invalid input " << argv[4] << endl;
        exit (1);
      }          
    }
    if (argc > (4 + 9) ) {
      double a[9];
        for (int j = 0; j < 9; j++) {
          if (Tcl_GetDouble(interp, argv[4+1+j], &a[j]) != TCL_OK) {
            cout << "Warning:  TclNewTemplate3Dep - invalid input " << argv[4+1+j] << endl;
            exit (1);
          }
        }
        initStre.val(1,1) = a[0]; initStre.val(1,2) = a[1]; initStre.val(1,3) = a[2];
        initStre.val(2,1) = a[3]; initStre.val(2,2) = a[4]; initStre.val(2,3) = a[5];
        initStre.val(3,1) = a[6]; initStre.val(3,2) = a[7]; initStre.val(3,3) = a[8];
    }
    if (argc > (4 + 18) ) {
      double a[9];
        for (int j = 0; j < 9; j++) {
          if (Tcl_GetDouble(interp, argv[4+10+j], &a[j]) != TCL_OK) {
            cout << "Warning:  TclNewTemplate3Dep - invalid input " << argv[4+10+j] << endl;
            exit (1);
          }
        }
        initStra.val(1,1) = a[0]; initStra.val(1,2) = a[1]; initStra.val(1,3) = a[2];
        initStra.val(2,1) = a[3]; initStra.val(2,2) = a[4]; initStra.val(2,3) = a[5];
        initStra.val(3,1) = a[6]; initStra.val(3,2) = a[7]; initStra.val(3,3) = a[8];
    }

    ElS = new elnp_Elastic(a1, a2, a3, a4, initStre, initStra);
  }

  // 3. PressureDependent elastic model
  else if (strcmp(argv[0],"PressureDependent") == 0) {
    int a1 = 0;
    int a2 = 0;
    int a3 = 0;
    int a4 = 0;
    int a5 = 0;
    stresstensor initStre;
    straintensor initStra;

    if (argc > 5 ) {
      if (Tcl_GetInt(interp, argv[1], &a1) != TCL_OK) {
        cout << "Warning:  TclNewTemplate3Dep - invalid input " << argv[1] << endl;
        exit (1);
      }
      if (Tcl_GetInt(interp, argv[2], &a2) != TCL_OK) {
        cout << "Warning:  TclNewTemplate3Dep - invalid input " << argv[2] << endl;
        exit (1);
      }
      if (Tcl_GetInt(interp, argv[3], &a3) != TCL_OK) {
        cout << "Warning:  TclNewTemplate3Dep - invalid input " << argv[3] << endl;
        exit (1);
      }    
      if (Tcl_GetInt(interp, argv[4], &a4) != TCL_OK) {
        cout << "Warning:  TclNewTemplate3Dep - invalid input " << argv[4] << endl;
        exit (1);
      }
      if (Tcl_GetInt(interp, argv[5], &a5) != TCL_OK) {
        cout << "Warning:  TclNewTemplate3Dep - invalid input " << argv[5] << endl;
        exit (1);
      } 
    }
    if (argc > (5 + 9) ) {
      double a[9];
        for (int j = 0; j < 9; j++) {
          if (Tcl_GetDouble(interp, argv[5+1+j], &a[j]) != TCL_OK) {
            cout << "Warning:  TclNewTemplate3Dep - invalid input " << argv[5+1+j] << endl;
            exit (1);
          }
        }
        initStre.val(1,1) = a[0]; initStre.val(1,2) = a[1]; initStre.val(1,3) = a[2];
        initStre.val(2,1) = a[3]; initStre.val(2,2) = a[4]; initStre.val(2,3) = a[5];
        initStre.val(3,1) = a[6]; initStre.val(3,2) = a[7]; initStre.val(3,3) = a[8];
    }
    if (argc > (5 + 18) ) {
      double a[9];
        for (int j = 0; j < 9; j++) {
          if (Tcl_GetDouble(interp, argv[4+10+j], &a[j]) != TCL_OK) {
            cout << "Warning:  TclNewTemplate3Dep - invalid input " << argv[5+10+j] << endl;
            exit (1);
          }
        }
        initStra.val(1,1) = a[0]; initStra.val(1,2) = a[1]; initStra.val(1,3) = a[2];
        initStra.val(2,1) = a[3]; initStra.val(2,2) = a[4]; initStra.val(2,3) = a[5];
        initStra.val(3,1) = a[6]; initStra.val(3,2) = a[7]; initStra.val(3,3) = a[8];
    }

    ElS = new PressureDependent_Elastic(a1, a2, a3, a4, a5, initStre, initStra);
  }

  // 4. DM04 elastic model
  else if ( strcmp(argv[0],"Dafalias-Manzari") == 0 ) {
    int a1 = 0;
    int a2 = 0;
    int a3 = 0;
    int a4 = 0;
    int a5 = 0;
    stresstensor initStre;
    straintensor initStra;

    if (argc > 5 ) {
      if (Tcl_GetInt(interp, argv[1], &a1) != TCL_OK) {
        cout << "Warning:  TclNewTemplate3Dep - invalid input " << argv[1] << endl;
        exit (1);
      }
      if (Tcl_GetInt(interp, argv[2], &a2) != TCL_OK) {
        cout << "Warning:  TclNewTemplate3Dep - invalid input " << argv[2] << endl;
        exit (1);
      }
      if (Tcl_GetInt(interp, argv[3], &a3) != TCL_OK) {
        cout << "Warning:  TclNewTemplate3Dep - invalid input " << argv[3] << endl;
        exit (1);
      }    
      if (Tcl_GetInt(interp, argv[4], &a4) != TCL_OK) {
        cout << "Warning:  TclNewTemplate3Dep - invalid input " << argv[4] << endl;
        exit (1);
      }
      if (Tcl_GetInt(interp, argv[5], &a5) != TCL_OK) {
        cout << "Warning:  TclNewTemplate3Dep - invalid input " << argv[5] << endl;
        exit (1);
      }      
    }
    if (argc > (5 + 9) ) {
      double a[9];
        for (int j = 0; j < 9; j++) {
          if (Tcl_GetDouble(interp, argv[5+1+j], &a[j]) != TCL_OK) {
            cout << "Warning:  TclNewTemplate3Dep - invalid input " << argv[5+1+j] << endl;
            exit (1);
          }
        }
        initStre.val(1,1) = a[0]; initStre.val(1,2) = a[1]; initStre.val(1,3) = a[2];
        initStre.val(2,1) = a[3]; initStre.val(2,2) = a[4]; initStre.val(2,3) = a[5];
        initStre.val(3,1) = a[6]; initStre.val(3,2) = a[7]; initStre.val(3,3) = a[8];
    }
    if (argc > (5 + 18) ) {
      double a[9];
        for (int j = 0; j < 9; j++) {
          if (Tcl_GetDouble(interp, argv[5+10+j], &a[j]) != TCL_OK) {
            cout << "Warning:  TclNewTemplate3Dep - invalid input " << argv[5+10+j] << endl;
            exit (1);
          }
        }
        initStra.val(1,1) = a[0]; initStra.val(1,2) = a[1]; initStra.val(1,3) = a[2];
        initStra.val(2,1) = a[3]; initStra.val(2,2) = a[4]; initStra.val(2,3) = a[5];
        initStra.val(3,1) = a[6]; initStra.val(3,2) = a[7]; initStra.val(3,3) = a[8];
    }

    ElS = new DM04_Elastic(a1, a2, a3, a4, a5, initStre, initStra);
  }

  // if others
  else {
    cout << "Warning: invalid elastic state object: " << argv[0] << endl;
    exit(1);
  }

  cleanup(argv);
  return ElS;
}


//**************************************************************************************
// Function - to create a Yield Function
YieldFunction *EvaluateYieldFunction(ClientData clientData, Tcl_Interp *interp, TCL_Char *tclString)
{
  int argc;
  TCL_Char **argv;

  // split the list
  if (Tcl_SplitList(interp, tclString, &argc, &argv) != TCL_OK) {
    exit (1);
  }

  if (argc == 0)
    exit (1);

  // now parse the list & construct the required object
  YieldFunction *YF = NULL;

  // 1. von Mises Yield Function
  if ((strcmp(argv[0],"Von-Mises") == 0) || (strcmp(argv[0],"von-Mises") == 0)) {
    int a1 = -1;  int b1 = 0;
    int a2 = -1;  int b2 = 0;

    if (argc <= 2)  {
        cout << "Warning:  TclNewTemplate3Dep - Yield Function (VM) - input parameters < 2" << endl;
        exit (1);
    }

    if (argc > 2 ) {
      if (Tcl_GetInt(interp, argv[1], &a1) != TCL_OK) {
        cout << "Warning:  TclNewTemplate3Dep - invalid input " << argv[1] << endl;
        exit (1);
      }
      if (Tcl_GetInt(interp, argv[2], &b1) != TCL_OK) {
        cout << "Warning:  TclNewTemplate3Dep - invalid input " << argv[2] << endl;
        exit (1);
      }
    }
    if (argc > 4 ) {
      if (Tcl_GetInt(interp, argv[3], &a2) != TCL_OK) {
        cout << "Warning:  TclNewTemplate3Dep - invalid input " << argv[3] << endl;
        exit (1);
      }
      if (Tcl_GetInt(interp, argv[4], &b2) != TCL_OK) {
        cout << "Warning:  TclNewTemplate3Dep - invalid input " << argv[4] << endl;
        exit (1);
      }
    }
  
    YF = new VM_YF(a1,b1, a2,b2);
  }

  // 2. Drucker-Prager Yield Function
  else if ( (strcmp(argv[0],"Drucker-Prager") == 0) ) {
    int a1 = -1;  int b1 = 0;
    int a2 = -1;  int b2 = 0;
    int a3 = -1;  int b3 = 0;

    if (argc <= 2)  {
        cout << "Warning:  TclNewTemplate3Dep - Yield Function (DP) - input parameters < 2" << endl;
        exit (1);
    }

    if (argc > 2 ) {
      if (Tcl_GetInt(interp, argv[1], &a1) != TCL_OK) {
        cout << "Warning:  TclNewTemplate3Dep - invalid input " << argv[1] << endl;
        exit (1);
      }
      if (Tcl_GetInt(interp, argv[2], &b1) != TCL_OK) {
        cout << "Warning:  TclNewTemplate3Dep - invalid input " << argv[2] << endl;
        exit (1);
      }
    }
    if (argc > 4 ) {
      if (Tcl_GetInt(interp, argv[3], &a2) != TCL_OK) {
        cout << "Warning:  TclNewTemplate3Dep - invalid input " << argv[3] << endl;
        exit (1);
      }
      if (Tcl_GetInt(interp, argv[4], &b2) != TCL_OK) {
        cout << "Warning:  TclNewTemplate3Dep - invalid input " << argv[4] << endl;
        exit (1);
      }
    }   
    if (argc > 6 ) {
      if (Tcl_GetInt(interp, argv[5], &a3) != TCL_OK) {
        cout << "Warning:  TclNewTemplate3Dep - invalid input " << argv[5] << endl;
        exit (1);
      }
      if (Tcl_GetInt(interp, argv[6], &b3) != TCL_OK) {
        cout << "Warning:  TclNewTemplate3Dep - invalid input " << argv[6] << endl;
        exit (1);
      }
    }

    YF = new DP_YF(a1,b1, a2,b2, a3, b3);
  }

  // 3. Cam-Clay Yield Function
  else if ( (strcmp(argv[0],"Cam-Clay") == 0) ) {
    int a1 = -1;  int b1 = 0;
    int a2 = -1;  int b2 = 0;

    if (argc <= 4)  {
        cout << "Warning:  TclNewTemplate3Dep - Yield Function (CC) - input parameters < 4" << endl;
        exit (1);
    }

    if (argc > 4 ) {
      if (Tcl_GetInt(interp, argv[1], &a1) != TCL_OK) {
        cout << "Warning:  TclNewTemplate3Dep - invalid input " << argv[1] << endl;
        exit (1);
      }
      if (Tcl_GetInt(interp, argv[2], &b1) != TCL_OK) {
        cout << "Warning:  TclNewTemplate3Dep - invalid input " << argv[2] << endl;
        exit (1);
      }
      if (Tcl_GetInt(interp, argv[3], &a2) != TCL_OK) {
        cout << "Warning:  TclNewTemplate3Dep - invalid input " << argv[3] << endl;
        exit (1);
      }
      if (Tcl_GetInt(interp, argv[4], &b2) != TCL_OK) {
        cout << "Warning:  TclNewTemplate3Dep - invalid input " << argv[4] << endl;
        exit (1);
      }
    }  

    YF = new CC_YF(a1,b1, a2,b2);
  }

  // 4. Dafalias-Manzari Yield Function
  else if ( strcmp(argv[0],"Dafalias-Manzari") == 0 ) {
    int a1 = -1;  int b1 = 0;
    int a2 = -1;  int b2 = 0;

    if (argc <= 2)  {
        cout << "Warning:  TclNewTemplate3Dep - Yield Function (DM) - input parameters < 2" << endl;
        exit (1);
    }

    if (argc > 2 ) {
      if (Tcl_GetInt(interp, argv[1], &a1) != TCL_OK) {
        cout << "Warning:  TclNewTemplate3Dep - invalid input " << argv[1] << endl;
        exit (1);
      }
      if (Tcl_GetInt(interp, argv[2], &b1) != TCL_OK) {
        cout << "Warning:  TclNewTemplate3Dep - invalid input " << argv[2] << endl;
        exit (1);
      }
      if (Tcl_GetInt(interp, argv[3], &a2) != TCL_OK) {
        cout << "Warning:  TclNewTemplate3Dep - invalid input " << argv[3] << endl;
        exit (1);
      }
      if (Tcl_GetInt(interp, argv[4], &b2) != TCL_OK) {
        cout << "Warning:  TclNewTemplate3Dep - invalid input " << argv[4] << endl;
        exit (1);
      }
    }   

    YF = new DM04_YF(a1,b1, a2,b2);
  }                                            

  // if others
  else {
    cout << "Warning: invalid yield function: " << argv[0] << endl;
    exit (1);
  }

  cleanup(argv);
  return YF;
}


//**************************************************************************************
// Function - to create a Plastic Flow object
PlasticFlow *EvaluatePlasticFlow(ClientData clientData, Tcl_Interp *interp, TCL_Char *tclString)
{
  int argc;
  TCL_Char **argv;

  // split the list
  if (Tcl_SplitList(interp, tclString, &argc, &argv) != TCL_OK) {
    exit (1);
  }

  if (argc == 0)
    exit (1);

  // now parse the list & construct the required object
  PlasticFlow *PF = NULL;

  // 1. von Mises Plastic Flow
  if ((strcmp(argv[0],"Von-Mises") == 0) || (strcmp(argv[0],"von-Mises") == 0)) {
    int a1 = -1;  int b1 = 0;

    //if (argc <= 2)  {
    //    cout << "Warning:  TclNewTemplate3Dep - Plastic Flow (VM) - input parameters < 2" << endl;
    //    exit (1);
    //}

    if (argc > 2 ) {
      if (Tcl_GetInt(interp, argv[1], &a1) != TCL_OK) {
        cout << "Warning:  TclNewTemplate3Dep - invalid input " << argv[1] << endl;
        exit (1);
      }
      if (Tcl_GetInt(interp, argv[2], &b1) != TCL_OK) {
        cout << "Warning:  TclNewTemplate3Dep - invalid input " << argv[2] << endl;
        exit (1);
      }
    } 

    PF = new VM_PF(a1,b1);
  }

  // 2. Drucker-Parger Plastic Flow
  else if ( (strcmp(argv[0],"Drucker-Prager") == 0) ) {
    int a1 = -1;  int b1 = 0;
    int a2 = -1;  int b2 = 0;

    if (argc <= 2)  {
        cout << "Warning:  TclNewTemplate3Dep - Plastic Flow (DP) - input parameters < 2" << endl;
        exit (1);
    }

    if (argc > 2 ) {
      if (Tcl_GetInt(interp, argv[1], &a1) != TCL_OK) {
        cout << "Warning:  TclNewTemplate3Dep - invalid input " << argv[1] << endl;
        exit (1);
      }
      if (Tcl_GetInt(interp, argv[2], &b1) != TCL_OK) {
        cout << "Warning:  TclNewTemplate3Dep - invalid input " << argv[2] << endl;
        exit (1);
      }
    }
    if (argc > 4 ) {
      if (Tcl_GetInt(interp, argv[3], &a2) != TCL_OK) {
        cout << "Warning:  TclNewTemplate3Dep - invalid input " << argv[3] << endl;
        exit (1);
      }
      if (Tcl_GetInt(interp, argv[4], &b2) != TCL_OK) {
        cout << "Warning:  TclNewTemplate3Dep - invalid input " << argv[4] << endl;
        exit (1);
      }
    }   
    
    PF = new DP_PF(a1,b1, a2,b2);
  }

  // 3. Cam-Clay Plastic Flow
  else if ( (strcmp(argv[0],"Cam-Clay") == 0) ) {
    int a1 = -1;  int b1 = 0;
    int a2 = -1;  int b2 = 0;

    if (argc <= 4)  {
        cout << "Warning:  TclNewTemplate3Dep - Plastic Flow (CC) - input parameters < 4" << endl;
        exit (1);
    }

    if (argc > 2 ) {
      if (Tcl_GetInt(interp, argv[1], &a1) != TCL_OK) {
        cout << "Warning:  TclNewTemplate3Dep - invalid input " << argv[1] << endl;
        exit (1);
      }
      if (Tcl_GetInt(interp, argv[2], &b1) != TCL_OK) {
        cout << "Warning:  TclNewTemplate3Dep - invalid input " << argv[2] << endl;
        exit (1);
      }
    }
    if (argc > 4 ) {
      if (Tcl_GetInt(interp, argv[3], &a2) != TCL_OK) {
        cout << "Warning:  TclNewTemplate3Dep - invalid input " << argv[3] << endl;
        exit (1);
      }
      if (Tcl_GetInt(interp, argv[4], &b2) != TCL_OK) {
        cout << "Warning:  TclNewTemplate3Dep - invalid input " << argv[4] << endl;
        exit (1);
      }
    }   

    PF = new CC_PF(a1,b1, a2,b2);
  }

  // 4. Dafalias--Manzari Plastic Flow
  else if ( strcmp(argv[0],"Dafalias-Manzari") == 0 ) {
    int a[12] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};  
    int b[12] = {0,0,0,0,0,0,0,0,0,0,0,0};

    if (argc <= 22)  {
        cout << "Warning:  TclNewTemplate3Dep - Plastic Flow (DM) - input parameters < 22" << endl;
        exit (1);
    }

    if (argc > 24 ) {
      for (int i = 0; i < 12; i++) {
        if (Tcl_GetInt(interp, argv[i*2+1], &a[i]) != TCL_OK) {
          cout << "Warning:  TclNewTemplate3Dep - invalid input " << argv[i*2+1] << endl;
          exit (1);
        }
        if (Tcl_GetInt(interp, argv[i*2+2], &b[i]) != TCL_OK) {
          cout << "Warning:  TclNewTemplate3Dep - invalid input " << argv[i*2+2] << endl;
          exit (1);
        }
      }
    }  

    PF = new DM04_PF( a[0],b[0], a[1],b[1], a[2],b[2], a[3],b[3], a[4],b[4], a[5],b[5], a[6],b[6], 
                      a[7],b[7], a[8],b[8], a[9],b[9], a[10],b[10], a[11],b[11]);
  }

  // if others
  else {
    cout << "Warning: invalid plastic flow: " << argv[0] << endl;
    exit (1);
  }

  cleanup(argv);
  return PF;
}


// Function - to create scalar evolution array
ScalarEvolution** EvaluateSE(ClientData clientData, Tcl_Interp *interp, TCL_Char *tclString)
{
  int argc;
  TCL_Char **argv;

  // split the list
  if (Tcl_SplitList(interp, tclString, &argc, &argv) != TCL_OK) {
    exit (1);
  }

  int loc = 0;
  ScalarEvolution** ISS = NULL;

  if (n_is > 0) {
    ISS = new ScalarEvolution* [n_is];
    if ( ISS == NULL) {
      cout << "Warning:  TclNewTemplate3Dep - invalid input " << endl;
      exit (1);
    }
  }
  else
    return NULL;

     for(int i = 0; i < n_is; i++) {

       // Linear
       if ( strcmp(argv[loc],"Linear") == 0 ) {
         int a1 = 0;          
         loc++; 
         if (Tcl_GetInt(interp, argv[loc], &a1) != TCL_OK) {
           cout << "Warning:  TclNewTemplate3Dep - invalid input " << argv[loc] << endl;
           exit (1);
         }
         ISS[i] = new Linear_Eeq(a1);
         loc++;
       }

       // CC
       else if ( strcmp(argv[0],"Cam-Clay") == 0 ) {
         int a[4] = {0,0,0,0};
         for (int j = 0; j < 4; j++) {
           loc++;
           if (Tcl_GetInt(interp, argv[loc], &a[j]) != TCL_OK) {
             cout << "Warning:  TclNewTemplate3Dep - invalid input " << argv[loc] << endl;
             exit (1);
           }
         }   
         ISS[i] = new CC_Ev(a[0], a[1], a[2], a[4]);
         loc++;
       }

       // other
       else {
         cout << "Warning:  TclNewTemplate3Dep - unknown scalar evolution: " << argv[loc] << endl;
         exit (1);       
       }
           
     }

  cleanup(argv);
  return ISS;
}

// Function - to create tensor evolution array
TensorEvolution** EvaluateTE(ClientData clientData, Tcl_Interp *interp, TCL_Char *tclString)
{
  int argc;
  TCL_Char **argv;

  // split the list
  if (Tcl_SplitList(interp, tclString, &argc, &argv) != TCL_OK) {
    exit (1);
  }

  int loc = 0;
  TensorEvolution** ITT = NULL;

  if (n_it > 0) {
    ITT = new TensorEvolution* [n_is];
    if ( ITT == NULL) {
      cout << "Warning:  TclNewTemplate3Dep - invalid input " << endl;
      exit (1);
    }
  }
  else
    return NULL;

     for(int i = 0; i < n_it; i++) {

       // Linear
       if ( strcmp(argv[loc],"Linear") == 0 ) {
         int a1 = 0;          
         loc++; 
         if (Tcl_GetInt(interp, argv[loc], &a1) != TCL_OK) {
           cout << "Warning:  TclNewTemplate3Dep - invalid input " << argv[loc] << endl;
           exit (1);
         }
         ITT[i] = new Linear_Eij(a1);
         loc++;
       }

       // AF
       else if ( strcmp(argv[0],"Armstrong-Frederick") == 0 ) {
         int a[3] = {0,0,0};
         for (int j = 0; j < 3; j++) {
           loc++;
           if (Tcl_GetInt(interp, argv[loc], &a[j]) != TCL_OK) {
             cout << "Warning:  TclNewTemplate3Dep - invalid input " << argv[loc] << endl;
             exit (1);
           }
         }  
         ITT[i] = new AF_Eij(a[0], a[1], a[2]);
         loc++;
       }

       // DM04-alpha
       else if ( strcmp(argv[loc],"Dafalias-Manzari") == 0 ) {
         int a[14] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0};
         for (int j = 0; j < 14; j++) {
           loc++;
           if (Tcl_GetInt(interp, argv[loc], &a[j]) != TCL_OK) {
             cout << "Warning:  TclNewTemplate3Dep - invalid input " << argv[loc] << endl;
             exit (1);
           }
         }   
         ITT[i] = new DM04_alpha_Eij(a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7], a[8], a[9], a[10], a[11], a[12], a[13]);
         loc++;
       }

       // DM04-z
       else if ( strcmp(argv[loc],"Dafalias-Manzari-fabric") == 0 ) {
         int a[5] = {0,0,0,0,0};
         for (int j = 0; j < 5; j++) {
           loc++;
           if (Tcl_GetInt(interp, argv[loc], &a[j]) != TCL_OK) {
             cout << "Warning:  TclNewTemplate3Dep - invalid input " << argv[loc] << endl;
             exit (1);
           }
         }   
         ITT[i] = new DM04_z_Eij(a[0], a[1], a[2], a[3], a[4]);
         loc++;
       }

       // other
       else {
         cout << "Warning:  TclNewTemplate3Dep - unknown tensor evolution: " << argv[loc] << endl;
         exit (1);              
       }
    
     }

  cleanup(argv);
  return ITT;
}
