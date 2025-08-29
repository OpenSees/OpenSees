//===----------------------------------------------------------------------===//
//
//                                   xara
//                              https://xara.so
//
//===----------------------------------------------------------------------===//
//
// Copyright (c) 2025, Claudio M. Perez
// All rights reserved.  No warranty, explicit or implicit, is provided.
//
// This source code is licensed under the BSD 2-Clause License.
// See LICENSE file or https://opensource.org/licenses/BSD-2-Clause
//
//===----------------------------------------------------------------------===//
//
// Written: rms, MHS, cmp
// Created: 07/99
//
// Description: This file contains the function invoked when the user invokes
// the section command in the interpreter.
//
//
#include <tcl.h>
#include <string.h>

#include <Parsing.h>
#include <Logging.h>
#include <BasicModelBuilder.h>
#include <PlateFiberMaterial.h>
#include <ElasticMembranePlateSection.h>
#include <MembranePlateFiberSection.h>
#include <LayeredShellFiberSection.h>

typedef SectionForceDeformation ShellSection;


int
TclCommand_addElasticShellSection(ClientData clientData, Tcl_Interp* interp,
                                  int argc, TCL_Char** const argv)
{

    BasicModelBuilder* builder = static_cast<BasicModelBuilder*>(clientData);

    if (argc < 5) {
      opserr << OpenSees::PromptValueError
             << "insufficient arguments\n";
      opserr << "Want: section ElasticMembranePlateSection tag? E? nu? h? "
                "<rho?> <Ep_mod?>"
             << endln;
      return TCL_ERROR;
    }

    int tag;
    double E, nu, h;
    double rho = 0.0;
    double Ep_mod = 1.0;

    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << OpenSees::PromptValueError 
             << "invalid section ElasticMembranePlateSection tag"
             << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[3], &E) != TCL_OK) {
      opserr << OpenSees::PromptValueError 
             << "invalid E" << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[4], &nu) != TCL_OK) {
      opserr << OpenSees::PromptValueError 
             << "invalid nu" << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[5], &h) != TCL_OK) {
      opserr << OpenSees::PromptValueError 
             << "invalid h" << endln;
      return TCL_ERROR;
    }

    if (argc > 6 && Tcl_GetDouble(interp, argv[6], &rho) != TCL_OK) {
      opserr << OpenSees::PromptValueError 
             << "invalid rho" << endln;
      return TCL_ERROR;
    }

    if (argc > 7 && Tcl_GetDouble(interp, argv[7], &Ep_mod) != TCL_OK) {
      opserr << OpenSees::PromptValueError 
             << "invalid Ep_mod" << endln;
      return TCL_ERROR;
    }

    builder->addTaggedObject<SectionForceDeformation>(*new ElasticMembranePlateSection(tag, E, nu, h, rho, Ep_mod));
    return TCL_OK;
}


int
TclCommand_ShellSection(ClientData clientData, Tcl_Interp* interp, 
                        Tcl_Size argc, TCL_Char** const argv)
{
  // Pointer to a section that will be added to the model builder
  SectionForceDeformation* theSection = nullptr;
  BasicModelBuilder *builder = static_cast<BasicModelBuilder*>(clientData);

  if ((strcmp(argv[1], "PlateFiber") == 0) ||
      (strcmp(argv[1], "PlateFiberThermal") == 0)) { // TODO: add thermal
    
    if (argc < 5) {
      opserr << "WARNING insufficient arguments\n";
      opserr << "Want: section PlateFiber tag? matTag? h? \n";
      return TCL_ERROR;
    }

    int tag;
    double h;

    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid section PlateFiber tag\n";
      return TCL_ERROR;
    }

    int matTag;
    if (Tcl_GetInt(interp, argv[3], &matTag) != TCL_OK) {
      opserr << "WARNING invalid material tag " << argv[3] << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[4], &h) != TCL_OK) {
      opserr << "WARNING invalid h\n";
      return TCL_ERROR;
    }

    NDMaterial* theMaterial = builder->getTypedObject<NDMaterial>(matTag);
    if (theMaterial == nullptr)
      return TCL_ERROR;

    theSection = new MembranePlateFiberSection(tag, h, *theMaterial);
  }

  else if ((strcmp(argv[1], "LayeredShell") == 0) ||
           (strcmp(argv[1], "LayeredShellThermal") == 0)) { // TODO: add thermal

    // section LayeredShell tag? nLayers? <mat1? h1? ... matTagn? hn?> -or- <matTag? thickness?> 

    int status = TCL_ERROR;
    if (argc < 6) {
      opserr << OpenSees::PromptValueError 
             << "insufficient arguments " << "\n";
      opserr << "Want: section LayeredShell tag? nLayers? mat1? h1? ... matn? hn? "
             << endln;
      return TCL_ERROR;
    }

    int tag, nLayers; // , matTag;
    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid section tag" << "\n";
      return TCL_ERROR;
    }

    double h, *thickness;
    NDMaterial **theMats;
    if (Tcl_GetInt(interp, argv[3], &nLayers) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid nLayers" << "\n";
      opserr << "LayeredShell section: " << tag << "\n";
      return TCL_ERROR;
    }

    if (nLayers < 3) {
      opserr << "ERROR number of layers must be larger than 2" << endln;
      return TCL_ERROR;
    }

    if (argc < 3+2*nLayers) {
      opserr << OpenSees::PromptValueError << "Must provide " << 2*nLayers << " layers\n";
      return TCL_ERROR;
    }

    theMats = new NDMaterial *[nLayers];
    thickness = new double[nLayers];

    for (int iLayer = 0; iLayer < nLayers; iLayer++) {
      int mat;
      if (Tcl_GetInt(interp, argv[4 + 2 * iLayer], &mat) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "invalid material tag" << endln;
        status = TCL_ERROR;
        goto cleanup;
      }

      NDMaterial* material = builder->getTypedObject<NDMaterial>(mat);
      if (material == nullptr) {
        status = TCL_ERROR;
        goto cleanup;
      }

      theMats[iLayer] = material->getCopy("PlateFiber");
      if (theMats[iLayer] == nullptr) {
        theMats[iLayer] = new PlateFiberMaterial(mat, *material);
      }

      if (Tcl_GetDouble(interp, argv[5 + 2 * iLayer], &h) != TCL_OK) {
        opserr << OpenSees::PromptValueError 
               << "invalid h at layer " << iLayer << "\n";
        status = TCL_ERROR;
        goto cleanup;
      }

      if (h <= 0) {
        opserr << OpenSees::PromptValueError 
               << "invalid h at layer " << iLayer << "\n";
        status = TCL_ERROR;
        goto cleanup;
      }

      thickness[iLayer] = h;
    }

    theSection = new LayeredShellFiberSection(tag, nLayers, thickness, theMats);


    if (builder->addTaggedObject<ShellSection>(*theSection) == TCL_OK) {
      status = TCL_OK;
    }

cleanup:
    if (thickness != nullptr)
      delete[] thickness;

    if (theMats != 0) {
      for (int iLayer = 0; iLayer < nLayers; iLayer++) {
        if (theMats[iLayer] != nullptr)
          delete theMats[iLayer];
      }
      delete[] theMats;
    }
    return status;
  }


  // Now add the material to the modelBuilder
  if (builder->addTaggedObject<ShellSection>(*theSection) < 0) {
    opserr << "WARNING could not add section to the domain\n";
    opserr << *theSection << "\n";
    delete theSection;
    return TCL_ERROR;
  }

  return TCL_OK;
}

