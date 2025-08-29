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
// standard library
#include <array>
#include <string>
#include <unordered_map>
// framework
#include <runtimeAPI.h>
#include <Parsing.h>
#include <packages.h>

// system of eqn and solvers
#include <BandSPDLinSOE.h>
#include <BandSPDLinLapackSolver.h>
//
#include <BandGenLinSOE.h>
#include <BandGenLinLapackSolver.h>
//
#include <ConjugateGradientSolver.h>
//
#include <FullGenLinSOE.h>
#include <FullGenLinLapackSolver.h>
//
#include <ProfileSPDLinSOE.h>
#include <ProfileSPDLinDirectSolver.h>
#include <DistributedProfileSPDLinSOE.h>
//
#include <DiagonalSOE.h>
#include <DiagonalDirectSolver.h>
//
#include <SProfileSPDLinSolver.h>
#include <SProfileSPDLinSOE.h>
//
#include <SparseGenColLinSOE.h>
//
#include <SparseGenRowLinSOE.h>
// #include <SymSparseLinSOE.h>
// #include <SymSparseLinSolver.h>
#include <ArpackSOE.h>
#include <ArpackSolver.h>

#if 1 || defined(_PETSC)
LinearSOE *TclCommand_newPetscSOE(int, TCL_Char**);
#endif

#ifdef _CUDA
#  include <BandGenLinSOE_Single.h>
#  include <BandGenLinLapackSolver_Single.h>
#endif

#if defined(_PARALLEL_PROCESSING)
//  parallel soe & solvers
#  include <DistributedBandSPDLinSOE.h>
#  include <DistributedSparseGenColLinSOE.h>
#  include <DistributedSparseGenRowLinSOE.h>
#  include <DistributedBandGenLinSOE.h>
#  include <DistributedDiagonalSOE.h>
#  include <DistributedDiagonalSolver.h>

#  include <MPIDiagonalSOE.h>
#  include <MPIDiagonalSolver.h>

#  define MPIPP_H
#  include <DistributedSuperLU.h>
#  include <DistributedProfileSPDLinSOE.h>
#endif

template <typename T>
using TclDispatch = T(*)(ClientData, Tcl_Interp*, int, const char**);

typedef LinearSOE*(G3_SysOfEqnSpecifier)(G3_Runtime*, int, G3_Char**);

// Specifiers defined in solver.cpp
G3_SysOfEqnSpecifier specify_SparseSPD;
G3_SysOfEqnSpecifier specifySparseGen;
LinearSOE* TclDispatch_newMumpsLinearSOE(ClientData, Tcl_Interp*, Tcl_Size, const char** const);
LinearSOE* TclDispatch_newUmfpackLinearSOE(ClientData, Tcl_Interp*, Tcl_Size, const char** const);
LinearSOE* TclDispatch_newItpackLinearSOE(ClientData, Tcl_Interp*, Tcl_Size, const char** const);

// Helpers to automatically create constructors for systems/solvers 
// that do not take arguments when they are constructed.
template <typename Solver, typename SOE>
LinearSOE *simple_soe(G3_Runtime*, int, G3_Char**) {return new SOE(*(new Solver()));}
#define G3_SOE(Solver, SOE) simple_soe<Solver, SOE>
#define SP_SOE(Solver, SOE) nullptr
#define MP_SOE(Solver, SOE) nullptr

typedef LinearSOE*(*fn)(G3_Runtime*, int, G3_Char**);
struct soefps {fn ss, sp, mp;};

std::unordered_map<std::string, struct soefps> soe_table = {
  {"bandspd", {
     G3_SOE(BandSPDLinLapackSolver,      BandSPDLinSOE),
#if 1
     SP_SOE(BandSPDLinLapackSolver,      BandSPDLinSOE),
     MP_SOE(BandSPDLinLapackSolver,      BandSPDLinSOE)}},
#else
     SP_SOE(BandSPDLinLapackSolver,      DistributedBandSPDLinSOE),
     MP_SOE(BandSPDLinLapackSolver,      DistributedBandSPDLinSOE)}},
#endif
  {"bandgeneral", { // BandGen, BandGEN
     G3_SOE(BandGenLinLapackSolver,      BandGenLinSOE),
     SP_SOE(BandGenLinLapackSolver,      DistributedBandGenLinSOE),
     MP_SOE(BandGenLinLapackSolver,      DistributedBandGenLinSOE)}},
  {"bandgen", { // BandGen, BandGEN
     G3_SOE(BandGenLinLapackSolver,      BandGenLinSOE),
     SP_SOE(BandGenLinLapackSolver,      DistributedBandGenLinSOE),
     MP_SOE(BandGenLinLapackSolver,      DistributedBandGenLinSOE)}},

  {"sparsegen",     {specifySparseGen, nullptr, nullptr}},
  {"sparsegeneral", {specifySparseGen, nullptr, nullptr}},
  {"superlu",       {specifySparseGen, nullptr, nullptr}},

  {"sparsesym", {
     specify_SparseSPD, nullptr, nullptr}},

  {"sparsespd", {
     // Legacy specifier
     specify_SparseSPD, nullptr, nullptr}},

  {"diagonal", {
     G3_SOE(DiagonalDirectSolver,        DiagonalSOE),
     SP_SOE(DistributedDiagonalSolver,   DistributedDiagonalSOE),
     MP_SOE(DistributedDiagonalSolver,   DistributedDiagonalSOE)}},

  {"mpidiagonal", {
     G3_SOE(DiagonalDirectSolver,        DiagonalSOE),
     SP_SOE(MPIDiagonalSolver,           MPIDiagonalSOE),
     MP_SOE(MPIDiagonalSolver,           MPIDiagonalSOE)}},

  {"sprofilespd", {
     G3_SOE(SProfileSPDLinSolver,        SProfileSPDLinSOE),
     SP_SOE(SProfileSPDLinSolver,        SProfileSPDLinSOE),
     MP_SOE(SProfileSPDLinSolver,        SProfileSPDLinSOE)}},

  {"profilespd", {
     G3_SOE(ProfileSPDLinDirectSolver,   ProfileSPDLinSOE),
#ifndef PARALLEL_PROFILESPD
     SP_SOE(ProfileSPDLinDirectSolver,   ProfileSPDLinSOE),
     MP_SOE(ProfileSPDLinDirectSolver,   ProfileSPDLinSOE)}},
#else
     SP_SOE(ProfileSPDLinDirectSolver,   DistributedProfileSPDLinSOE),
     MP_SOE(ProfileSPDLinDirectSolver,   DistributedProfileSPDLinSOE)}},

  {"parallelprofilespd", {
     nullptr, nullptr,
     MP_SOE(ProfileSPDLinDirectSolver,   DistributedProfileSPDLinSOE)}},
#endif
  {"fullgeneral", {
     G3_SOE(FullGenLinLapackSolver,      FullGenLinSOE),
     SP_SOE(FullGenLinLapackSolver,      FullGenLinSOE),
     MP_SOE(FullGenLinLapackSolver,      FullGenLinSOE)}},
  {"fullgen", {
     G3_SOE(FullGenLinLapackSolver,      FullGenLinSOE),
     SP_SOE(FullGenLinLapackSolver,      FullGenLinSOE),
     MP_SOE(FullGenLinLapackSolver,      FullGenLinSOE)}},

#ifdef _CUDA
  {"BandGeneral_Single", {  // "BandGEN_Single", "BandGen_Single"
     G3_SOE(BandGenLinLapackSolver_Single,    BandGenLinSOE_Single),
     SP_SOE(BandGenLinLapackSolver_Single,    BandGenLinSOE_Single),
     MP_SOE(BandGenLinLapackSolver_Single,    BandGenLinSOE_Single)}}
  }
#endif
};

