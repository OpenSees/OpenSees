/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */

// Written: Minjie

// Description: A tcl wrapper for OpenSees commands
//

#include "TclWrapper.h"
#include "OpenSeesCommands.h"
#include <OPS_Globals.h>

static TclWrapper* wrapper = 0;

TclWrapper::TclWrapper()
    :currentArgv(0), currentArg(0), numberArgs(0)
{
    wrapper = this;
}

TclWrapper::~TclWrapper()
{
}

void
TclWrapper::resetCommandLine(int nArgs, int cArg, TCL_Char** argv)
{
    numberArgs = nArgs;
    currentArg = cArg;
    currentArgv = argv;
}

void
TclWrapper::resetCommandLine(int cArg)
{
    if (cArg < 0) {
	currentArg += cArg;
	if (currentArg < 0) currentArg = 0;
    } else {
	currentArg = cArg;
    }
}

void
TclWrapper::addCommand(Tcl_Interp* interp, const char* name, Tcl_CmdProc* proc)
{
    Tcl_CreateCommand(interp, name, proc,NULL,NULL);
}

void
TclWrapper::setOutputs(Tcl_Interp* interp, int* data, int numArgs)
{
    char buffer[40];
    for (int i=0; i<numArgs; i++) {
	sprintf(buffer, "%d ", data[i]);
	Tcl_AppendResult(interp, buffer, NULL);
    }
}

void
TclWrapper::setOutputs(Tcl_Interp* interp, double* data, int numArgs)
{
    char buffer[40];
    for (int i=0; i<numArgs; i++) {
	sprintf(buffer, "%35.20f ", data[i]);
	Tcl_AppendResult(interp, buffer, NULL);
    }
}

void
TclWrapper::setOutputs(Tcl_Interp* interp, const char* str)
{
    Tcl_SetResult(interp, (char*)str, TCL_VOLATILE);
}

///////////////////////////////////////////
/////// Tcl wrapper functions  ////////////
///////////////////////////////////////////
static int Tcl_ops_UniaxialMaterial(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_UniaxialMaterial() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_testUniaxialMaterial(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_testUniaxialMaterial() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_setStrain(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_setStrain() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_getStrain(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_getStrain() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_getStress(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_getStress() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_getTangent(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_getTangent() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_getDampTangent(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_getDampTangent() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_wipe(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_wipe() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_model(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_model() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_node(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_Node() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_fix(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_HomogeneousBC() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_element(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_Element() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_timeSeries(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_TimeSeries() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_pattern(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_Pattern() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_nodalLoad(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_NodalLoad() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_system(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_System() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_numberer(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_Numberer() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_constraints(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_ConstraintHandler() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_integrator(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_Integrator() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_algorithm(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_Algorithm() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_analysis(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_Analysis() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_analyze(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_analyze() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_nodeDisp(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_nodeDisp() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_test(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_CTest() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_section(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_Section() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_fiber(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_Fiber() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_patch(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_Patch() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_layer(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_Layer() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_geomTransf(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_CrdTransf() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_beamIntegration(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_BeamIntegration() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_loadConst(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_loadConst() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_eleLoad(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_ElementalLoad() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_reactions(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_calculateNodalReactions() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_nodeReaction(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_nodeReaction() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_eigen(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_eigenAnalysis() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_nDMaterial(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_NDMaterial() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_block2d(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_doBlock2D() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_block3d(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_doBlock3D() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_rayleigh(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_rayleighDamping() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_wipeAnalysis(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_wipeAnalysis() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_setTime(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_setTime() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_remove(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_removeObject() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_mass(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_addNodalMass() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_equalDOF(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_EqualDOF() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_nodeEigenvector(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_nodeEigenvector() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_getTime(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_getTime() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_eleResponse(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_eleResponse() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_SP(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_SP() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_fixX(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_HomogeneousBC_X() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_fixY(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_HomogeneousBC_Y() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_fixZ(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_HomogeneousBC_Z() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_reset(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_resetModel() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_initialize(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_initializeAnalysis() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_getLoadFactor(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_getLoadFactor() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_build(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_buildModel() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_print(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_printModel() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_printA(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_printA() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_printB(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_printB() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_printGID(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_printModelGID() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_getCTestNorms(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_getCTestNorms() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_getCTestIter(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_getCTestIter() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_recorder(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_Recorder() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_database(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_Database() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_record(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_record() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_save(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_save() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_restore(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_restore() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_eleForce(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_eleForce() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_eleDynamicalForce(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_eleDynamicalForce() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_nodeUnbalance(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_nodeUnbalance() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_nodeVel(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_nodeVel() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_setNodeVel(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_setNodeVel() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_setNodeDisp(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_setNodeDisp() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_nodeAccel(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_nodeAccel() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_setNodeAccel(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_setNodeAccel() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_nodeResponse(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_nodeResponse() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_nodeCoord(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_nodeCoord() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_setNodeCoord(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_setNodeCoord() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_updateElementDomain(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_updateElementDomain() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_eleNodes(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_eleNodes() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_nodeDOFs(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_nodeDOFs() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_nodeMass(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_nodeMass() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_nodePressure(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_nodePressure() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_nodeBounds(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_nodeBounds() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_startTimer(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_startTimer() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_stopTimer(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_stopTimer() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_modalDamping(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_modalDamping() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_modalDampingQ(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_modalDampingQ() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_setElementRayleighDampingFactors(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_setElementRayleighDampingFactors() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_region(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_MeshRegion() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_setPrecision(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_setPrecision() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_searchPeerNGA(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_peerNGA() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_domainChange(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_domainChange() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_metaData(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_neesMetaData() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_defaultUnits(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_defaultUnits() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_neesUpload(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_neesUpload() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_stripXML(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_stripOpenSeesXML() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_convertBinaryToText(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_convertBinaryToText() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_convertTextToBinary(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_convertTextToBinary() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_getEleTags(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_getEleTags() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_getNodeTags(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_getNodeTags() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_getParamTags(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_getParamTags() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_getParamValue(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_getParamValue() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_sectionForce(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_sectionForce() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_sectionDeformation(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_sectionDeformation() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_sectionStiffness(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_sectionStiffness() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_sectionFlexibility(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_sectionFlexibility() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_sectionLocation(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_sectionLocation() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_sectionWeight(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_sectionWeight() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_basicDeformation(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_basicDeformation() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_basicForce(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_basicForce() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_basicStiffness(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_basicStiffness() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_InitialStateAnalysis(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_InitialStateAnalysis() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_totalCPU(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_totalCPU() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_solveCPU(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_solveCPU() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_accelCPU(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_accelCPU() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_numFact(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_numFact() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_numIter(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_numIter() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_systemSize(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_systemSize() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_version(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_version() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_setMaxOpenFiles(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_maxOpenFiles() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_limitCurve(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_LimitCurve() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_imposedMotion(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_ImposedMotionSP() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_groundMotion(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_groundMotion() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_equalDOF_Mixed(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_EqualDOF_Mixed() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_rigidLink(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_RigidLink() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_rigidDiaphragm(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_RigidDiaphragm() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_ShallowFoundationGen(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_ShallowFoundationGen() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_setElementRayleighFactors(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_addElementRayleigh() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_mesh(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_mesh() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_remesh(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_remesh() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_parameter(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_Parameter() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_addToParameter(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_addToParameter() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_updateParameter(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_updateParameter() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_setParameter(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv) {
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_setParameter() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_getPID(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv)
{
    wrapper->resetCommandLine(argc, 1, argv);

    //if (OPS_getPID() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_getNP(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv)
{
    wrapper->resetCommandLine(argc, 1, argv);

    //if (OPS_getNP() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_barrier(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv)
{
    wrapper->resetCommandLine(argc, 1, argv);

    //if (OPS_barrier() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_send(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv)
{
    wrapper->resetCommandLine(argc, 1, argv);

    //if (OPS_send() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_recv(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv)
{
    wrapper->resetCommandLine(argc, 1, argv);

    //if (OPS_recv() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_frictionModel(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv)
{
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_FrictionModel() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_computeGradients(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv)
{
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_computeGradients() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_sensitivityAlgorithm(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv)
{
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_sensitivityAlgorithm() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_sensNodeDisp(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv)
{
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_sensNodeDisp() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_sensNodeVel(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv)
{
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_sensNodeVel() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_sensNodeAccel(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv)
{
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_sensNodeAccel() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_sensLambda(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv)
{
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_sensLambda() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_sensSectionForce(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv)
{
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_sensSectionForce() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_sensNodePressure(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv)
{
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_sensNodePressure() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_randomVariable(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv)
{
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_randomVariable() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_getRVTags(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv)
{
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_getRVTags() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_getRVMean(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv)
{
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_getRVMean() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_getRVStdv(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv)
{
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_getRVStdv() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_getRVPDF(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv)
{
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_getRVPDF() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_getRVCDF(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv)
{
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_getRVCDF() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_getRVInverseCDF(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv)
{
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_getRVInverseCDF() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_addCorrelate(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv)
{
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_addCorrelate() < 0) return TCL_ERROR;

    return TCL_OK;
}
static int Tcl_ops_updateMaterialStage(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv)
{
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_updateMaterialStage() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_sdfResponse(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv)
{
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_sdfResponse() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_getNumThreads(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv)
{
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_getNumThreads() < 0) return TCL_ERROR;

    return TCL_OK;
}

static int Tcl_ops_setNumThreads(ClientData clientData, Tcl_Interp *interp, int argc,   TCL_Char **argv)
{
    wrapper->resetCommandLine(argc, 1, argv);

    if (OPS_setNumThreads() < 0) return TCL_ERROR;

    return TCL_OK;
}

//////////////////////////////////////////////
////////////// Add Tcl commands //////////////
//////////////////////////////////////////////
void
TclWrapper::addOpenSeesCommands(Tcl_Interp* interp)
{
    addCommand(interp,"uniaxialMaterial", &Tcl_ops_UniaxialMaterial);
    addCommand(interp,"testUniaxialMaterial", &Tcl_ops_testUniaxialMaterial);
    addCommand(interp,"setStrain", &Tcl_ops_setStrain);
    addCommand(interp,"getStrain", &Tcl_ops_getStrain);
    addCommand(interp,"getStress", &Tcl_ops_getStress);
    addCommand(interp,"getTangent", &Tcl_ops_getTangent);
    addCommand(interp,"getDampTangent", &Tcl_ops_getDampTangent);
    addCommand(interp,"wipe", &Tcl_ops_wipe);
    addCommand(interp,"model", &Tcl_ops_model);
    addCommand(interp,"node", &Tcl_ops_node);
    addCommand(interp,"fix", &Tcl_ops_fix);
    addCommand(interp,"element", &Tcl_ops_element);
    addCommand(interp,"timeSeries", &Tcl_ops_timeSeries);
    addCommand(interp,"pattern", &Tcl_ops_pattern);
    addCommand(interp,"load", &Tcl_ops_nodalLoad);
    addCommand(interp,"system", &Tcl_ops_system);
    addCommand(interp,"numberer", &Tcl_ops_numberer);
    addCommand(interp,"constraints", &Tcl_ops_constraints);
    addCommand(interp,"integrator", &Tcl_ops_integrator);
    addCommand(interp,"algorithm", &Tcl_ops_algorithm);
    addCommand(interp,"analysis", &Tcl_ops_analysis);
    addCommand(interp,"analyze", &Tcl_ops_analyze);
    addCommand(interp,"test", &Tcl_ops_test);
    addCommand(interp,"section", &Tcl_ops_section);
    addCommand(interp,"fiber", &Tcl_ops_fiber);
    addCommand(interp,"patch", &Tcl_ops_patch);
    addCommand(interp,"layer", &Tcl_ops_layer);
    addCommand(interp,"geomTransf", &Tcl_ops_geomTransf);
    addCommand(interp,"beamIntegration", &Tcl_ops_beamIntegration);
    addCommand(interp,"loadConst", &Tcl_ops_loadConst);
    addCommand(interp,"eleLoad", &Tcl_ops_eleLoad);
    addCommand(interp,"reactions", &Tcl_ops_reactions);
    addCommand(interp,"nodeReaction", &Tcl_ops_nodeReaction);
    addCommand(interp,"eigen", &Tcl_ops_eigen);
    addCommand(interp,"nDMaterial", &Tcl_ops_nDMaterial);
    addCommand(interp,"block2D", &Tcl_ops_block2d);
    addCommand(interp,"block3D", &Tcl_ops_block3d);
    addCommand(interp,"rayleigh", &Tcl_ops_rayleigh);
    addCommand(interp,"wipeAnalysis", &Tcl_ops_wipeAnalysis);
    addCommand(interp,"setTime", &Tcl_ops_setTime);
    addCommand(interp,"remove", &Tcl_ops_remove);
    addCommand(interp,"mass", &Tcl_ops_mass);
    addCommand(interp,"equalDOF", &Tcl_ops_equalDOF);
    addCommand(interp,"nodeEigenvector", &Tcl_ops_nodeEigenvector);
    addCommand(interp,"getTime", &Tcl_ops_getTime);
    addCommand(interp,"eleResponse", &Tcl_ops_eleResponse);
    addCommand(interp,"sp", &Tcl_ops_SP);
    addCommand(interp,"fixX", &Tcl_ops_fixX);
    addCommand(interp,"fixY", &Tcl_ops_fixY);
    addCommand(interp,"fixZ", &Tcl_ops_fixZ);
    addCommand(interp,"reset", &Tcl_ops_reset);
    addCommand(interp,"initialize", &Tcl_ops_initialize);
    addCommand(interp,"getLoadFactor", &Tcl_ops_getLoadFactor);
    addCommand(interp,"build", &Tcl_ops_build);
    addCommand(interp,"printModel", &Tcl_ops_print);
    addCommand(interp,"printA", &Tcl_ops_printA);
    addCommand(interp,"printB", &Tcl_ops_printB);
    addCommand(interp,"printGID", &Tcl_ops_printGID);
    addCommand(interp,"getCTestNorms", &Tcl_ops_getCTestNorms);
    addCommand(interp,"getCTestIter", &Tcl_ops_getCTestIter);
    addCommand(interp,"recorder", &Tcl_ops_recorder);
    addCommand(interp,"database", &Tcl_ops_database);
    addCommand(interp,"save", &Tcl_ops_save);
    addCommand(interp,"restore", &Tcl_ops_restore);
    addCommand(interp,"eleForce", &Tcl_ops_eleForce);
    addCommand(interp,"eleDynamicalForce", &Tcl_ops_eleDynamicalForce);
    addCommand(interp,"nodeUnbalance", &Tcl_ops_nodeUnbalance);
    addCommand(interp,"nodeDisp", &Tcl_ops_nodeDisp);
    addCommand(interp,"setNodeDisp", &Tcl_ops_setNodeDisp);
    addCommand(interp,"nodeVel", &Tcl_ops_nodeVel);
    addCommand(interp,"setNodeVel", &Tcl_ops_setNodeVel);
    addCommand(interp,"nodeAccel", &Tcl_ops_nodeAccel);
    addCommand(interp,"setNodeAccel", &Tcl_ops_setNodeAccel);
    addCommand(interp,"nodeResponse", &Tcl_ops_nodeResponse);
    addCommand(interp,"nodeCoord", &Tcl_ops_nodeCoord);
    addCommand(interp,"setNodeCoord", &Tcl_ops_setNodeCoord);
    addCommand(interp,"updateElementDomain", &Tcl_ops_updateElementDomain);
    addCommand(interp,"eleNodes", &Tcl_ops_eleNodes);
    addCommand(interp,"nodeMass", &Tcl_ops_nodeMass);
    addCommand(interp,"nodeDOFs", &Tcl_ops_nodeDOFs);
    addCommand(interp,"nodePressure", &Tcl_ops_nodePressure);
    addCommand(interp,"nodeBounds", &Tcl_ops_nodeBounds);
    addCommand(interp,"start", &Tcl_ops_startTimer);
    addCommand(interp,"stop", &Tcl_ops_stopTimer);
    addCommand(interp,"modalDamping", &Tcl_ops_modalDamping);
    addCommand(interp,"modalDampingQ", &Tcl_ops_modalDampingQ);
    addCommand(interp,"setElementRayleighDampingFactors", &Tcl_ops_setElementRayleighDampingFactors);
    addCommand(interp,"region", &Tcl_ops_region);
    addCommand(interp,"setPrecision", &Tcl_ops_setPrecision);
    addCommand(interp,"searchPeerNGA", &Tcl_ops_searchPeerNGA);
    addCommand(interp,"domainChange", &Tcl_ops_domainChange);
    addCommand(interp,"record", &Tcl_ops_record);
    addCommand(interp,"metaData", &Tcl_ops_metaData);
    addCommand(interp,"defaultUnits", &Tcl_ops_defaultUnits);
    addCommand(interp,"neesUpload", &Tcl_ops_neesUpload);
    addCommand(interp,"stripXML", &Tcl_ops_stripXML);
    addCommand(interp,"convertBinaryToText", &Tcl_ops_convertBinaryToText);
    addCommand(interp,"convertTextToBinary", &Tcl_ops_convertTextToBinary);
    addCommand(interp,"getEleTags", &Tcl_ops_getEleTags);
    addCommand(interp,"getNodeTags", &Tcl_ops_getNodeTags);
    addCommand(interp,"getParamTags", &Tcl_ops_getParamTags);
    addCommand(interp,"getParamValue", &Tcl_ops_getParamValue);
    addCommand(interp,"sectionForce", &Tcl_ops_sectionForce);
    addCommand(interp,"sectionDeformation", &Tcl_ops_sectionDeformation);
    addCommand(interp,"sectionStiffness", &Tcl_ops_sectionStiffness);
    addCommand(interp,"sectionFlexibility", &Tcl_ops_sectionFlexibility);
    addCommand(interp,"sectionLocation", &Tcl_ops_sectionLocation);
    addCommand(interp,"sectionWeight", &Tcl_ops_sectionWeight);
    addCommand(interp,"basicDeformation", &Tcl_ops_basicDeformation);
    addCommand(interp,"basicForce", &Tcl_ops_basicForce);
    addCommand(interp,"basicStiffness", &Tcl_ops_basicStiffness);
    addCommand(interp,"InitialStateAnalysis", &Tcl_ops_InitialStateAnalysis);
    addCommand(interp,"totalCPU", &Tcl_ops_totalCPU);
    addCommand(interp,"solveCPU", &Tcl_ops_solveCPU);
    addCommand(interp,"accelCPU", &Tcl_ops_accelCPU);
    addCommand(interp,"numFact", &Tcl_ops_numFact);
    addCommand(interp,"numIter", &Tcl_ops_numIter);
    addCommand(interp,"systemSize", &Tcl_ops_systemSize);
    addCommand(interp,"version", &Tcl_ops_version);
    addCommand(interp,"setMaxOpenFiles", &Tcl_ops_setMaxOpenFiles);
    addCommand(interp,"limitCurve", &Tcl_ops_limitCurve);
    addCommand(interp,"imposedMotion", &Tcl_ops_imposedMotion);
    addCommand(interp,"imposedSupportMotion", &Tcl_ops_imposedMotion);
    addCommand(interp,"groundMotion", &Tcl_ops_groundMotion);
    addCommand(interp,"equalDOF_Mixed", &Tcl_ops_equalDOF_Mixed);
    addCommand(interp,"rigidLink", &Tcl_ops_rigidLink);
    addCommand(interp,"rigidDiaphragm", &Tcl_ops_rigidDiaphragm);
    addCommand(interp,"ShallowFoundationGen", &Tcl_ops_ShallowFoundationGen);
    addCommand(interp,"setElementRayleighFactors", &Tcl_ops_setElementRayleighFactors);
    addCommand(interp,"mesh", &Tcl_ops_mesh);
    addCommand(interp,"remesh", &Tcl_ops_remesh);
    addCommand(interp,"parameter", &Tcl_ops_parameter);
    addCommand(interp,"addToParameter", &Tcl_ops_addToParameter);
    addCommand(interp,"updateParameter", &Tcl_ops_updateParameter);
    addCommand(interp,"setParameter", &Tcl_ops_setParameter);
    addCommand(interp,"getPID", &Tcl_ops_getPID);
    addCommand(interp,"getNP", &Tcl_ops_getNP);
    addCommand(interp,"barrier", &Tcl_ops_barrier);
    addCommand(interp,"send", &Tcl_ops_send);
    addCommand(interp,"recv", &Tcl_ops_recv);
    addCommand(interp,"frictionModel", &Tcl_ops_frictionModel);
    addCommand(interp,"computeGradients", &Tcl_ops_computeGradients);
    addCommand(interp,"sensitivityAlgorithm", &Tcl_ops_sensitivityAlgorithm);
    addCommand(interp,"sensNodeDisp", &Tcl_ops_sensNodeDisp);
    addCommand(interp,"sensNodeVel", &Tcl_ops_sensNodeVel);
    addCommand(interp,"sensNodeAccel", &Tcl_ops_sensNodeAccel);
    addCommand(interp,"sensLambda", &Tcl_ops_sensLambda);
    addCommand(interp,"sensSectionForce", &Tcl_ops_sensSectionForce);
    addCommand(interp,"sensNodePressure", &Tcl_ops_sensNodePressure);
    addCommand(interp,"randomVariable", &Tcl_ops_randomVariable);
    addCommand(interp,"updateMaterialStage", &Tcl_ops_updateMaterialStage);
    addCommand(interp,"sdfResponse", &Tcl_ops_sdfResponse);
    addCommand(interp,"getNumThreads", &Tcl_ops_getNumThreads);
    addCommand(interp,"setNumThreads", &Tcl_ops_setNumThreads);

