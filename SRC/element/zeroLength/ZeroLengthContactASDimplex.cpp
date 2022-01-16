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

// $Source: /usr/local/cvs/OpenSees/SRC/element/zeroLength/ZeroLengthContactASDimplex.cpp,v $
// $Revision: 1.0 $
// $Date: 2020-May $

// Written: Onur Deniz Akan         (onur.akan@iusspavia.it)
//          Dr. Massimo Petracca    
//          Prof. Guido Camata      
//          Prof. Enrico Spacone
//          Prof. Carlo G. Lai
//
// Created: May 2020
//
// Description: This file contains the implementation for the ZeroLengthContactASDimplex class.

#include "ZeroLengthContactASDimplex.h"
#include <Information.h>
#include <Domain.h>
#include <Node.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Renderer.h>
#include <limits>
#include <algorithm>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <ElementResponse.h>
#include <elementAPI.h>
#include <map>

namespace
{
    class GlobalStorage
    {
    public:
        int size = 0;
        Matrix K; // stiffness
        Matrix K0; // initial stiffness
        Matrix M; // mass 
        Matrix D; // damping
        Vector U; // displacement
        Vector R; // residual

    public:
        GlobalStorage() = default;
        GlobalStorage& resize(int N) {
            if (N != size) {
                K.resize(N, N);
                K0.resize(N, N);
                M.resize(N, N);
                D.resize(N, N);
                U.resize(N);
                R.resize(N);
            }
            return *this;
        }
    };

    static GlobalStorage& getGlobalStorage(int N)
    {
        static std::map<int, GlobalStorage> gsmap;
        return gsmap[N].resize(N);
    }

    inline void cross(const Vector& A, const Vector& B, Vector& C) {

        C(0) = A(1) * B(2) - A(2) * B(1);
        C(1) = A(2) * B(0) - A(0) * B(2);
        C(2) = A(0) * B(1) - A(1) * B(0);
    }
}

void * OPS_ADD_RUNTIME_VPV(OPS_ZeroLengthContactASDimplex) {

    double SmallNumber = 1.0e-6;
    Element* theElement = nullptr;

    // some kudos
    static int counter = 0;
    if (++counter == 1)
        opserr << "ZeroLengthContactASDimplex element - Implemented: Akan, OD., Petracca, M., Camata, G., Spacone, E. & Lai, CG. (2020)\n";

    // model dimension
    int ndm = OPS_GetNDM();
    if (ndm < 2 || ndm > 3) {
        opserr << "ZeroLengthContactASDimplex: Unsupported NDM (" << ndm << "). It should be 2 or 3\n";
        return theElement;
    }

    // a quick check on number of args
    if (OPS_GetNumRemainingInputArgs() < 6) {
        opserr << "ZeroLengthContactASDimplex: WARNING: too few arguments \n" <<
            "want - element zeroLengthContactASDimplex eleTag? iNode? jNode? Kn? Kt? mu? <-orient $x1 $x2 $x3> <-intType type?>\n";
        return theElement;
    }

    // start with mandatory inputs
    // read eleTag, iNode, jNode
    int idata[3];
    int numdata = 3;
    if (OPS_GetIntInput(&numdata, idata) < 0) {
        opserr << "ZeroLengthContactASDimplex: WARNING: invalid int inputs\n";
        return theElement;
    }

    // read Kn, Kt, mu
    double ddata[3];
    numdata = 3;
    if (OPS_GetDoubleInput(&numdata, ddata) < 0) {
        opserr << "ZeroLengthContactASDimplex: WARNING: invalid double inputs\n";
        return theElement;
    }

    // continue with optional inputs
    Vector x_e(3); x_e(0) = 1.0; x_e(1) = 0.0; x_e(2) = 0.0;                        // initialize orientation vector
    int integrationType = 0;                                                        // implicit by default
    while (OPS_GetNumRemainingInputArgs() > 0) {
        const char* inputstring = OPS_GetString();
        if (strcmp(inputstring, "-orient") == 0) {                                  // #1 read global element orientation
            if (ndm == 2) {
                if (OPS_GetNumRemainingInputArgs() < 2) {
                    opserr << "ZeroLengthContactASDimplex: WARNING: insufficient orient values in 2D\n";
                    return theElement;
                }
                numdata = 3;
                if (OPS_GetDoubleInput(&numdata, &x_e(0)) < 0) {
                    opserr << "ZeroLengthContactASDimplex: WARNING: invalid double input after -orient\n";
                    return theElement;
                }
            }
            else if (ndm == 3) {
                if (OPS_GetNumRemainingInputArgs() < 3) {
                    opserr << "ZeroLengthContactASDimplex: WARNING: insufficient orient values in 3D\n";
                    return theElement;
                }
                numdata = 3;
                if (OPS_GetDoubleInput(&numdata, &x_e(0)) < 0) {
                    opserr << "ZeroLengthContactASDimplex: WARNING: invalid double input after -orient\n";
                    return theElement;
                }
            }
            else {
                opserr << "ZeroLengthContactASDimplex: WARNING: -orient: model dimension is invalid! \n";
                return theElement;
            }
        }
        else if (strcmp(inputstring, "-intType") == 0) {                             // #2 read type of integration 
            numdata = 1;
            if (OPS_GetIntInput(&numdata, &integrationType) < 0) {
                opserr << "ZeroLengthContactASDimplex: WARNING: invalid integer after -intType\n";
                return theElement;
            }
        }
    }
    // input reading stage is complete

    // check integration type and pick implicit if neither 1 or 0
    if (integrationType != 1 && integrationType != 0) {
        opserr << "ZeroLengthContactASDimplex: WARNING: type of integration is set to IMPLICIT due to invalid flag\n";
        integrationType = false;
    }
    // check the normal vector and normalize
    if (x_e.Norm() < SmallNumber) {
        opserr << "ZeroLengthContactASDimplex: WARNING: normal vector is NOT valid!: -orient $x1 $x2 $x3 cannot be < 0, 0, 0 >\n";
        return theElement;
    }
    x_e.Normalize(); // normalized it on input!

    // finally, create the element
    theElement = new ZeroLengthContactASDimplex(idata[0], idata[1], idata[2], ddata[0], ddata[1],
        ddata[2], ndm, integrationType, x_e[0], x_e[1], x_e[2]);

    if (theElement == 0) {
        opserr << "WARNING: out of memory: element zeroLengthContactASDimplex " << idata[0] <<
            " iNode? jNode? Kn? Kt? mu? <-orient $x1 $x2 $x3> <-intType type?>\n";
    }

    return theElement;
}

ZeroLengthContactASDimplex::ZeroLengthContactASDimplex(int tag, int Nd1, int Nd2,
    double Kn, double Kt, double fcoeff,
    int ndm, bool itype, double xN, double yN, double zN)
    : Element(tag, ELE_TAG_ZeroLengthContactASDimplex)
    , connectedExternalNodes(2)
    , Knormal(Kn)
    , Kfriction(Kt)
    , mu(fcoeff)
    , numDIM(ndm)
    , use_implex(itype)
{
    connectedExternalNodes(0) = Nd1;
    connectedExternalNodes(1) = Nd2;
    Xorient[0] = xN;
    Xorient[1] = yN;
    Xorient[2] = zN;
}

ZeroLengthContactASDimplex::ZeroLengthContactASDimplex()
    : Element(0, ELE_TAG_ZeroLengthContactASDimplex)
{
}

ZeroLengthContactASDimplex::~ZeroLengthContactASDimplex()
{
}

int ZeroLengthContactASDimplex::getNumExternalNodes() const
{
    return 2;
}

const ID& ZeroLengthContactASDimplex::getExternalNodes()
{
    return connectedExternalNodes;
}

Node** ZeroLengthContactASDimplex::getNodePtrs()
{
    return theNodes.data();
}

int ZeroLengthContactASDimplex::getNumDOF()
{
    return numDOF[0] + numDOF[1];
}

void ZeroLengthContactASDimplex::setDomain(Domain* theDomain)
{
    // check Domain is not null - invoked when object removed from a domain
    if (theDomain == 0) {
        theNodes[0] = nullptr;
        theNodes[1] = nullptr;
        return;
    }

    // first set the node pointers
    int Nd1 = connectedExternalNodes(0);
    int Nd2 = connectedExternalNodes(1);
    theNodes[0] = theDomain->getNode(Nd1);
    theNodes[1] = theDomain->getNode(Nd2);

    // check nodes
    if ((theNodes[0] == nullptr) || (theNodes[1] == nullptr)) {
        opserr <<
            "FATAL ERROR ZeroLengthContactASDimplex::setDomain() - Nd1: " << Nd1 <<
            " and/or Nd2: " << Nd2 << " do not exist in the model.\n";
        exit(-1);
    }

    // check NDM
    if (theNodes[0]->getCrds().Size() != numDIM || theNodes[1]->getCrds().Size() != numDIM) {
        opserr <<
            "FATAL ERROR ZeroLengthContactASDimplex::setDomain() - Nd1: " << Nd1 <<
            " and/or Nd2: " << Nd2 << " have an incorrect number of coordinates.\n"
            "Element NDM = " << numDIM << "\n"
            "NDM at Nd1: " << theNodes[0]->getCrds().Size() << "\n"
            "NDM at Nd2: " << theNodes[1]->getCrds().Size() << "\n";
        exit(-1);
    }

    // now determine the number of dof and the dimension
    numDOF[0] = theNodes[0]->getNumberDOF();
    numDOF[1] = theNodes[1]->getNumberDOF();

    // check dofs
    if (numDIM == 2) {
        for (int i = 0; i < 2; ++i) {
            int idof = numDOF[i];
            if (idof < 2 || idof > 3) {
                opserr <<
                    "FATAL ERROR ZeroLengthContactASDimplex::setDomain() - #DOFs ("
                    << idof << ") at node " << i + 1
                    << " is not supported! it can be either 2 or 3\n";
                exit(-1);
            }
        }
    }
    else {
        for (int i = 0; i < 2; ++i) {
            int idof = numDOF[i];
            if ((idof != 3) && (idof != 4) && (idof != 6)) {
                opserr <<
                    "FATAL ERROR ZeroLengthContactASDimplex::setDomain() - #DOFs ("
                    << idof << ") at node " << i + 1
                    << " is not supported! it can be either 3, 4 or 6\n";
                exit(-1);
            }
        }
    }

    // compute the initial gap vector in global coordinates
    // accounting for the geometrical gap and initial displacement
    if (!gap0_initialized) {
        const Vector& P0 = theNodes[0]->getCrds();
        const Vector& P1 = theNodes[1]->getCrds();
        const Vector& U0 = theNodes[0]->getTrialDisp();
        const Vector& U1 = theNodes[1]->getTrialDisp();
        gap0.Zero();
        for (int i = 0; i < numDIM; ++i)
            gap0(i) = P1(i) - U1(i) - P0(i) + U0(i);
        gap0_initialized = true;
    }

    // call the base class method
    DomainComponent::setDomain(theDomain);
}

int ZeroLengthContactASDimplex::commitState(void)
{
    // do the implicit correction if impl-ex
    if (use_implex) {
        // update material internal variables
        updateInternal(false, false);  // explicit_phase?, do_tangent?
    }

    // commit internal variables
    sv.eps_commit = sv.eps;
    sv.shear_commit = sv.shear;
    sv.xs_commit = sv.xs;
    sv.rs_commit_old = sv.rs_commit;
    sv.rs_commit = sv.rs;
    sv.cres_commit_old = sv.cres_commit;
    sv.cres_commit = sv.cres;
    sv.PC_commit = sv.PC;
    sv.dtime_n_commit = sv.dtime_n;

    // done
    return 0;
}

int ZeroLengthContactASDimplex::revertToLastCommit(void)
{
    // restore committed internal variables
    sv.eps = sv.eps_commit;
    sv.shear = sv.shear_commit;
    sv.xs = sv.xs_commit;
    sv.rs = sv.rs_commit;
    sv.cres = sv.cres_commit;
    sv.PC = sv.PC_commit;
    sv.dtime_n = sv.dtime_n_commit;

    //done
    return 0;
}

int ZeroLengthContactASDimplex::revertToStart()
{
    // reset state variables
    sv = StateVariables();
    // done
    return 0;
}

int ZeroLengthContactASDimplex::update()
{
    if (!sv.dtime_is_user_defined) {
        sv.dtime_n = ops_Dt;
        if (!sv.dtime_first_set) {
            sv.dtime_n_commit = sv.dtime_n;
            sv.dtime_first_set = true;
        }
    }
    computeStrain();
    if (use_implex) {
        updateInternal(true, true);
        sv.sig_implex = sv.sig;
    }
    else {
        // for the implicit case we use the numerical tangent... seems more stable
        static Vector strain(3);
        static Matrix Cnum(3, 3);
        constexpr double pert = 1.0e-9;
        strain = sv.eps;
        for (int j = 0; j < 3; ++j) {
            sv.eps(j) = strain(j) + pert;
            updateInternal(true, false);
            for (int i = 0; i < 3; ++i)
                Cnum(i, j) = sv.sig(i);
            sv.eps(j) = strain(j) - pert;
            updateInternal(true, false);
            for (int i = 0; i < 3; ++i)
                Cnum(i, j) = (Cnum(i, j) - sv.sig(i)) / 2.0 / pert;
            sv.eps(j) = strain(j);
        }
        updateInternal(true, false);
        sv.C = Cnum;
    }

    return 0;
}

const Matrix& ZeroLengthContactASDimplex::getTangentStiff()
{
    auto& gs = getGlobalStorage(numDOF[0] + numDOF[1]);
    auto& stiff = gs.K;
    const auto& C = sv.C;
    formStiffnessMatrix(C, stiff);
    return stiff;
}

const Matrix& ZeroLengthContactASDimplex::getInitialStiff()
{
    auto& gs = getGlobalStorage(numDOF[0] + numDOF[1]);
    auto& stiff = gs.K0;
    static Matrix C0(3, 3);
    C0.Zero();
    const Vector& LGap = getInitialGap();
    double Un = LGap(0); // gap parallel to normal vector
    if (Un <= 1.0e-10) { // send elastic stiffness, if contact
        C0(0, 0) = Knormal;
        C0(1, 1) = C0(2, 2) = Kfriction;
    }
    formStiffnessMatrix(C0, stiff);
    return stiff;
}

const Matrix& ZeroLengthContactASDimplex::getDamp()
{
    // get global storage for gloabl DOFset
    auto& gs = getGlobalStorage(numDOF[0] + numDOF[1]);
    gs.D.Zero();
    return gs.D;
}

const Matrix& ZeroLengthContactASDimplex::getMass()
{
    // get global storage for global DOFset
    auto& gs = getGlobalStorage(numDOF[0] + numDOF[1]);
    gs.M.Zero();
    return gs.M;
}

const Vector& ZeroLengthContactASDimplex::getResistingForce() {

    auto& gs = getGlobalStorage(numDOF[0] + numDOF[1]);
    auto& R = gs.R;

    // stress vector in local coordinate system & local dof-set
    const auto& F = sv.sig;

    // residual vector in local coordinate system & global dof-set
    static Vector RL(6);
    const Matrix& B = theBMatrix();
    RL.addMatrixTransposeVector(0.0, B, F, 1.0);

    // residual vector in global system, local DOFset
    static Vector RG(6);
    const Matrix& T2 = getRotationMatrix66();
    RG.addMatrixTransposeVector(0.0, T2, RL, 1.0);

    // compute global residual in global DOFset
    R.Zero();
    int index = numDOF[0];
    for (int i = 0; i < numDIM; i++) {
        R(i) = RG(i);
        R(i + index) = RG(i + 3);
    }
    return R;
}

const Vector& ZeroLengthContactASDimplex::getResistingForceIncInertia()
{
    return getResistingForce();
}

int ZeroLengthContactASDimplex::sendSelf(int commitTag, Channel& theChannel) {

    int res = 0;
    int dataTag = this->getDbTag();

    // int data
    static ID idata(10);
    idata(0) = getTag();
    idata(1) = numDIM;
    idata(2) = numDOF[0];
    idata(3) = numDOF[1];
    idata(4) = connectedExternalNodes[0];
    idata(5) = connectedExternalNodes[1];
    idata(6) = use_implex ? 1 : 0;
    idata(7) = sv.dtime_is_user_defined ? 1 : 0;
    idata(8) = sv.dtime_first_set ? 1 : 0;
    idata(9) = gap0_initialized ? 1 : 0;
    res = theChannel.sendID(dataTag, commitTag, idata);
    if (res < 0) {
        opserr << "WARNING ZeroLengthContactASDimplex::sendSelf() - " << this->getTag() << " failed to send ID\n";
        return -1;
    }

    // double data
    static Vector ddata(31);
    ddata(0) = Knormal;
    ddata(1) = Kfriction;
    ddata(2) = mu;
    ddata(3) = Xorient(0);
    ddata(4) = Xorient(1);
    ddata(5) = Xorient(2);
    ddata(6) = sv.eps(0);
    ddata(7) = sv.eps(1);
    ddata(8) = sv.eps(2);
    ddata(9) = sv.eps_commit(0);
    ddata(10) = sv.eps_commit(1);
    ddata(11) = sv.eps_commit(2);
    ddata(12) = sv.shear(0);
    ddata(13) = sv.shear(1);
    ddata(14) = sv.shear_commit(0);
    ddata(15) = sv.shear_commit(1);
    ddata(16) = sv.xs;
    ddata(17) = sv.xs_commit;
    ddata(18) = sv.rs;
    ddata(19) = sv.rs_commit;
    ddata(20) = sv.rs_commit_old;
    ddata(21) = sv.cres;
    ddata(22) = sv.cres_commit;
    ddata(23) = sv.cres_commit_old;
    ddata(24) = sv.PC;
    ddata(25) = sv.PC_commit;
    ddata(26) = sv.dtime_n;
    ddata(27) = sv.dtime_n_commit;
    ddata(28) = gap0(0);
    ddata(29) = gap0(1);
    ddata(30) = gap0(2);
    res = theChannel.sendVector(dataTag, commitTag, ddata);
    if (res < 0) {
        opserr << "WARNING ZeroLengthContactASDimplex::sendSelf() - " << this->getTag() << " failed to send Vector\n";
        return -1;
    }
    
    // done
    return 0;
}

int ZeroLengthContactASDimplex::recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker) {

    int res;
    int dataTag = this->getDbTag();

    // int data
    static ID idata(10);
    res = theChannel.recvID(dataTag, commitTag, idata);
    if (res < 0) {
        opserr << "WARNING ZeroLengthContactASDimplex::recvSelf() - failed to receive ID\n";
        return -1;
    }
    setTag(idata(0));
    numDIM = idata(1);
    numDOF[0] = idata(2);
    numDOF[1] = idata(3);
    connectedExternalNodes[0] = idata(4);
    connectedExternalNodes[1] = idata(5);
    use_implex = idata(6) == 1;
    sv.dtime_is_user_defined = idata(7) == 1;
    sv.dtime_first_set = idata(8) == 1;
    gap0_initialized = idata(9) == 1;

    // double data
    static Vector ddata(31);
    res = theChannel.recvVector(dataTag, commitTag, ddata);
    if (res < 0) {
        opserr << "WARNING ZeroLengthContactASDimplex::recvSelf() - failed to receive Vector\n";
        return -1;
    }
    Knormal = ddata(0);
    Kfriction = ddata(1);
    mu = ddata(2);
    Xorient(0) = ddata(3);
    Xorient(1) = ddata(4);
    Xorient(2) = ddata(5);
    sv.eps(0) = ddata(6);
    sv.eps(1) = ddata(7);
    sv.eps(2) = ddata(8);
    sv.eps_commit(0) = ddata(9);
    sv.eps_commit(1) = ddata(10);
    sv.eps_commit(2) = ddata(11);
    sv.shear(0) = ddata(12);
    sv.shear(1) = ddata(13);
    sv.shear_commit(0) = ddata(14);
    sv.shear_commit(1) = ddata(15);
    sv.xs = ddata(16);
    sv.xs_commit = ddata(17);
    sv.rs = ddata(18);
    sv.rs_commit = ddata(19);
    sv.rs_commit_old = ddata(20);
    sv.cres = ddata(21);
    sv.cres_commit = ddata(22);
    sv.cres_commit_old = ddata(23);
    sv.PC = ddata(24);
    sv.PC_commit = ddata(25);
    sv.dtime_n = ddata(26);
    sv.dtime_n_commit = ddata(27);
    gap0(0) = ddata(28);
    gap0(1) = ddata(39);
    gap0(2) = ddata(30);

    return 0;
}

int ZeroLengthContactASDimplex::displaySelf(Renderer& theViewer, int displayMode, float fact, const char** modes, int numMode)
{
    if (theNodes[0] == 0 || theNodes[1] == 0)
        return 0;

    static Vector v1(3);
    static Vector v2(3);
    float d1 = 1.0;

    theNodes[0]->getDisplayCrds(v1, 0.);
    theNodes[1]->getDisplayCrds(v2, 0.);
    return theViewer.drawPoint(v1, d1, 10);
}

void ZeroLengthContactASDimplex::Print(OPS_Stream& strm, int flag) {

    if (flag == 0) {
        strm << "Element: " << this->getTag();
        strm << " type: ZeroLengthContactASDimplex  iNode: " << connectedExternalNodes(0);
        strm << " jNode: " << connectedExternalNodes(1) << endln;
    }
    else if (flag == 1) {
        strm << this->getTag() << endln;
    }
}

Response* ZeroLengthContactASDimplex::setResponse(const char** argv, int argc, OPS_Stream& output)
{
    Response* theResponse = nullptr;

    output.tag("ElementOutput");
    output.attr("eleType", "zeroLengthContactASDimplex");
    output.attr("eleTag", this->getTag());
    output.attr("node1", connectedExternalNodes[0]);
    output.attr("node2", connectedExternalNodes[1]);

    auto lam_open_gauss = [&output]() {
        output.tag("GaussPoint");
        output.attr("number", 1);
        output.attr("eta", 0.0);

        output.tag("NdMaterialOutput");
        output.attr("classType", 0);
        output.attr("tag", 0);
    };
    auto lam_close_gauss = [&output]() {
        output.endTag(); // GaussPoint
        output.endTag(); // NdMaterialOutput
    };

    // results which depend on the dimension
    if (numDIM == 2) {
        if (strcmp(argv[0], "force") == 0 || strcmp(argv[0], "forces") == 0) {
            output.tag("ResponseType", "Px_1");
            output.tag("ResponseType", "Py_1");
            output.tag("ResponseType", "Px_2");
            output.tag("ResponseType", "Py_2");
            theResponse = new ElementResponse(this, 1, Vector(4));
        }
        else if (strcmp(argv[0], "displacement") == 0 || strcmp(argv[0], "dispJump") == 0) {
            lam_open_gauss();
            output.tag("ResponseType", "dUx");
            output.tag("ResponseType", "dUy");
            lam_close_gauss();
            theResponse = new ElementResponse(this, 2, Vector(2));
        }
        else if (strcmp(argv[0], "localForce") == 0 || strcmp(argv[0], "localForces") == 0) {
            lam_open_gauss();
            output.tag("ResponseType", "N");
            output.tag("ResponseType", "Tx");
            lam_close_gauss();
            theResponse = new ElementResponse(this, 3, Vector(2));
        }
        else if (strcmp(argv[0], "localForceImplex") == 0 || strcmp(argv[0], "localForcesImplex") == 0) {
            lam_open_gauss();
            output.tag("ResponseType", "N");
            output.tag("ResponseType", "Tx");
            lam_close_gauss();
            theResponse = new ElementResponse(this, 33, Vector(3));
        }
        else if (strcmp(argv[0], "localDisplacement") == 0 || strcmp(argv[0], "localDispJump") == 0) {
            lam_open_gauss();
            output.tag("ResponseType", "dUN");
            output.tag("ResponseType", "dUTx");
            lam_close_gauss();
            theResponse = new ElementResponse(this, 4, Vector(2));
        }
    }
    else {
        if (strcmp(argv[0], "force") == 0 || strcmp(argv[0], "forces") == 0) {
            output.tag("ResponseType", "Px_1");
            output.tag("ResponseType", "Py_1");
            output.tag("ResponseType", "Pz_1");
            output.tag("ResponseType", "Px_2");
            output.tag("ResponseType", "Py_2");
            output.tag("ResponseType", "Pz_2");

            theResponse = new ElementResponse(this, 1, Vector(6));
        }
        else if (strcmp(argv[0], "displacement") == 0 || strcmp(argv[0], "dispJump") == 0) {
            lam_open_gauss();
            output.tag("ResponseType", "dUx");
            output.tag("ResponseType", "dUy");
            output.tag("ResponseType", "dUz");
            lam_close_gauss();
            theResponse = new ElementResponse(this, 2, Vector(3));
        }
        else if (strcmp(argv[0], "localForce") == 0 || strcmp(argv[0], "localForces") == 0) {
            lam_open_gauss();
            output.tag("ResponseType", "N");
            output.tag("ResponseType", "Tx");
            output.tag("ResponseType", "Ty");
            lam_close_gauss();
            theResponse = new ElementResponse(this, 3, Vector(3));
        }
        else if (strcmp(argv[0], "localForceImplex") == 0 || strcmp(argv[0], "localForcesImplex") == 0) {
            lam_open_gauss();
            output.tag("ResponseType", "N");
            output.tag("ResponseType", "Tx");
            output.tag("ResponseType", "Ty");
            lam_close_gauss();
            theResponse = new ElementResponse(this, 33, Vector(3));
        }
        else if (strcmp(argv[0], "localDisplacement") == 0 || strcmp(argv[0], "localDispJump") == 0) {
            lam_open_gauss();
            output.tag("ResponseType", "dUN");
            output.tag("ResponseType", "dUTx");
            output.tag("ResponseType", "dUTy");
            lam_close_gauss();
            theResponse = new ElementResponse(this, 4, Vector(3));
        }
    }

    // other results
    if (theResponse == nullptr) {
        if ((strcmp(argv[0], "slip") == 0) || (strcmp(argv[0], "slipMultiplier") == 0)) {
            lam_open_gauss();
            output.tag("ResponseType", "lambda");
            lam_close_gauss();
            theResponse = new ElementResponse(this, 5, Vector(1));
        }
        else if ((strcmp(argv[0], "NormalContactForce") == 0) || (strcmp(argv[0], "normalContactForce") == 0)) {
            lam_open_gauss();
            output.tag("ResponseType", "N");
            lam_close_gauss();
            theResponse = new ElementResponse(this, 6, Vector(1));
        }
        else if ((strcmp(argv[0], "TangentialContactForce") == 0) || (strcmp(argv[0], "tangentialContactForce") == 0)) {
            lam_open_gauss();
            output.tag("ResponseType", "|T|");
            lam_close_gauss();
            theResponse = new ElementResponse(this, 7, Vector(1));
        }
        else if (strcmp(argv[0], "cres") == 0) {
            lam_open_gauss();
            output.tag("ResponseType", "cres(n+1)");
            output.tag("ResponseType", "cres(n)");
            output.tag("ResponseType", "cres(n-1)");
            lam_close_gauss();
            theResponse = new ElementResponse(this, 8, Vector(3));
        }
    }

    output.endTag(); // ElementOutput
    return theResponse;
}

int ZeroLengthContactASDimplex::getResponse(int responseID, Information& eleInfo) {

    auto& gs = getGlobalStorage(numDOF[0] + numDOF[1]);
    static Vector small(numDIM);
    static Vector large(2 * numDIM);
    static Vector scalar(1);

    if (responseID == 1) {
        // global contact forces
        const Vector& nodeForces = this->getResistingForce();
        for (int i = 0; i < numDIM; i++)
        {
            large(i) = nodeForces(i);
            large(i + numDIM) = nodeForces(i + numDOF[0]);
        }
        return eleInfo.setVector(large);
    }
    else if (responseID == 2) {
        // global displacement jump
        const Matrix& T1 = getRotationMatrix33();
        static Vector dU(3);
        dU.addMatrixTransposeVector(0.0, T1, sv.eps, 1.0);
        for (int i = 0; i < numDIM; i++) {
            small(i) = dU(i);
        }
        return eleInfo.setVector(small);
    }
    else if (responseID == 3) {
        // local contact forces
        for (int i = 0; i < numDIM; i++) {
            small(i) = sv.sig(i);
        }
        return eleInfo.setVector(small);
    }
    else if (responseID == 33) {
        // local contact forces
        for (int i = 0; i < numDIM; i++) {
            small(i) = sv.sig_implex(i);
        }
        return eleInfo.setVector(small);
    }
    else if (responseID == 4) {
        // local displacement jump
        for (int i = 0; i < numDIM; i++) {
            small(i) = sv.eps(i);
        }
        return eleInfo.setVector(small);
    }
    else if (responseID == 5) {
        // material slip multiplier
        scalar(0) = sv.xs;
        return eleInfo.setVector(scalar);
    }
    else if (responseID == 6) {
        // normal contanct force
        scalar(0) = sv.sig(0);
        return eleInfo.setVector(scalar);
    }
    else if (responseID == 7) {
        // tangential contanct force norm
        scalar(0) = sqrt(sv.sig(1)*sv.sig(1) + sv.sig(2)*sv.sig(2));
        return eleInfo.setVector(scalar);
    }
    else if (responseID == 8) {
        // cres internal variable hist
        static Vector cres(3);
        cres(0) = sv.cres; cres(1) = sv.cres_commit; cres(2) = sv.cres_commit_old;
        return eleInfo.setVector(cres);
    }
    else {
        return -1;
    }
}

int ZeroLengthContactASDimplex::updateParameter(int parameterID, double value)
{
    if (parameterID == 1) {
        // set user defined current time increment
        // this is useful for rate dependency in implicit mode and for the implex model
        // when using arc length or displacement control methods, where the pseudo time step
        // is actually the load factor.
        // if when this variable is first set or when it is set before the first commit
        // we set the committed variable to the same value
        sv.dtime_n = value;
        if (!sv.dtime_first_set) {
            sv.dtime_n_commit = sv.dtime_n;
            sv.dtime_first_set = true;
        }
        sv.dtime_is_user_defined = true;
    }
    // done
    return 0;
}

const Matrix& ZeroLengthContactASDimplex::getRotationMatrix33()
{
    static Matrix T(3, 3);

    // initialize orthogonal vectors to Xorient
    static Vector rY(3);
    static Vector rZ(3);

    // create global Y and Z axes
    // static and initialized only once
    auto make3DVector = [](double x, double y, double z) {
        Vector v(3);
        v(0) = x;
        v(1) = y;
        v(2) = z;
        return v;
    };
    static Vector gY = make3DVector(0.0, 1.0, 0.0);   // global Y axis
    static Vector gZ = make3DVector(0.0, 0.0, 1.0);   // global Z axis

    // compute two orthonal vectors to Xorient
    if (fabs(Xorient ^ gY) < 0.99) { // assume X is not parallel to global Y
        cross(Xorient, gY, rZ);
        rZ.Normalize();
        cross(rZ, Xorient, rY);
        rY.Normalize();
    }
    else { // if not, assume X is not parallel to global Z
        cross(Xorient, gZ, rY);
        rY.Normalize();
        cross(rY, Xorient, rZ);
        rZ.Normalize();
    }

    // fill the local 3x3 orientation matrix (transpoed of rotation matrix)
    for (int j = 0; j < 3; j++) { 
        T(0, j) = Xorient(j);
        T(1, j) = rY(j);
        T(2, j) = rZ(j);
    }

    return T;
}

const Matrix& ZeroLengthContactASDimplex::getRotationMatrix66()
{
    static Matrix T2(6, 6);
    T2.Zero();

    const Matrix& T1 = getRotationMatrix33();

    for (int q = 0; q < 2; ++q) {
        int index = q * 3;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                T2(index + i, index + j) = T1(i, j);
    }

    return T2;
}

const Vector& ZeroLengthContactASDimplex::getInitialGap()
{
    // compute local gap vector
    static Vector LGap(3);
    const Matrix& T1 = getRotationMatrix33();
    LGap.addMatrixVector(0.0, T1, gap0, 1.0);
    return LGap;
}

void ZeroLengthContactASDimplex::computeStrain()
{
    // get global displacement for global DOFset
    const Vector& U1 = theNodes[0]->getTrialDisp();
    const Vector& U2 = theNodes[1]->getTrialDisp();

    // set global displacement for local DOFset 
    static Vector UGL(6);
    for (int i = 0; i < numDIM; i++) {
        UGL(i) = U1(i);
        UGL(i + 3) = U2(i);
    }

    // compute local displacement for local DOFset
    static Vector UL(6);
    const Matrix& T2 = getRotationMatrix66();
    UL.addMatrixVector(0.0, T2, UGL, 1.0);

    // compute displacement jump
    const Matrix& B = theBMatrix();
    sv.eps.addMatrixVector(0.0, B, UL, 1.0);

    // add the initial gap
    const Vector& LGap = getInitialGap();
    sv.eps.addVector(1.0, LGap, 1.0);
}

void ZeroLengthContactASDimplex::updateInternal(bool do_implex, bool do_tangent)
{
    // strain layout
    // [Normal, Tangential1, Tangential2]

    // get committed values
    sv.rs = sv.rs_commit;
    sv.xs = sv.xs_commit;
    sv.shear = sv.shear_commit;
    sv.cres = sv.cres_commit;

    // time factor for explicit extrapolation
    double time_factor = 1.0;
    if (do_implex && use_implex && (sv.dtime_n_commit > 0.0))
        time_factor = sv.dtime_n / sv.dtime_n_commit;
    // note: the implex method just wants the ratio of the new to the old time step
    // not the real time step, so it is just fine to assume it to 1.
    // otherwise we have to deal with the problem of the opensees pseudo-time step
    // being the load multiplier in continuation methods...
    time_factor = 1.0;

    // elastic trial
    double SN = Knormal * sv.eps(0);
    double T1 = sv.shear(0) + Kfriction * (sv.eps(1) - sv.eps_commit(1));
    double T2 = sv.shear(1) + Kfriction * (sv.eps(2) - sv.eps_commit(2));

    // tangential stress norm
    double SS = sqrt(T1 * T1 + T2 * T2);

    // compute residual stress in shear
    if (do_implex && use_implex) {
        // explicit extrapolation
        sv.cres = std::max(0.0, sv.cres_commit + time_factor * (sv.cres_commit - sv.cres_commit_old));
    }
    else {
        // implicit
        if (SN < 0.0)
            sv.cres = -mu * SN;
        else {
            if(!use_implex && sv.eps(0) < 1.0e-6)
                sv.cres = 1.0e-10;
        }
    }

    // shear stress moved to total stress space
    double SS_bar = SS + sv.xs * Kfriction;

    // compute equivalent shear stress (exceeding cres)
    if (do_implex && use_implex) {
        // explicit extrapolation
        sv.rs = sv.rs_commit + time_factor * (sv.rs_commit - sv.rs_commit_old);
    }
    else {
        // implicit
        double rs_trial = SS_bar - sv.cres;
        sv.rs = std::max(sv.rs, rs_trial);
    }

    // update shear plastic multipler and apparent plastic damage
    double equivalent_total_strain = sv.rs / Kfriction;
    double rs_effective = (equivalent_total_strain - sv.xs) * Kfriction + sv.cres;
    sv.xs = equivalent_total_strain;
    double damage = 0.0;
    if (sv.xs > std::numeric_limits<double>::epsilon()) {
        damage = 1.0;
        if (rs_effective > std::numeric_limits<double>::epsilon())
            damage = 1.0 - sv.cres / rs_effective;
    }

    // extract compressive normal stress
    if (do_implex && use_implex) {
        // explicit extrapolation (actually keep the commited one... since it's either 0 or 1)
        sv.PC = sv.PC_commit;
    }
    else {
        sv.PC = SN <= 0.0 ? 1.0 : 0.0;
    }
    double SC = sv.PC * SN;

    // update effective shear stress
    sv.shear(0) = (1.0 - damage) * T1;
    sv.shear(1) = (1.0 - damage) * T2;

    // compute nominal stress
    sv.sig(0) = SC;
    sv.sig(1) = sv.shear(0);
    sv.sig(2) = sv.shear(1);

    // compute tangent tensor
    if (do_tangent) {
        sv.C.Zero();
        // the normal part is the same for both implicit and implex
        sv.C(0, 0) = sv.PC * Knormal;
        // the secant part is also the same
        sv.C(1, 1) = sv.C(2, 2) = (1.0 - damage) * Kfriction;
        if (!use_implex) {
            // implicit algorithmic tangent
            // may be non symmetric due to the derivative of damage
            // and the dependency of cres on normal stress
            if (sv.cres > std::numeric_limits<double>::epsilon()) {
                if (sv.rs > sv.rs_commit) {
                    // plastic loading
                    double denom_n = SN * mu + SS;
                    double denom_t = SS - sv.cres;
                    double UTsqnorm = sv.eps(1) * sv.eps(1) + sv.eps(2) * sv.eps(2);
                    double dE1 = 0.0;
                    double dE2 = 0.0;
                    if (UTsqnorm > std::numeric_limits<double>::epsilon()) {
                        double dE1 = sv.eps(1) / UTsqnorm;
                        double dE2 = sv.eps(2) / UTsqnorm;
                    }
                    double dDdEn = Knormal * SS * mu / (denom_n * denom_n);
                    double dDdE1 = SS * sv.cres * dE1 / (denom_t * denom_t);
                    double dDdE2 = SS * sv.cres * dE2 / (denom_t * denom_t);
                    sv.C(1, 0) = -T1 * dDdEn;
                    sv.C(2, 0) = -T2 * dDdEn;
                    sv.C(1, 1) = -T1 * dDdE1;
                    sv.C(2, 1) = -T2 * dDdE1;
                    sv.C(1, 2) = -T1 * dDdE2;
                    sv.C(2, 2) = -T2 * dDdE2;
                }
            }
        }
    }
}

void ZeroLengthContactASDimplex::formStiffnessMatrix(const Matrix& C, Matrix& K)
{
    // element stiffness in local system
    static Matrix KL(6, 6);
    const Matrix& B = theBMatrix();
    KL.addMatrixTripleProduct(0.0, B, C, 1.0);

    // element stiffness in global system, local DOFset
    static Matrix KG(6, 6);
    const Matrix& T2 = getRotationMatrix66();
    KG.addMatrixTripleProduct(0.0, T2, KL, 1.0);

    // element stiffness in global DOFset
    K.Zero();
    int index = numDOF[0];
    for (int i = 0; i < numDIM; i++) {
        for (int j = 0; j < numDIM; j++) {
            K(i, j) = KG(i, j);
            K(i + index, j) = KG(i + 3, j);
            K(i, j + index) = KG(i, j + 3);
            K(i + index, j + index) = KG(i + 3, j + 3);
        }
    }
}

const Matrix& ZeroLengthContactASDimplex::theBMatrix()
{
    // this is the strain-displacement matrix.
    // the strain is a actually the displacement jump U2-U1
    static Matrix B(3, 6);
    B.Zero();
    for (int i = 0; i < 3; i++) {
        B(i, i + 3) = 1;
        B(i, i) = -1;
    }
    return B;
}
