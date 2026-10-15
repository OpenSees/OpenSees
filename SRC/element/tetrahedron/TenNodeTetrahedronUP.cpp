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

// ============================================================================
// 2026 By Marlio Molina Manrique @ Politecnico di Torino, Turin, Italy
// and Pedro Arduino @ University of Washington, WA.
// marlio.molina@polito.it | parduino@uw.edu
// ============================================================================
// Please read detailed description in TenNodeTetrahedronUP.h.
// ============================================================================

// Ten-node displacement / four-node pressure tetrahedron.

#include <TenNodeTetrahedronUP.h>

#include <Channel.h>
#include <Domain.h>
#include <ElementResponse.h>
#include <ElementalLoad.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <Node.h>
#include <OPS_Globals.h>
#include <OPS_Stream.h>
#include <Parameter.h>
#include <Renderer.h>
#include <elementAPI.h>

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>

namespace {
const double tetAlpha = (5.0 + 3.0 * std::sqrt(5.0)) / 20.0;
const double tetBeta = (5.0 - std::sqrt(5.0)) / 20.0;
const double tetWeight = 1.0 / 24.0;
const double solidGaussPoints[4][3] = {
    {tetAlpha, tetBeta, tetBeta},
    {tetBeta, tetAlpha, tetBeta},
    {tetBeta, tetBeta, tetAlpha},
    {tetBeta, tetBeta, tetBeta}
};

// Keast 11-point tetrahedral rule, exact through total polynomial degree 4.
// The negative centroid weight is part of the rule. Weights sum to 1/6.
const double massGaussPoints[11][3] = {
    {0.25, 0.25, 0.25},
    {11.0 / 14.0, 1.0 / 14.0, 1.0 / 14.0},
    {1.0 / 14.0, 11.0 / 14.0, 1.0 / 14.0},
    {1.0 / 14.0, 1.0 / 14.0, 11.0 / 14.0},
    {1.0 / 14.0, 1.0 / 14.0, 1.0 / 14.0},
    {0.3994035761667992, 0.3994035761667992, 0.1005964238332008},
    {0.3994035761667992, 0.1005964238332008, 0.3994035761667992},
    {0.3994035761667992, 0.1005964238332008, 0.1005964238332008},
    {0.1005964238332008, 0.3994035761667992, 0.3994035761667992},
    {0.1005964238332008, 0.3994035761667992, 0.1005964238332008},
    {0.1005964238332008, 0.1005964238332008, 0.3994035761667992}
};
const double massGaussWeights[11] = {
    -74.0 / 5625.0,
    343.0 / 45000.0, 343.0 / 45000.0,
    343.0 / 45000.0, 343.0 / 45000.0,
    56.0 / 2250.0, 56.0 / 2250.0, 56.0 / 2250.0,
    56.0 / 2250.0, 56.0 / 2250.0, 56.0 / 2250.0
};
}

void *
OPS_TenNodeTetrahedronUP()
{
    static int num_TenNodeTetrahedronUP = 0;
    if (num_TenNodeTetrahedronUP == 0) {
        num_TenNodeTetrahedronUP++;
        opserr << "TenNodeTetrahedronUP element - Written: M.Molina_Manrique, and P.Arduino, U.Washington\n";
    }

    if (OPS_GetNDM() != 3) {
        opserr << "WARNING TenNodeTetrahedronUP requires ndm = 3" << endln;
        return nullptr;
    }

    // eleTag, ten node tags, matTag, bulk, fluid density, kx, ky, kz
    if (OPS_GetNumRemainingInputArgs() < 17) {
        opserr << "WARNING insufficient arguments for TenNodeTetrahedronUP"
               << endln;
        opserr << "Want: element TenNodeTetrahedronUP eleTag? "
               << "n1? n2? n3? n4? n5? n6? n7? n8? n9? n10? "
               << "matTag? bulk? fluidDensity? permX? permY? permZ? "
               << "<b1? b2? b3?>" << endln;
        return nullptr;
    }

    int integerData[12];
    int numberData = 12;
    if (OPS_GetIntInput(&numberData, integerData) < 0) {
        opserr << "WARNING TenNodeTetrahedronUP received invalid integer data"
               << endln;
        return nullptr;
    }

    NDMaterial *material = OPS_getNDMaterial(integerData[11]);
    if (material == nullptr) {
        opserr << "WARNING TenNodeTetrahedronUP material " << integerData[11]
               << " was not found for element " << integerData[0] << endln;
        return nullptr;
    }

    double physicalData[5];
    numberData = 5;
    if (OPS_GetDoubleInput(&numberData, physicalData) < 0) {
        opserr << "WARNING TenNodeTetrahedronUP received invalid fluid data"
               << endln;
        return nullptr;
    }

    if (physicalData[0] <= 0.0 || physicalData[1] < 0.0 ||
        physicalData[2] < 0.0 || physicalData[3] < 0.0 ||
        physicalData[4] < 0.0) {
        opserr << "WARNING TenNodeTetrahedronUP requires bulk > 0, "
               << "fluidDensity >= 0, and permeability >= 0" << endln;
        return nullptr;
    }

    double bodyForce[3] = {0.0, 0.0, 0.0};
    numberData = OPS_GetNumRemainingInputArgs();
    if (numberData > 3) {
        numberData = 3;
    }
    if (numberData > 0 &&
        OPS_GetDoubleInput(&numberData, bodyForce) < 0) {
        opserr << "WARNING TenNodeTetrahedronUP received invalid body forces"
               << endln;
        return nullptr;
    }

    return new TenNodeTetrahedronUP(
        integerData[0], integerData[1], integerData[2], integerData[3],
        integerData[4], integerData[5], integerData[6], integerData[7],
        integerData[8], integerData[9], integerData[10], *material,
        physicalData[0], physicalData[1], physicalData[2], physicalData[3],
        physicalData[4], bodyForce[0], bodyForce[1], bodyForce[2]);
}

Matrix TenNodeTetrahedronUP::stiffness(NumDOF, NumDOF);
Matrix TenNodeTetrahedronUP::damping(NumDOF, NumDOF);
Matrix TenNodeTetrahedronUP::mass(NumDOF, NumDOF);
Vector TenNodeTetrahedronUP::residual(NumDOF);
Matrix TenNodeTetrahedronUP::B(NumStressComponents, 3);

TenNodeTetrahedronUP::TenNodeTetrahedronUP()
    : Element(0, ELE_TAG_TenNodeTetrahedronUP),
      connectedExternalNodes(NumNodes), applyLoad(0), fluidDensity(0.0),
      combinedBulkModulus(0.0), load(nullptr), initialStiffness(nullptr)
{
    for (int node = 0; node < NumNodes; ++node) {
        nodePointers[node] = nullptr;
        for (int direction = 0; direction < 3; ++direction) {
            coordinates[direction][node] = 0.0;
        }
    }
    for (int point = 0; point < NumSolidGaussPoints; ++point) {
        materialPointers[point] = nullptr;
    }
    for (int direction = 0; direction < 3; ++direction) {
        bodyForce[direction] = 0.0;
        appliedBodyForce[direction] = 0.0;
        permeability[direction] = 0.0;
    }
}

TenNodeTetrahedronUP::TenNodeTetrahedronUP(
    int tag, int node1, int node2, int node3, int node4, int node5,
    int node6, int node7, int node8, int node9, int node10,
    NDMaterial &theMaterial, double bulk, double rhoFluid, double permX,
    double permY, double permZ, double b1, double b2, double b3)
    : Element(tag, ELE_TAG_TenNodeTetrahedronUP),
      connectedExternalNodes(NumNodes), applyLoad(0),
      fluidDensity(rhoFluid), combinedBulkModulus(bulk), load(nullptr),
      initialStiffness(nullptr)
{
    const int nodeTags[NumNodes] = {node1, node2, node3, node4, node5,
                                    node6, node7, node8, node9, node10};
    for (int node = 0; node < NumNodes; ++node) {
        connectedExternalNodes(node) = nodeTags[node];
        nodePointers[node] = nullptr;
        for (int direction = 0; direction < 3; ++direction) {
            coordinates[direction][node] = 0.0;
        }
    }

    for (int point = 0; point < NumSolidGaussPoints; ++point) {
        materialPointers[point] =
            theMaterial.getCopy("ThreeDimensional");
        if (materialPointers[point] == nullptr) {
            opserr << "TenNodeTetrahedronUP::constructor -- material "
                   << theMaterial.getTag()
                   << " cannot provide a ThreeDimensional copy" << endln;
            std::exit(-1);
        }
    }

    bodyForce[0] = b1;
    bodyForce[1] = b2;
    bodyForce[2] = b3;
    appliedBodyForce[0] = 0.0;
    appliedBodyForce[1] = 0.0;
    appliedBodyForce[2] = 0.0;
    permeability[0] = permX;
    permeability[1] = permY;
    permeability[2] = permZ;
}

TenNodeTetrahedronUP::~TenNodeTetrahedronUP()
{
    for (int point = 0; point < NumSolidGaussPoints; ++point) {
        delete materialPointers[point];
    }
    delete load;
    delete initialStiffness;
}

int
TenNodeTetrahedronUP::displacementDOF(int node, int component)
{
    // Corner nodes store [ux, uy, uz, p].  Midside nodes store [ux,uy,uz].
    return node < NumPressureNodes ? 4 * node + component
                                   : 16 + 3 * (node - NumPressureNodes) + component;
}

int
TenNodeTetrahedronUP::pressureDOF(int cornerNode)
{
    return 4 * cornerNode + 3;
}

void
TenNodeTetrahedronUP::setDomain(Domain *theDomain)
{
    if (theDomain == nullptr) {
        for (int node = 0; node < NumNodes; ++node) {
            nodePointers[node] = nullptr;
        }
        this->DomainComponent::setDomain(nullptr);
        return;
    }

    for (int node = 0; node < NumNodes; ++node) {
        nodePointers[node] = theDomain->getNode(connectedExternalNodes(node));
        if (nodePointers[node] == nullptr) {
            opserr << "TenNodeTetrahedronUP::setDomain -- node "
                   << connectedExternalNodes(node) << " does not exist for element "
                   << this->getTag() << endln;
            return;
        }

        const int requiredDOF = node < NumPressureNodes ? 4 : 3;
        if (nodePointers[node]->getNumberDOF() != requiredDOF) {
            opserr << "TenNodeTetrahedronUP::setDomain -- node "
                   << connectedExternalNodes(node) << " has "
                   << nodePointers[node]->getNumberDOF() << " DOFs; expected "
                   << requiredDOF << " for element " << this->getTag() << endln;
            return;
        }

        const Vector &nodeCoordinates = nodePointers[node]->getCrds();
        if (nodeCoordinates.Size() < 3) {
            opserr << "TenNodeTetrahedronUP::setDomain -- node "
                   << connectedExternalNodes(node)
                   << " does not have three coordinates" << endln;
            return;
        }
        for (int direction = 0; direction < 3; ++direction) {
            coordinates[direction][node] = nodeCoordinates(direction);
        }
    }

    this->DomainComponent::setDomain(theDomain);
}

int TenNodeTetrahedronUP::getNumExternalNodes() const { return NumNodes; }

const ID &TenNodeTetrahedronUP::getExternalNodes()
{
    return connectedExternalNodes;
}

Node **TenNodeTetrahedronUP::getNodePtrs() { return nodePointers; }

int TenNodeTetrahedronUP::getNumDOF() { return NumDOF; }

int
TenNodeTetrahedronUP::commitState()
{
    int result = this->Element::commitState();
    for (int point = 0; point < NumSolidGaussPoints; ++point) {
        result += materialPointers[point]->commitState();
    }
    return result;
}

int
TenNodeTetrahedronUP::revertToLastCommit()
{
    int result = 0;
    for (int point = 0; point < NumSolidGaussPoints; ++point) {
        result += materialPointers[point]->revertToLastCommit();
    }
    return result;
}

int
TenNodeTetrahedronUP::revertToStart()
{
    int result = 0;
    for (int point = 0; point < NumSolidGaussPoints; ++point) {
        result += materialPointers[point]->revertToStart();
    }
    return result;
}

int
TenNodeTetrahedronUP::update()
{
    int result = 0;
    double shape[4][NumNodes];

    for (int point = 0; point < NumSolidGaussPoints; ++point) {
        double determinant = 0.0;
        this->shapeFunction(solidGaussPoints[point], determinant, shape);
        if (std::abs(determinant) <= 1.0e-18) {
            opserr << "TenNodeTetrahedronUP::update -- singular Jacobian "
                   << determinant << " in element " << this->getTag()
                   << ", integration point " << point + 1 << endln;
            return -1;
        }

        Vector strain(NumStressComponents);
        strain.Zero();
        for (int node = 0; node < NumNodes; ++node) {
            const Vector &displacement = nodePointers[node]->getTrialDisp();
            const Matrix &nodeB = this->computeB(node, shape);
            for (int component = 0; component < NumStressComponents;
                 ++component) {
                for (int direction = 0; direction < 3; ++direction) {
                    strain(component) +=
                        nodeB(component, direction) * displacement(direction);
                }
            }
        }
        result += materialPointers[point]->setTrialStrain(strain);
    }
    return result;
}

const Matrix &
TenNodeTetrahedronUP::getTangentStiff()
{
    this->formResidAndTangent(1);
    return stiffness;
}

const Matrix &
TenNodeTetrahedronUP::getInitialStiff()
{
    if (initialStiffness == nullptr) {
        this->formResidAndTangent(1, true);
        initialStiffness = new Matrix(stiffness);
    }
    return *initialStiffness;
}

const Matrix &
TenNodeTetrahedronUP::getDamp()
{
    this->formDampingTerms(1);
    return damping;
}

const Matrix &
TenNodeTetrahedronUP::getMass()
{
    this->formInertiaTerms(1);
    return mass;
}

void
TenNodeTetrahedronUP::zeroLoad()
{
    if (load != nullptr) {
        load->Zero();
    }
    applyLoad = 0;
    appliedBodyForce[0] = 0.0;
    appliedBodyForce[1] = 0.0;
    appliedBodyForce[2] = 0.0;
}

int
TenNodeTetrahedronUP::addLoad(ElementalLoad *elementLoad, double loadFactor)
{
    int loadType = 0;
    const Vector &loadData = elementLoad->getData(loadType, loadFactor);
    if (loadType != LOAD_TAG_SelfWeight) {
        opserr << "TenNodeTetrahedronUP::addLoad -- unsupported load type "
               << loadType << " for element " << this->getTag() << endln;
        return -1;
    }
    if (loadData.Size() < 3) {
        opserr << "TenNodeTetrahedronUP::addLoad -- selfWeight data are "
               << "incomplete for element " << this->getTag() << endln;
        return -1;
    }

    applyLoad = 1;
    for (int direction = 0; direction < 3; ++direction) {
        appliedBodyForce[direction] +=
            loadFactor * loadData(direction) * bodyForce[direction];
    }
    return 0;
}

int
TenNodeTetrahedronUP::addInertiaLoadToUnbalance(const Vector &acceleration)
{
    Vector transformedAcceleration(NumDOF);
    transformedAcceleration.Zero();

    for (int node = 0; node < NumNodes; ++node) {
        const Vector &nodeAcceleration =
            nodePointers[node]->getRV(acceleration);
        if (nodeAcceleration.Size() < 3) {
            opserr << "TenNodeTetrahedronUP::addInertiaLoadToUnbalance -- "
                   << "incompatible acceleration vector at node "
                   << connectedExternalNodes(node) << endln;
            return -1;
        }
        for (int direction = 0; direction < 3; ++direction) {
            transformedAcceleration(displacementDOF(node, direction)) =
                nodeAcceleration(direction);
        }
    }

    this->formInertiaTerms(1);
    if (load == nullptr) {
        load = new Vector(NumDOF);
        load->Zero();
    }
    load->addMatrixVector(1.0, mass, transformedAcceleration, -1.0);
    return 0;
}

const Vector &
TenNodeTetrahedronUP::getResistingForce()
{
    this->formResidAndTangent(0);
    if (load != nullptr) {
        residual -= *load;
    }
    return residual;
}

const Vector &
TenNodeTetrahedronUP::getResistingForceIncInertia()
{
    Vector acceleration(NumDOF);
    Vector velocity(NumDOF);
    acceleration.Zero();
    velocity.Zero();

    for (int node = 0; node < NumNodes; ++node) {
        const Vector &nodeAcceleration = nodePointers[node]->getTrialAccel();
        const Vector &nodeVelocity = nodePointers[node]->getTrialVel();
        const int requiredSize = node < NumPressureNodes ? 4 : 3;
        if (nodeAcceleration.Size() != requiredSize ||
            nodeVelocity.Size() != requiredSize) {
            opserr << "TenNodeTetrahedronUP::getResistingForceIncInertia -- "
                   << "incompatible trial-vector size at node "
                   << connectedExternalNodes(node) << endln;
            return residual;
        }
        for (int direction = 0; direction < 3; ++direction) {
            const int index = displacementDOF(node, direction);
            acceleration(index) = nodeAcceleration(direction);
            velocity(index) = nodeVelocity(direction);
        }
        if (node < NumPressureNodes) {
            const int index = pressureDOF(node);
            acceleration(index) = nodeAcceleration(3);
            velocity(index) = nodeVelocity(3);
        }
    }

    // Assemble damping before the residual.  With stiffness-proportional
    // Rayleigh damping, formDampingTerms() calls getTangentStiff(), which
    // rebuilds and clears residual.
    this->formDampingTerms(1);
    this->getResistingForce();
    this->formInertiaTerms(1);
    residual.addMatrixVector(1.0, mass, acceleration, 1.0);
    residual.addMatrixVector(1.0, damping, velocity, 1.0);
    return residual;
}

int
TenNodeTetrahedronUP::sendSelf(int commitTag, Channel &channel)
{
    const int dataTag = this->getDbTag();

    Vector data(17);
    data(0) = this->getTag();
    data(1) = combinedBulkModulus;
    data(2) = fluidDensity;
    data(3) = permeability[0];
    data(4) = permeability[1];
    data(5) = permeability[2];
    data(6) = bodyForce[0];
    data(7) = bodyForce[1];
    data(8) = bodyForce[2];
    data(9) = appliedBodyForce[0];
    data(10) = appliedBodyForce[1];
    data(11) = appliedBodyForce[2];
    data(12) = applyLoad;
    data(13) = alphaM;
    data(14) = betaK;
    data(15) = betaK0;
    data(16) = betaKc;
    if (channel.sendVector(dataTag, commitTag, data) < 0) {
        opserr << "TenNodeTetrahedronUP::sendSelf -- failed to send data for "
               << "element " << this->getTag() << endln;
        return -1;
    }

    ID identifiers(18);
    for (int node = 0; node < NumNodes; ++node) {
        identifiers(node) = connectedExternalNodes(node);
    }
    for (int point = 0; point < NumSolidGaussPoints; ++point) {
        identifiers(10 + point) = materialPointers[point]->getClassTag();
        int materialDbTag = materialPointers[point]->getDbTag();
        if (materialDbTag == 0) {
            materialDbTag = channel.getDbTag();
            materialPointers[point]->setDbTag(materialDbTag);
        }
        identifiers(14 + point) = materialDbTag;
    }
    if (channel.sendID(dataTag, commitTag, identifiers) < 0) {
        opserr << "TenNodeTetrahedronUP::sendSelf -- failed to send IDs for "
               << "element " << this->getTag() << endln;
        return -2;
    }

    for (int point = 0; point < NumSolidGaussPoints; ++point) {
        if (materialPointers[point]->sendSelf(commitTag, channel) < 0) {
            opserr << "TenNodeTetrahedronUP::sendSelf -- material " << point
                   << " failed in element " << this->getTag() << endln;
            return -3;
        }
    }
    return 0;
}

int
TenNodeTetrahedronUP::recvSelf(int commitTag, Channel &channel,
                               FEM_ObjectBroker &broker)
{
    const int dataTag = this->getDbTag();

    Vector data(17);
    if (channel.recvVector(dataTag, commitTag, data) < 0) {
        opserr << "TenNodeTetrahedronUP::recvSelf -- failed to receive data"
               << endln;
        return -1;
    }

    this->setTag(static_cast<int>(data(0)));
    combinedBulkModulus = data(1);
    fluidDensity = data(2);
    permeability[0] = data(3);
    permeability[1] = data(4);
    permeability[2] = data(5);
    bodyForce[0] = data(6);
    bodyForce[1] = data(7);
    bodyForce[2] = data(8);
    appliedBodyForce[0] = data(9);
    appliedBodyForce[1] = data(10);
    appliedBodyForce[2] = data(11);
    applyLoad = static_cast<int>(data(12));
    alphaM = data(13);
    betaK = data(14);
    betaK0 = data(15);
    betaKc = data(16);

    ID identifiers(18);
    if (channel.recvID(dataTag, commitTag, identifiers) < 0) {
        opserr << "TenNodeTetrahedronUP::recvSelf -- failed to receive IDs "
               << "for element " << this->getTag() << endln;
        return -2;
    }
    for (int node = 0; node < NumNodes; ++node) {
        connectedExternalNodes(node) = identifiers(node);
        nodePointers[node] = nullptr;
    }

    for (int point = 0; point < NumSolidGaussPoints; ++point) {
        const int materialClassTag = identifiers(10 + point);
        const int materialDbTag = identifiers(14 + point);
        if (materialPointers[point] == nullptr ||
            materialPointers[point]->getClassTag() != materialClassTag) {
            delete materialPointers[point];
            materialPointers[point] =
                broker.getNewNDMaterial(materialClassTag);
            if (materialPointers[point] == nullptr) {
                opserr << "TenNodeTetrahedronUP::recvSelf -- broker could not "
                       << "create NDMaterial class " << materialClassTag
                       << " for element " << this->getTag() << endln;
                return -3;
            }
        }
        materialPointers[point]->setDbTag(materialDbTag);
        if (materialPointers[point]->recvSelf(commitTag, channel, broker) < 0) {
            opserr << "TenNodeTetrahedronUP::recvSelf -- material " << point
                   << " failed in element " << this->getTag() << endln;
            return -4;
        }
    }

    delete initialStiffness;
    initialStiffness = nullptr;
    delete load;
    load = nullptr;
    return 0;
}

void
TenNodeTetrahedronUP::Print(OPS_Stream &stream, int)
{
    stream << "TenNodeTetrahedronUP, element id: " << this->getTag() << endln;
    stream << "  nodes: " << connectedExternalNodes << endln;
    stream << "  combined bulk modulus: " << combinedBulkModulus << endln;
    stream << "  fluid density: " << fluidDensity << endln;
    stream << "  permeability: " << permeability[0] << " " << permeability[1]
           << " " << permeability[2] << endln;
    stream << "  status: mixed Tet10-Tet4 u-p formulation active" << endln;
    for (int point = 0; point < NumSolidGaussPoints; ++point) {
        stream << "  integration point " << point + 1
               << " strain: " << materialPointers[point]->getStrain();
        stream << "  integration point " << point + 1
               << " stress: " << materialPointers[point]->getStress();
    }
}

int
TenNodeTetrahedronUP::displaySelf(Renderer &, int, float, const char **, int)
{
    return 0;
}

Response *
TenNodeTetrahedronUP::setResponse(const char **argv, int argc,
                                  OPS_Stream &output)
{
    if (argc < 1) {
        return nullptr;
    }

    output.tag("ElementOutput");
    output.attr("eleType", "TenNodeTetrahedronUP");
    output.attr("eleTag", this->getTag());
    char label[32];
    for (int node = 0; node < NumNodes; ++node) {
        std::snprintf(label, sizeof(label), "node%d", node + 1);
        output.attr(label, connectedExternalNodes(node));
    }

    Response *response = nullptr;
    if (std::strcmp(argv[0], "damping") == 0 ||
        std::strcmp(argv[0], "damp") == 0 ||
        std::strcmp(argv[0], "coupling") == 0) {
        response = new ElementResponse(this, 1, Matrix(NumDOF, NumDOF));
    } else if (std::strcmp(argv[0], "mass") == 0) {
        response = new ElementResponse(this, 2, Matrix(NumDOF, NumDOF));
    } else if (std::strcmp(argv[0], "force") == 0 ||
               std::strcmp(argv[0], "forces") == 0 ||
               std::strcmp(argv[0], "globalForce") == 0 ||
               std::strcmp(argv[0], "globalForces") == 0) {
        for (int node = 0; node < NumNodes; ++node) {
            for (int direction = 1; direction <= 3; ++direction) {
                std::snprintf(label, sizeof(label), "P%d_%d", direction,
                              node + 1);
                output.tag("ResponseType", label);
            }
            if (node < NumPressureNodes) {
                std::snprintf(label, sizeof(label), "Pp_%d", node + 1);
                output.tag("ResponseType", label);
            }
        }
        response = new ElementResponse(this, 3, Vector(NumDOF));
    } else if (std::strcmp(argv[0], "forceIncInertia") == 0 ||
               std::strcmp(argv[0], "forcesIncInertia") == 0) {
        response = new ElementResponse(this, 4, Vector(NumDOF));
    } else if (std::strcmp(argv[0], "stiff") == 0 ||
               std::strcmp(argv[0], "stiffness") == 0) {
        response = new ElementResponse(this, 5, Matrix(NumDOF, NumDOF));
    } else if (std::strcmp(argv[0], "stresses") == 0 ||
               std::strcmp(argv[0], "stress") == 0) {
        const char *stressLabels[NumStressComponents] = {
            "sigma11", "sigma22", "sigma33",
            "sigma12", "sigma23", "sigma13"};
        for (int point = 0; point < NumSolidGaussPoints; ++point) {
            output.tag("GaussPoint");
            output.attr("number", point + 1);
            output.tag("NdMaterialOutput");
            output.attr("classType", materialPointers[point]->getClassTag());
            output.attr("tag", materialPointers[point]->getTag());
            for (int component = 0; component < NumStressComponents;
                 ++component) {
                output.tag("ResponseType", stressLabels[component]);
            }
            output.endTag();
            output.endTag();
        }
        response = new ElementResponse(
            this, 6,
            Vector(NumSolidGaussPoints * NumStressComponents));
    } else if (std::strcmp(argv[0], "strains") == 0 ||
               std::strcmp(argv[0], "strain") == 0) {
        const char *strainLabels[NumStressComponents] = {
            "eps11", "eps22", "eps33", "eps12", "eps23", "eps13"};
        for (int point = 0; point < NumSolidGaussPoints; ++point) {
            output.tag("GaussPoint");
            output.attr("number", point + 1);
            for (int component = 0; component < NumStressComponents;
                 ++component) {
                output.tag("ResponseType", strainLabels[component]);
            }
            output.endTag();
        }
        response = new ElementResponse(
            this, 7,
            Vector(NumSolidGaussPoints * NumStressComponents));
    } else if (std::strcmp(argv[0], "pressure") == 0 ||
               std::strcmp(argv[0], "pressures") == 0 ||
               std::strcmp(argv[0], "porePressure") == 0 ||
               std::strcmp(argv[0], "porePressures") == 0) {
        for (int node = 0; node < NumPressureNodes; ++node) {
            std::snprintf(label, sizeof(label), "p_%d", node + 1);
            output.tag("ResponseType", label);
        }
        response = new ElementResponse(this, 8, Vector(NumPressureNodes));
    } else if ((std::strcmp(argv[0], "material") == 0 ||
                std::strcmp(argv[0], "integrPoint") == 0 ||
                std::strcmp(argv[0], "integrationPoint") == 0) &&
               argc >= 3) {
        const int point = std::atoi(argv[1]);
        if (point >= 1 && point <= NumSolidGaussPoints) {
            output.tag("GaussPoint");
            output.attr("number", point);
            response = materialPointers[point - 1]->setResponse(
                &argv[2], argc - 2, output);
            output.endTag();
        }
    }

    output.endTag();
    return response;
}

int
TenNodeTetrahedronUP::getResponse(int responseID, Information &information)
{
    if (responseID == 1) {
        return information.setMatrix(this->getDamp());
    }
    if (responseID == 2) {
        return information.setMatrix(this->getMass());
    }
    if (responseID == 3) {
        return information.setVector(this->getResistingForce());
    }
    if (responseID == 4) {
        return information.setVector(this->getResistingForceIncInertia());
    }
    if (responseID == 5) {
        return information.setMatrix(this->getTangentStiff());
    }
    if (responseID == 6 || responseID == 7) {
        static Vector responseValues(NumSolidGaussPoints
                                     * NumStressComponents);
        int index = 0;
        for (int point = 0; point < NumSolidGaussPoints; ++point) {
            const Vector &values = responseID == 6
                                       ? materialPointers[point]->getStress()
                                       : materialPointers[point]->getStrain();
            for (int component = 0; component < NumStressComponents;
                 ++component) {
                responseValues(index++) = values(component);
            }
        }
        return information.setVector(responseValues);
    }
    if (responseID == 8) {
        static Vector pressures(NumPressureNodes);
        for (int node = 0; node < NumPressureNodes; ++node) {
            // OpenSees u-p elements store pore pressure in the fourth entry
            // of the nodal trial-velocity vector.
            pressures(node) = nodePointers[node]->getTrialVel()(3);
        }
        return information.setVector(pressures);
    }
    return -1;
}

void
TenNodeTetrahedronUP::computeBasis()
{
    for (int node = 0; node < NumNodes; ++node) {
        const Vector &nodeCoordinates = nodePointers[node]->getCrds();
        for (int direction = 0; direction < 3; ++direction) {
            coordinates[direction][node] = nodeCoordinates(direction);
        }
    }
}

void
TenNodeTetrahedronUP::shapeFunction(const double naturalCoordinates[3],
                                    double &determinant,
                                    double shape[4][NumNodes]) const
{
    const double l1 = naturalCoordinates[0];
    const double l2 = naturalCoordinates[1];
    const double l3 = naturalCoordinates[2];
    const double l4 = 1.0 - l1 - l2 - l3;

    shape[3][0] = l1 * (2.0 * l1 - 1.0);
    shape[3][1] = l2 * (2.0 * l2 - 1.0);
    shape[3][2] = l3 * (2.0 * l3 - 1.0);
    shape[3][3] = l4 * (2.0 * l4 - 1.0);
    shape[3][4] = 4.0 * l1 * l2;
    shape[3][5] = 4.0 * l2 * l3;
    shape[3][6] = 4.0 * l3 * l1;
    shape[3][7] = 4.0 * l1 * l4;
    shape[3][8] = 4.0 * l3 * l4; // OpenSees node 9: edge 3-4
    shape[3][9] = 4.0 * l2 * l4; // OpenSees node 10: edge 2-4

    const double derivative[NumNodes][3] = {
        {4.0 * l1 - 1.0, 0.0, 0.0},
        {0.0, 4.0 * l2 - 1.0, 0.0},
        {0.0, 0.0, 4.0 * l3 - 1.0},
        {1.0 - 4.0 * l4, 1.0 - 4.0 * l4, 1.0 - 4.0 * l4},
        {4.0 * l2, 4.0 * l1, 0.0},
        {0.0, 4.0 * l3, 4.0 * l2},
        {4.0 * l3, 0.0, 4.0 * l1},
        {4.0 * (l4 - l1), -4.0 * l1, -4.0 * l1},
        {-4.0 * l3, -4.0 * l3, 4.0 * (l4 - l3)},
        {-4.0 * l2, 4.0 * (l4 - l2), -4.0 * l2}
    };

    // J(a,b) = dx_a / dxi_b, where a is a physical direction and b is a
    // natural-coordinate direction.
    double jacobian[3][3] = {{0.0, 0.0, 0.0},
                             {0.0, 0.0, 0.0},
                             {0.0, 0.0, 0.0}};
    for (int physical = 0; physical < 3; ++physical) {
        for (int natural = 0; natural < 3; ++natural) {
            for (int node = 0; node < NumNodes; ++node) {
                jacobian[physical][natural] +=
                    coordinates[physical][node] * derivative[node][natural];
            }
        }
    }

    const double a = jacobian[0][0];
    const double b = jacobian[0][1];
    const double c = jacobian[0][2];
    const double d = jacobian[1][0];
    const double e = jacobian[1][1];
    const double f = jacobian[1][2];
    const double g = jacobian[2][0];
    const double h = jacobian[2][1];
    const double i = jacobian[2][2];

    determinant = a * (e * i - f * h) - b * (d * i - f * g)
                  + c * (d * h - e * g);
    if (determinant == 0.0) {
        for (int direction = 0; direction < 3; ++direction) {
            for (int node = 0; node < NumNodes; ++node) {
                shape[direction][node] = 0.0;
            }
        }
        return;
    }

    // inverseJ(natural,physical) = dxi_natural / dx_physical
    const double inverseJ[3][3] = {
        {(e * i - f * h) / determinant, (c * h - b * i) / determinant,
         (b * f - c * e) / determinant},
        {(f * g - d * i) / determinant, (a * i - c * g) / determinant,
         (c * d - a * f) / determinant},
        {(d * h - e * g) / determinant, (b * g - a * h) / determinant,
         (a * e - b * d) / determinant}
    };

    for (int node = 0; node < NumNodes; ++node) {
        for (int physical = 0; physical < 3; ++physical) {
            shape[physical][node] = 0.0;
            for (int natural = 0; natural < 3; ++natural) {
                shape[physical][node] +=
                    derivative[node][natural] * inverseJ[natural][physical];
            }
        }
    }
}

void
TenNodeTetrahedronUP::pressureShapeFunction(
    const double naturalCoordinates[3],
    double pressureShape[NumPressureNodes],
    double pressureGradient[3][NumPressureNodes]) const
{
    const double l1 = naturalCoordinates[0];
    const double l2 = naturalCoordinates[1];
    const double l3 = naturalCoordinates[2];
    pressureShape[0] = l1;
    pressureShape[1] = l2;
    pressureShape[2] = l3;
    pressureShape[3] = 1.0 - l1 - l2 - l3;

    // The natural coordinates of the ten displacement nodes. Because the
    // quadratic interpolation exactly reproduces r, s, and t, their physical
    // gradients can be recovered from the already-verified Tet10 gradients.
    static const double naturalNodeCoordinates[NumNodes][3] = {
        {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0},
        {0.0, 0.0, 0.0}, {0.5, 0.5, 0.0}, {0.0, 0.5, 0.5},
        {0.5, 0.0, 0.5}, {0.5, 0.0, 0.0}, {0.0, 0.0, 0.5},
        {0.0, 0.5, 0.0}
    };

    double determinant = 0.0;
    double displacementShape[4][NumNodes];
    this->shapeFunction(naturalCoordinates, determinant, displacementShape);

    for (int physical = 0; physical < 3; ++physical) {
        for (int pressureNode = 0; pressureNode < 3; ++pressureNode) {
            pressureGradient[physical][pressureNode] = 0.0;
            for (int displacementNode = 0; displacementNode < NumNodes;
                 ++displacementNode) {
                pressureGradient[physical][pressureNode] +=
                    naturalNodeCoordinates[displacementNode][pressureNode]
                    * displacementShape[physical][displacementNode];
            }
        }
        pressureGradient[physical][3] =
            -pressureGradient[physical][0]
            -pressureGradient[physical][1]
            -pressureGradient[physical][2];
    }
}

const Matrix &
TenNodeTetrahedronUP::computeB(int node,
                               const double shape[4][NumNodes]) const
{
    B.Zero();
    B(0, 0) = shape[0][node];
    B(1, 1) = shape[1][node];
    B(2, 2) = shape[2][node];
    B(3, 0) = shape[1][node];
    B(3, 1) = shape[0][node];
    B(4, 1) = shape[2][node];
    B(4, 2) = shape[1][node];
    B(5, 0) = shape[2][node];
    B(5, 2) = shape[0][node];
    return B;
}

void
TenNodeTetrahedronUP::formResidAndTangent(int tangentFlag,
                                          bool useInitialTangent)
{
    stiffness.Zero();
    residual.Zero();
    this->computeBasis();

    double shape[4][NumNodes];
    for (int point = 0; point < NumSolidGaussPoints; ++point) {
        double determinant = 0.0;
        this->shapeFunction(solidGaussPoints[point], determinant, shape);
        if (std::abs(determinant) <= 1.0e-18) {
            opserr << "TenNodeTetrahedronUP::formResidAndTangent -- "
                   << "singular Jacobian " << determinant << " in element "
                   << this->getTag() << ", integration point " << point + 1
                   << endln;
            return;
        }

        // The established OpenSees/Abaqus Tet10 connectivity may have the
        // opposite natural-coordinate orientation. Keep the signed Jacobian
        // in the inverse mapping, but integrate over positive physical volume.
        const double differentialVolume = tetWeight * std::abs(determinant);
        const Vector &stress = materialPointers[point]->getStress();

        for (int nodeI = 0; nodeI < NumNodes; ++nodeI) {
            const Matrix bi(this->computeB(nodeI, shape));
            for (int directionI = 0; directionI < 3; ++directionI) {
                const int row = displacementDOF(nodeI, directionI);
                for (int stressComponent = 0;
                     stressComponent < NumStressComponents;
                     ++stressComponent) {
                    residual(row) += bi(stressComponent, directionI)
                                     * stress(stressComponent)
                                     * differentialVolume;
                }
            }

            if (tangentFlag == 0) {
                continue;
            }

            const Matrix &materialTangent =
                useInitialTangent
                    ? materialPointers[point]->getInitialTangent()
                    : materialPointers[point]->getTangent();
            for (int nodeJ = 0; nodeJ < NumNodes; ++nodeJ) {
                const Matrix bj(this->computeB(nodeJ, shape));
                for (int directionI = 0; directionI < 3; ++directionI) {
                    const int row = displacementDOF(nodeI, directionI);
                    for (int directionJ = 0; directionJ < 3; ++directionJ) {
                        const int column = displacementDOF(nodeJ, directionJ);
                        double entry = 0.0;
                        for (int m = 0; m < NumStressComponents; ++m) {
                            for (int n = 0; n < NumStressComponents; ++n) {
                                entry += bi(m, directionI)
                                         * materialTangent(m, n)
                                         * bj(n, directionJ);
                            }
                        }
                        stiffness(row, column) += entry * differentialVolume;
                    }
                }
            }
        }
    }

    const double *activeBodyForce =
        applyLoad == 0 ? bodyForce : appliedBodyForce;
    const double density = this->mixtureDensity(0);

    // Consistent Tet10 solid body-force vector. The same degree-four rule as
    // the displacement mass integrates N_i*rho*b over the element.
    for (int point = 0; point < NumMassGaussPoints; ++point) {
        double determinant = 0.0;
        this->shapeFunction(massGaussPoints[point], determinant, shape);
        if (std::abs(determinant) <= 1.0e-18) {
            opserr << "TenNodeTetrahedronUP::formResidAndTangent -- "
                   << "singular body-force Jacobian " << determinant
                   << " in element " << this->getTag()
                   << ", integration point " << point + 1 << endln;
            return;
        }
        const double differentialVolume =
            massGaussWeights[point] * std::abs(determinant);
        for (int node = 0; node < NumNodes; ++node) {
            for (int direction = 0; direction < 3; ++direction) {
                residual(displacementDOF(node, direction)) -=
                    differentialVolume * density
                    * shape[3][node] * activeBodyForce[direction];
            }
        }
    }

    // Fluid body-force contribution to the continuity equation:
    // rho_f * integral(grad(Np)^T k b)dV.
    double pressureShape[NumPressureNodes];
    double pressureGradient[3][NumPressureNodes];
    for (int point = 0; point < NumSolidGaussPoints; ++point) {
        double determinant = 0.0;
        this->shapeFunction(solidGaussPoints[point], determinant, shape);
        if (std::abs(determinant) <= 1.0e-18) {
            return;
        }
        this->pressureShapeFunction(solidGaussPoints[point], pressureShape,
                                    pressureGradient);
        const double differentialVolume = tetWeight * std::abs(determinant);
        for (int pressureNode = 0; pressureNode < NumPressureNodes;
             ++pressureNode) {
            double contribution = 0.0;
            for (int direction = 0; direction < 3; ++direction) {
                contribution += pressureGradient[direction][pressureNode]
                                * permeability[direction]
                                * activeBodyForce[direction];
            }
            residual(pressureDOF(pressureNode)) +=
                differentialVolume * fluidDensity * contribution;
        }
    }
}

void
TenNodeTetrahedronUP::formDampingTerms(int)
{
    damping.Zero();
    this->computeBasis();

    double displacementShape[4][NumNodes];
    double pressureShape[NumPressureNodes];
    double pressureGradient[3][NumPressureNodes];

    for (int point = 0; point < NumSolidGaussPoints; ++point) {
        double determinant = 0.0;
        this->shapeFunction(solidGaussPoints[point], determinant,
                            displacementShape);
        if (std::abs(determinant) <= 1.0e-18) {
            opserr << "TenNodeTetrahedronUP::formDampingTerms -- "
                   << "singular Jacobian " << determinant << " in element "
                   << this->getTag() << ", integration point " << point + 1
                   << endln;
            return;
        }

        this->pressureShapeFunction(solidGaussPoints[point], pressureShape,
                                    pressureGradient);
        const double differentialVolume = tetWeight * std::abs(determinant);

        // D_up = -Q and D_pu = -Q^T. The fourth nodal velocity is pore
        // pressure in the OpenSees u-p convention.
        for (int displacementNode = 0; displacementNode < NumNodes;
             ++displacementNode) {
            for (int direction = 0; direction < 3; ++direction) {
                const int displacementIndex =
                    displacementDOF(displacementNode, direction);
                for (int pressureNode = 0;
                     pressureNode < NumPressureNodes; ++pressureNode) {
                    const int pressureIndex = pressureDOF(pressureNode);
                    const double coupling =
                        -differentialVolume
                        * displacementShape[direction][displacementNode]
                        * pressureShape[pressureNode];
                    damping(displacementIndex, pressureIndex) += coupling;
                    damping(pressureIndex, displacementIndex) += coupling;
                }
            }
        }

        // D_pp = -H, with H_ij = integral(grad(Np_i)^T k grad(Np_j)) dV.
        // The permeability tensor is diagonal in the global coordinate system.
        for (int pressureNodeI = 0;
             pressureNodeI < NumPressureNodes; ++pressureNodeI) {
            const int row = pressureDOF(pressureNodeI);
            for (int pressureNodeJ = 0;
                 pressureNodeJ < NumPressureNodes; ++pressureNodeJ) {
                const int column = pressureDOF(pressureNodeJ);
                double hydraulicTerm = 0.0;
                for (int direction = 0; direction < 3; ++direction) {
                    hydraulicTerm +=
                        permeability[direction]
                        * pressureGradient[direction][pressureNodeI]
                        * pressureGradient[direction][pressureNodeJ];
                }
                damping(row, column) -=
                    differentialVolume * hydraulicTerm;
            }
        }
    }

    // Add Rayleigh damping to the solid-skeleton displacement block only.
    // The pressure-storage block in getMass() is part of the u-p time
    // formulation and must not receive mass-proportional Rayleigh damping.
    if (alphaM != 0.0) {
        this->formInertiaTerms(1);
        for (int nodeI = 0; nodeI < NumNodes; ++nodeI) {
            for (int nodeJ = 0; nodeJ < NumNodes; ++nodeJ) {
                for (int directionI = 0; directionI < 3; ++directionI) {
                    const int row = displacementDOF(nodeI, directionI);
                    for (int directionJ = 0; directionJ < 3; ++directionJ) {
                        const int column =
                            displacementDOF(nodeJ, directionJ);
                        damping(row, column) += alphaM * mass(row, column);
                    }
                }
            }
        }
    }

    if (betaK != 0.0) {
        damping.addMatrix(1.0, this->getTangentStiff(), betaK);
    }
    if (betaK0 != 0.0) {
        damping.addMatrix(1.0, this->getInitialStiff(), betaK0);
    }
    if (betaKc != 0.0 && Kc != nullptr) {
        damping.addMatrix(1.0, *Kc, betaKc);
    }
}

void
TenNodeTetrahedronUP::formInertiaTerms(int)
{
    mass.Zero();
    this->computeBasis();

    double displacementShape[4][NumNodes];
    double pressureShape[NumPressureNodes];
    double pressureGradient[3][NumPressureNodes];

    // Consistent Tet10 displacement mass. The degree-four Keast rule avoids
    // the rank-deficient mass produced by the four-point stiffness rule.
    const double density = this->mixtureDensity(0);
    for (int point = 0; point < NumMassGaussPoints; ++point) {
        double determinant = 0.0;
        this->shapeFunction(massGaussPoints[point], determinant,
                            displacementShape);
        if (std::abs(determinant) <= 1.0e-18) {
            opserr << "TenNodeTetrahedronUP::formInertiaTerms -- "
                   << "singular mass Jacobian " << determinant
                   << " in element " << this->getTag()
                   << ", integration point " << point + 1 << endln;
            return;
        }
        const double differentialMass =
            density * massGaussWeights[point] * std::abs(determinant);
        for (int nodeI = 0; nodeI < NumNodes; ++nodeI) {
            for (int nodeJ = 0; nodeJ < NumNodes; ++nodeJ) {
                const double entry = differentialMass
                                     * displacementShape[3][nodeI]
                                     * displacementShape[3][nodeJ];
                for (int direction = 0; direction < 3; ++direction) {
                    mass(displacementDOF(nodeI, direction),
                         displacementDOF(nodeJ, direction)) += entry;
                }
            }
        }
    }

    for (int point = 0; point < NumSolidGaussPoints; ++point) {
        double determinant = 0.0;
        this->shapeFunction(solidGaussPoints[point], determinant,
                            displacementShape);
        if (std::abs(determinant) <= 1.0e-18) {
            opserr << "TenNodeTetrahedronUP::formInertiaTerms -- "
                   << "singular Jacobian " << determinant << " in element "
                   << this->getTag() << ", integration point " << point + 1
                   << endln;
            return;
        }

        this->pressureShapeFunction(solidGaussPoints[point], pressureShape,
                                    pressureGradient);
        const double differentialVolume = tetWeight * std::abs(determinant);

        // M_pp = -S, where S_ij = integral(Np_i Np_j / Kc) dV.
        for (int pressureNodeI = 0;
             pressureNodeI < NumPressureNodes; ++pressureNodeI) {
            const int row = pressureDOF(pressureNodeI);
            for (int pressureNodeJ = 0;
                 pressureNodeJ < NumPressureNodes; ++pressureNodeJ) {
                const int column = pressureDOF(pressureNodeJ);
                mass(row, column) -=
                    differentialVolume
                    * pressureShape[pressureNodeI]
                    * pressureShape[pressureNodeJ]
                    / combinedBulkModulus;
            }
        }
    }
}

double
TenNodeTetrahedronUP::mixtureDensity(int) const
{
    // Match the established 20_8_BrickUP convention: NDMaterial::getRho()
    // supplies the bulk density used for solid inertia.
    return materialPointers[0]->getRho();
}

int
TenNodeTetrahedronUP::setParameter(const char **argv, int argc,
                                   Parameter &parameter)
{
    if (argc < 1) {
        return -1;
    }

    // Match the established SSPbrickUP permeability parameter names.
    if (std::strcmp(argv[0], "xPerm") == 0) {
        return parameter.addObject(1, this);
    }
    if (std::strcmp(argv[0], "yPerm") == 0) {
        return parameter.addObject(2, this);
    }
    if (std::strcmp(argv[0], "zPerm") == 0) {
        return parameter.addObject(3, this);
    }

    int result = -1;
    for (int point = 0; point < NumSolidGaussPoints; ++point) {
        const int materialResult =
            materialPointers[point]->setParameter(argv, argc, parameter);
        if (materialResult != -1) {
            result = materialResult;
        }
    }
    return result;
}

int
TenNodeTetrahedronUP::updateParameter(int parameterID, Information &info)
{
    if (parameterID >= 1 && parameterID <= 3) {
        permeability[parameterID - 1] = info.theDouble;
        return 0;
    }

    int result = -1;
    for (int point = 0; point < NumSolidGaussPoints; ++point) {
        const int materialResult =
            materialPointers[point]->updateParameter(parameterID, info);
        if (materialResult != -1) {
            result = materialResult;
        }
    }
    return result;
}
