/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */

// Ten-node displacement / four-node pressure tetrahedral element.
//
// Nodes 1-4  : ux, uy, uz, p
// Nodes 5-10 : ux, uy, uz
// Total      : 34 element DOFs

#ifndef TenNodeTetrahedronUP_H
#define TenNodeTetrahedronUP_H

#include <Element.h>
#include <ID.h>
#include <Matrix.h>
#include <NDMaterial.h>
#include <Vector.h>

class Channel;
class Domain;
class ElementalLoad;
class FEM_ObjectBroker;
class Information;
class Node;
class OPS_Stream;
class Parameter;
class Renderer;
class Response;

class TenNodeTetrahedronUP : public Element
{
public:
    TenNodeTetrahedronUP();

    TenNodeTetrahedronUP(int tag,
                         int node1, int node2, int node3, int node4,
                         int node5, int node6, int node7, int node8,
                         int node9, int node10,
                         NDMaterial &theMaterial,
                         double bulk, double fluidDensity,
                         double permX, double permY, double permZ,
                         double b1 = 0.0, double b2 = 0.0,
                         double b3 = 0.0);

    ~TenNodeTetrahedronUP() override;

    const char *getClassType() const override
    {
        return "TenNodeTetrahedronUP";
    }

    void setDomain(Domain *theDomain) override;

    int getNumExternalNodes() const override;
    const ID &getExternalNodes() override;
    Node **getNodePtrs() override;
    int getNumDOF() override;

    int commitState() override;
    int revertToLastCommit() override;
    int revertToStart() override;
    int update() override;

    const Matrix &getTangentStiff() override;
    const Matrix &getInitialStiff() override;
    const Matrix &getDamp() override;
    const Matrix &getMass() override;

    void zeroLoad() override;
    int addLoad(ElementalLoad *theLoad, double loadFactor) override;
    int addInertiaLoadToUnbalance(const Vector &accel) override;

    const Vector &getResistingForce() override;
    const Vector &getResistingForceIncInertia() override;

    int sendSelf(int commitTag, Channel &theChannel) override;
    int recvSelf(int commitTag, Channel &theChannel,
                 FEM_ObjectBroker &theBroker) override;

    void Print(OPS_Stream &stream, int flag = 0) override;
    int displaySelf(Renderer &theViewer, int displayMode, float fact,
                    const char **modes = nullptr,
                    int numModes = 0) override;

    Response *setResponse(const char **argv, int argc,
                          OPS_Stream &output) override;
    int getResponse(int responseID, Information &elementInfo) override;

    int setParameter(const char **argv, int argc, Parameter &param) override;
    int updateParameter(int parameterID, Information &info) override;

private:
    enum {
        NumNodes = 10,
        NumPressureNodes = 4,
        NumSolidGaussPoints = 4,
        NumMassGaussPoints = 11,
        NumStressComponents = 6,
        NumDisplacementDOF = 30,
        NumPressureDOF = 4,
        NumDOF = 34
    };

    // Tet10 shape functions and global derivatives. The OpenSees ordering is
    // 5:1-2, 6:2-3, 7:3-1, 8:1-4, 9:3-4, 10:2-4.
    void shapeFunction(const double naturalCoordinates[3], double &detJ,
                       double shape[4][NumNodes]) const;

    // Linear pressure interpolation at the four corner nodes.
    void pressureShapeFunction(const double naturalCoordinates[3],
                               double pressureShape[NumPressureNodes],
                               double pressureGradient[3][NumPressureNodes]) const;

    const Matrix &computeB(int node,
                           const double shape[4][NumNodes]) const;

    void computeBasis();
    void formResidAndTangent(int tangentFlag,
                             bool useInitialTangent = false);
    void formInertiaTerms(int tangentFlag);
    void formDampingTerms(int tangentFlag);

    double mixtureDensity(int integrationPoint) const;

    // Map a displacement component (0=x, 1=y, 2=z) to the 34-DOF vector.
    static int displacementDOF(int node, int component);
    static int pressureDOF(int cornerNode);

    ID connectedExternalNodes;
    Node *nodePointers[NumNodes];
    NDMaterial *materialPointers[NumSolidGaussPoints];

    double coordinates[3][NumNodes];
    double bodyForce[3];
    double appliedBodyForce[3];
    int applyLoad;

    double fluidDensity;
    double combinedBulkModulus;
    double permeability[3];

    Vector *load;
    Matrix *initialStiffness;

    static Matrix stiffness;
    static Matrix damping;
    static Matrix mass;
    static Vector residual;
    static Matrix B;
};

#endif // TenNodeTetrahedronUP_H
