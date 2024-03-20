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

// $Revision: 1.10 $
// $Date: 2024/03

// Original implementation: Massimo Petracca (ASDEA)
//
// A 3-node general shear-deformable (thick) shell element based on 
// the Cell-based Smoothed Discrete Shear Gap (CS-DSG) formulation 
// to avoid transverse-shear-locking in the thin shell limit.
// 
// It supports both linear and corotational kinematics.
//
//  References:
//
// 1 Nguyen-Thoi, T., Phung-Van, P., Thai-Hoang, C., & Nguyen-Xuan, H. 
//   "A cell-based smoothed discrete shear gap method (CS-DSG3) using triangular 
//   elements for static and free vibration analyses of shell structures"
//   International Journal of Mechanical Sciences, 74, 32-45 (2013).
//   https://www.sciencedirect.com/science/article/pii/S002074031300129X
//   https://biblio.ugent.be/publication/5662705/file/5752662
//
// 2 Hughes, Thomas JR, and F. Brezzi. 
//   "On drilling degrees of freedom." 
//   Computer methods in applied mechanics and engineering 72.1 (1989): 105-121.
//   https://www.sciencedirect.com/science/article/pii/0045782589901242
//

#ifndef ASDShellT3_h
#define ASDShellT3_h

#include <Element.h>
#include <ID.h>
#include <Vector.h>
#include <Matrix.h>
#include <Damping.h>
#include <vector>

class SectionForceDeformation;
class ASDShellT3Transformation;
class ASDShellT3LocalCoordinateSystem;

class ASDShellT3 : public Element
{
public:
    enum DrillingDOFMode {
        DrillingDOF_Elastic = 0,
        DrillingDOF_NonLinear = 1,
    };
    class NLDrillingData {
    public:
        std::vector<Vector> strain_comm = { Vector(8), Vector(8), Vector(8) };
        std::vector<Vector> stress_comm = { Vector(8), Vector(8), Vector(8) };
        double damage = 0.0;
        double damage_comm = 0.0;
    };

public:

    // life cycle
    ASDShellT3();
    ASDShellT3(
        int tag,
        int node1,
        int node2,
        int node3,
        SectionForceDeformation* section,
        const Vector& local_x,
        bool corotational = false,
        DrillingDOFMode drill_mode = DrillingDOF_Elastic,
        Damping* theDamping = 0);
    virtual ~ASDShellT3();

    const char* getClassType(void) const { return "ASDShellT3"; }

    // domain
    void setDomain(Domain* theDomain);
    int setDamping(Domain* theDomain, Damping* theDamping);

    // print
    void Print(OPS_Stream& s, int flag);

    // methods dealing with nodes and number of external dof
    int getNumExternalNodes() const;
    const ID& getExternalNodes();
    Node** getNodePtrs();
    int getNumDOF();

    // methods dealing with committed state and update
    int commitState();
    int revertToLastCommit();
    int revertToStart();
    int update();

    // methods to return the current linearized stiffness,
    // damping and mass matrices
    const Matrix& getTangentStiff();
    const Matrix& getInitialStiff();
    const Matrix& getMass();

    // methods for applying loads
    void zeroLoad();
    int addLoad(ElementalLoad* theLoad, double loadFactor);
    int addInertiaLoadToUnbalance(const Vector& accel);

    // methods for obtaining resisting force (force includes elemental loads)
    const Vector& getResistingForce();
    const Vector& getResistingForceIncInertia();

    // public methods for element output
    int sendSelf(int commitTag, Channel& theChannel);
    int recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker);

    Response* setResponse(const char** argv, int argc, OPS_Stream& output);
    int getResponse(int responseID, Information& eleInfo);

    int setParameter(const char** argv, int argc, Parameter& param);

    // calculate the characteristic length for this element
    double getCharacteristicLength(void);

    // display
    int displaySelf(Renderer&, int mode, float fact, const char** displayModes = 0, int numModes = 0);

private:

    // internal method to compute everything using switches...
    int calculateAll(Matrix& LHS, Vector& RHS, int options);

    void AGQIinitialize();
    void AGQIupdate(const Vector& UL);
    void AGQIbeginGaussLoop(const ASDShellT3LocalCoordinateSystem& reference_cs);

private:

    // cross section
    SectionForceDeformation* m_section = nullptr;

    // nodal ids
    ID m_node_ids = ID(3);
    Node* nodePointers[3] = { nullptr, nullptr, nullptr };

    // coordinate transformation
    ASDShellT3Transformation* m_transformation = nullptr;

    // vectors for applying load (allocated only if necessary)
    Vector* m_load = nullptr;

    // drilling strain for the independent rotation field (Hughes-Brezzi)
    DrillingDOFMode m_drill_mode = DrillingDOF_Elastic;
    double m_drill_strain = 0.0;
    double m_drill_stiffness = 0.0;
    NLDrillingData* m_nldrill = nullptr;

    // section orientation with respect to the local coordinate system
    Vector* m_local_x = nullptr;
    double m_angle = 0.0;

    // damping
    Damping* m_damping = nullptr;

    // initialization flag
    bool m_initialized = false;
};

#endif // ASDShellT3_h
