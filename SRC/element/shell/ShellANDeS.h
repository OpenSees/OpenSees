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
                                                                        
// $Revision: 1.19 $
// $Date: 2009-08-25 22:32:08 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/Element.h,v $
                                                                        
                                                                        
// Written: jaabell 
// Created: 06/2017
// Revision: A
//
// What: "@(#) Element.h, revA"
/*

Element is based on the following references by Carlos Felippa

Membrane and drilling parts:
(i)   Membrane Triangles with Corner Drilling Freedoms Part I: The EFF Element.
        Ken Alvin, Horacio M. de la Fuente, Bjorn Haugen, & Carlos A. Felippa.
        August 1991,
        University of Colorado, Boulder,
        Report No. CU-CDDC-91-24
(ii)  Membrane Triangles with Corner Drilling Freedoms Part II: The ANDES Element.
        Carlos A. Felippa & Carmello Militello,
        August 1991,
        University of Colorado, Boulder,
        Report No. CU-CDDC-91-24b
(iii) Membrane Triangles with Corner Drilling Freedoms Part III: Implementation and performance evaluation.
        Carlos A. Felippa & Scott Alexander,
        August 1991,
        University of Colorado, Boulder,
        Report No. CU-CDDC-91-24c
(iv)  A Study of Optimal Triangles with Drilling Freedoms.
        Carlos A. Felippa
        February 2003
        Report No. CU-CAS-03-02

For the bending component
(v)   The First ANDES Elements: 9-DOF Plate Bending Trianges
        Carmello Militello & Carlos A. Felippa
        December 1989
        Report No. CU-CSSC-89-22

Obtainable as of July 22 2012 at http://www.colorado.edu/engineering/CAS/Felippa.d/FelippaHome.d/Home.html

These documents all mirror published works in indexed journals.
*/

#ifndef ShellANDeS_H
#define ShellANDeS_H


#ifndef _bool_h
#include "bool.h"
#endif

#include <Information.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Element.h>
#include <Node.h>
#include <ID.h>
#include <Domain.h>
#include <OPS_Globals.h>
#include <NDMaterial.h>
#include <Matrix.h>
#include <Vector.h>
#include <ElementalLoad.h>
// #include <ElasticMembranePlateSection.h>
// #include <SectionForceDeformation.h>

class Node;

class ShellANDeS: public Element
{

public:

    ShellANDeS ();

    ShellANDeS(int element_number, int node_numb_1, int node_numb_2, int node_numb_3, double t, double E, double nu, double rho);
    ShellANDeS(int element_number, int node_numb_1, int node_numb_2, int node_numb_3, double t, double E11, double E22,
    double E33, double E12, double E13, double E23, double n1, double n2, double n3, double rho);
    
    // ShellANDeS(int element_number,
    //                     int node_numb_1, int node_numb_2, int node_numb_3,
    //                      SectionForceDeformation *elasticmembraneplatsection_ptr); //future

    ~ShellANDeS();

    int getNumExternalNodes () const;
    const ID &getExternalNodes ();
    Node **getNodePtrs(void);
    int getNumDOF ();
    void setDomain(Domain *theDomain);

    int commitState ();
    int revertToLastCommit ();
    int revertToStart ();
    int update(void);


    const Matrix &getTangentStiff ();
    const Matrix &getInitialStiff();
    const Matrix &getMass ();

    void zeroLoad ();
    int addLoad(ElementalLoad *theLoad, double loadFactor);
    int addInertiaLoadToUnbalance(const Vector &accel);

    const Vector &getResistingForce ();
    const Vector &getResistingForceIncInertia ();

    int sendSelf (int commitTag, Channel &theChannel);
    int recvSelf (int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

    void Print(OPS_Stream &s, int flag = 0);
    Response* setResponse (const char** argv, int argc, OPS_Stream& theHandler);
    int getResponse (int responseID, Information& eleInformation);

    int displaySelf(Renderer&, int mode, float fact, const char** displayModes = 0, int numModes = 0);

    Matrix returnMass(void);

    void useThisCoordinateSystem(Vector e1, Vector e2, Vector e3);

    double getArea() const;

    const Vector &getBodyForce(double loadFactor, const Vector &data);
    const Vector FormEquiBodyForce(void);

    //For stress recovery
    const Vector & get_bending_moment_field();

    //Utility
    bool gotMass() const;           // Got Mass?  (check for rho value)

    std::string getElementName() const
    {
        return "ShellANDeS";
    }



private:
    const Matrix &getBendingTangentStiffness();
    const Matrix &getMembraneTangentStiffness();
    const Matrix &getMembraneMass ();
    const Matrix &getBendingMass ();
    
    //Membrane functions
    void initializeGeometry(double n1, double n2, double n3);                              // Calculates all geometry related internal variables
    void calculate_E_planestress_and_beta0();
    void initializeBetaArrays();
    Matrix getMembraneForceLumpingMatrix();
    Matrix getMembraneHierarchicalProjectionMatrix();
    Matrix getMembraneNaturalStrainProjectionMatrix();
    Matrix getMembraneBasicStiffness();
    Matrix getMembraneHighOrderStiffness();

    //Bending functions
    void initializeMq();
    Matrix getBendingBasicStiffness();
    Matrix getBendingHighOrderStiffness();

    inline int sendAndCheckVector(int dataTag , int commitTag , Vector& v , Channel & channel , const std::string &name);
    inline int recvAndCheckVector(int dataTag , int commitTag , Vector& v , Channel & channel , const std::string &name);
    inline int sendAndCheckMatrix(int dataTag , int commitTag , Matrix& v , Channel & channel , const std::string &name);
    inline int recvAndCheckMatrix(int dataTag , int commitTag , Matrix& v , Channel & channel , const std::string &name);
    inline int sendAndCheckID(int dataTag     , int commitTag , ID& v     , Channel & channel , const std::string &name);
    inline int recvAndCheckID(int dataTag     , int commitTag , ID& v     , Channel & channel , const std::string &name);

    //===================================================================================
    // Internal variables
    //===================================================================================

    // Pointers to other Domain components
    ID  connectedExternalNodes;     // Tags of  nodes
    Node *theNodes[3];              // pointers to 3 nodes

    // Finite element matrices
    Matrix K;                       // Element stiffness Matrix
    Matrix M;                       // Element mass matrix
    Vector P;                       // Element resisting force vector
    Vector Q;                       // Applied nodal loads
    Vector bf;                      // Body forces

    // ThreeNodeAndesBending *bending_element;
    // ThreeNodeAndesMembrane *membrane_element;

    static unsigned int number_of_three_node_andes_shells;

    bool is_stiffness_calculated;
    bool is_mass_calculated;


    // Membrane stuff
      // Geometry related
    double thickness;               // Element thickness
    Vector xl1,  xl2, xl3;          // Node i's local x and y coordinates (z = 0)
    Vector x0;                      // Centroid in global coordinates
    Matrix T_lg;                    // Local-to-global transformation matrix T_lg = [x_l, y_l, z_l]
    double Area;                    // yes... its the area of the triangle :)
    double x12, x23, x31;           // Differences in local x coordinates on nodes ij, xij = xli[0] - xlj[0]
    double y12, y23, y31;           // Differences in local y coordinates on nodes ij, yij = yli[0] - ylj[0]


    // Material related
    double rho;
    double mE11, mE22, mE33, mE12, mE13, mE23;
    double n1, n2, n3;              // 1-1 direction
    Matrix E_planestress;           // Plane stress constitutive matrix

    // Element formulation related
    static double alpha_membrane;   // Scale factor relating boundary normal displacements to corner rotations. See paper (1) by Felippa.
    static Vector beta_membrane;    // beta_membrane parameter values
    double beta0;                   // Optimal beta0 depends on material Poisson's ratio as seen in ref (iv) page 17 --> beta0 = 1/2*(1-4*nu^2)

    static unsigned int number_of_three_node_andes_membrane;

    //Bending stuff
   // Element formulation related
    double disp_init[3][6];
    
    static Matrix Mq;

    int initialized_disps;
};


#endif // ShellANDeS_H
