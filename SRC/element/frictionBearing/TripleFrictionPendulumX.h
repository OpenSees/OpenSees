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

// $Revision$
// $Date$
// $URL$

#ifndef TripleFrictionPendulumX_h
#define TripleFrictionPendulumX_h

// Header file for TripleFrictionPendulumX element
// Extended from TripleFrictionPendulum element 
// Written by Hyunmyung Kim (hkim59@buffalo.edu) and Michael C. Constantinou (constan1@buffalo.edu)
// Created: 2021/11
// Last update: 2023/05
// Version: 2.0

#include <Element.h>
#include <Matrix.h>
#include <Vector.h>

class Node;
class Channel;
class FrictionModel;
class UniaxialMaterial;

class TripleFrictionPendulumX : public Element
{
public:
    // constructors
    TripleFrictionPendulumX(int tag,
        int Nd1, int Nd2, int tag1, int tag2,
        UniaxialMaterial** theMaterials,
        int kpFactor, int kTFactor, int kvFactor,
        double Mu_ref1,
        double Mu_ref2,
        double Mu_ref3,
        double L1,
        double L2,
        double L3,
        double Ubar1,
        double Ubar2,
        double Ubar3,
        double B1,
        double B2,
        double B3,
        double PLATE2,
        double PLATE3,
        double W,
        double Uy,
        double Kvt,
        double minFv,
        double tol,
        double refPressure1,
        double refPressure2,
        double refPressure3,
        double Diffusivity,
        double Conductivity,
        double Temperature0,
        double rateParam,
        double tempParam,
        double unit);



    TripleFrictionPendulumX();

    // destructor
    ~TripleFrictionPendulumX();

    // method to get class type
    const char* getClassType() const { return "TripleFrictionPendulumX"; };

    // public methods to obtain information about dof & connectivity
    int getNumExternalNodes() const;
    const ID& getExternalNodes();
    Node** getNodePtrs();
    int getNumDOF();
    void setDomain(Domain* theDomain);

    // public methods to set the state of the element
    int commitState();
    int revertToLastCommit();
    int revertToStart();
    int update();
    Element* getCopy();

    // public methods to obtain stiffness, mass, damping and residual information
    const Matrix& getTangentStiff();
    const Matrix& getInitialStiff();
    const Matrix& getDamp();
    const Matrix& getMass();

    const Vector& getResistingForce();

    // public methods for output
    int sendSelf(int commitTag, Channel& theChannel);
    int recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker);
    int displaySelf(Renderer& theViewer, int displayMode, float fact, const char** modes, int numMode);
    void Print(OPS_Stream& s, int flag = 0);

    Response* setResponse(const char** argv, int argc, OPS_Stream& s);
    int getResponse(int responseID, Information& eleInformation);

protected:

private:
    // private member functions - only available to objects of the class
    void CircularElasticGap(Matrix& kj, Vector& fj, double Ej, double Gapj, Vector di);
    void BidirectionalPlastic(Matrix& ki, Vector& fi, Vector& epitmp, Vector& qitmp, double Fyi, double Ei, double Hi, Vector epi, Vector qi, Vector di);
    void Segment(Vector& epitmp, Vector& qitmp, bool& conv, Matrix& kij, Vector& di, Vector epi, Vector qi, Vector f, Vector df, double Fyi, double Ei, double Hi, double Ej, double Gapj, double Tol, int Niter);
    void TFPElement(bool& Conv, Vector& ep1tmp, Vector& ep3tmp, Vector& ep5tmp, Vector& q1tmp, Vector& q3tmp, Vector& q5tmp, Matrix& K, Vector& f, Matrix& k12, Matrix& k34, Matrix& k56, Vector& d1, Vector& d3, Vector& d5, Vector ep1, Vector ep3, Vector ep5, Vector q1, Vector q3, Vector q5, Vector u, Vector dusub, double Fy1, double Fy3, double Fy5, double E1, double E3, double E5, double H1, double H3, double H5, double E2, double E4, double E6, double Gap2, double Gap4, double Gap6, double Tol, int Niter);
    void StiffnessForm(Matrix& K, Matrix k12, Matrix k34, Matrix k56);
    double sgn(double x);
    double dTdt_FINITE(double Diffu, double Conduc, double SlabL, double depthzz, double tauu); 
    UniaxialMaterial* theMaterials[4];  // array of uniaxial materials

    double kTF1, kpF1, kvF1; // Dependency factors on COF1
    double kTF2, kpF2, kvF2; // Dependency factors on COF2
    double kTF3, kpF3, kvF3; // Dependency factors on COF3
    double rContact1, rContact2, rContact3; // "B/2"
    double refPressure1, refPressure2, refPressure3; // Reference pressures
    double trialP1, trialP2, trialP3; // Pressure on the sliding surface
    double trialDisp1, trialDisp2, trialDisp3; // Displacements on the sliding surface
    double trialVel1, trialVel2, trialVel3; // Velocities on the sliding surface
    double Mu_ref1, Mu_ref2, Mu_ref3;       // Reference COFs
    double Mu_Adj1, Mu_Adj2, Mu_Adj3;           // current coefficient of friction (COF)
    double Mu_Adj1pr, Mu_Adj2pr, Mu_Adj3pr;           // stored coefficient of friction (COF)
    double Mu_Adj1ppr, Mu_Adj2ppr, Mu_Adj3ppr;           // stored coefficient of friction (COF)
    double trialT;
    double Diffusivity;
    double Conductivity;
    double Temperature0;
    double Temperature_Surface1, Temperature_Surface2, Temperature_Surface3;
    double Temperature_Change1, Temperature_Change2, Temperature_Change3;
    double Temperature_Depth2, Temperature_Depth3;
    double Temperature_Depth_Change2, Temperature_Depth_Change3;

    double DtAnalysis;
    double tau;
    double PI;
    int iCountTime;

    double p_Unit_Convert; //To convert to MPa for pressure factor
    double v_Unit_Convert; //To convert velocity for velocity factor
    double initialTemperature; // Initial temperature

    double disp1, disp2, disp3; // Displacements used for calculating temperature
    double vel1, vel2, vel3; // Velocities used for calculating temerature
    double disp1pr, disp2pr, disp3pr;
    double vel1pr, vel2pr, vel3pr;
    double v23sumpr, u23sumpr;

    Vector DomainTime; //Current time in the domain	
    Vector DomainTimeTemp; //To be used to save data while resizing DomainTime
    Vector DomainHeatFlux1, DomainHeatFlux2, DomainHeatFlux3; //History of heat flux
    Vector DomainHeatFluxTemp1, DomainHeatFluxTemp2, DomainHeatFluxTemp3; //To be used to save data while resizing DomainHeatFlux	

    Vector kTFTemp1, kTFTemp2, kTFTemp3;
    Vector kpFTemp1, kpFTemp2, kpFTemp3;
    Vector kvFTemp1, kvFTemp2, kvFTemp3;
    Vector TemperatureCenter1, TemperatureCenter2, TemperatureCenter3;

    Vector TemperatureDepth2, TemperatureDepth3; 
    double dTdt2_surface, dTdt2_depth, dTdt3_surface, dTdt3_depth; 

    Vector HeatFluxCenter1, HeatFluxCenter2, HeatFluxCenter3;

    int kpFactor;		//If friction is pressure dependent 1, 0 otherwise
    int kTFactor;		//If friction is temperature dependent 1, 0 otherwise
    int kvFactor;		//If friction is velocity dependent 1, 0 otherwise
    double rateParam;	//rate parameter for velocity dependence
    double tempParam;   //Parameter for friction model (1: kT = 1/3, 2: kT = 1/2, 3: kT = 2/3)
    double unit;			//Units of force,displacement etc. 1: N,m,s,C; 2: kN,m,s,C; 3: N,mm,s,C; 4: kN,mm,s,C; 5: lb,in,s,C; 6: kip,in,s,C; 7: lb,ft,s,C; 8: kip,ft,s,C           


    int tag1; // temperature dependent friction
    int tag2; // Heat conduction theories (Indefinite or Finite plate thickness)
    double L1;
    double L2;
    double L3;
    double Ubar1;
    double Ubar2;
    double Ubar3;
    double B1;
    double B2;
    double B3;
    double PLATE2;
    double PLATE3;
    double W;
    double Uy;
    double Kvt;
    double MinFv;
    double TOL;
    int Niter;


    Matrix K;
    Matrix Kpr;
    Vector f;
    Vector fpr;
    Matrix k12;
    Matrix k12pr;
    Matrix k34;
    Matrix k34pr;
    Matrix k56;
    Matrix k56pr;
    Vector d1;
    Vector d1pr;
    Vector d1ppr; //(for velocity calc)
    Vector d3;
    Vector d3pr;
    Vector d3ppr;//(for velocity calc)
    Vector d5;
    Vector d5pr;
    Vector d5ppr;//(for velocity calc)
    Vector v1;
    Vector v3;
    Vector v5;

    Vector ep1;
    Vector ep1pr;
    Vector ep3;
    Vector ep3pr;
    Vector ep5;
    Vector ep5pr;
    Vector q1;
    Vector q1pr;
    Vector q3;
    Vector q3pr;
    Vector q5;
    Vector q5pr;
    Vector ep1tmp;
    Vector ep3tmp;
    Vector ep5tmp;
    Vector q1tmp;
    Vector q3tmp;
    Vector q5tmp;

    double v1Fact, v3Fact, v5Fact;
    double Gap2, Gap4, Gap6;
    double Vel1Avg, Vel3Avg, Vel5Avg;
    double u23xx_storedpr, u23yy_storedpr, D1prAvg, D3prAvg, D5prAvg; // (for velocity calc)
    double Disp1Avg, Disp3Avg, Disp5Avg;
    double Vx, Vy, Vz;

    double Dx_stored, Dy_stored;
    double dispSlope_x, dispSlope_y;

    double uxx, uyy;
    double u23yy, u23xx, u23sum;
    double v23yy, v23xx, v23sum;
    double u23yy_stored, u23xx_stored;




    // Parameters for displacement and velocity histories
    double fx, fy; // Forces
    double u_star, u_star2, udr1, udr4; // Displacement limit states
    double F_f1, F_f2, F_f4, F_dr1, F_dr4; // Force limit states

    // Signs for capturing loading and unloading phases
    double forceSlope_x, forceSlope_y; // To capture the change of phases
    double forceSlope_x_stored, forceSlope_y_stored;

    double sign_fx, sign_fy, changeSignX, changeSignY; // To capture the first direction of loading 
    double loading_x, loading_y, unloading_x, unloading_y; // Tag for not allowing to go back to first cycle

    // Displacement parameters in the loop
    // X direction
    double u1_ref, u4_ref, F_ref; // Reference point when phase change
    double u1_tr, u4_tr, u1_tr_u, u4_tr_u, F_tr, F_tr_u; // Stored points for reference point
    double u1, u4, u1_stored, u4_stored, u1_stored2, u4_stored2; // Displacement histories

    double u2_ref, u3_ref; // Reference point when phase change
    double u2_tr, u2_tr_u, u3_tr, u3_tr_u; // Stored points for reference point
    double u2, u2_stored, u3, u3_stored, u2_stored2, u3_stored2; // Displacement histories
    double u23, v23;//for series element


    // Y direction
    double u1y_ref, u4y_ref, Fy_ref; // Reference point when phase change
    double u1y_tr, u4y_tr, u1y_tr_u, u4y_tr_u, Fy_tr, Fy_tr_u; // Stored points for reference point
    double u1y, u4y, u1y_stored, u4y_stored, u1y_stored2, u4y_stored2; // Displacement histories

    double u2y_ref, u3y_ref; // Reference point when phase change
    double u2y_tr, u2y_tr_u, u3y_tr, u3y_tr_u; // Stored points for reference point
    double u2y, u2y_stored, u3y, u3y_stored, u2y_stored2, u3y_stored2; // Displacement histories
    double u23y, v23y;//for series element

    // Velocity
    double v1x, v2x, v3x, v4x; // Velocity histories
    double v1y, v2y, v3y, v4y;

    // Resultant histories
    double u23t, u1t, u4t;
    double v23t, v1t, v4t;

    double Fy1pr, Fy3pr, Fy5pr;
    double Wpr, Wcr, Wavg;
    double Fy1, Fy3, Fy5;
    double E1, E3, E5;
    double E2, E4, E6;
    double H1, H3, H5;
    double Fvert, Kvert;
    double TorqX, KrotX;
    double TorqY, KrotY;
    double TorqZ, KrotZ;
    double Hisolator;
    double Dx, Dy, Dz;
    bool Conv;

    // private attributes - a copy for each object of the class
    ID  externalNodes;  // contains the id's of end nodes
    Matrix trans;       // hold the transformation matrix, could use a Vector
                        // if ever bother to change the Vector interface for
                        // x-product.

    Node* theNodes[2];  // node pointers

    // static data - single copy for all objects of the class
    static Matrix eleK;      // class wide matrix for returning stiffness
    static Matrix eleKinit;  // class wide matrix for returning initial stiffness
    static Matrix eleD;      // class wide matrix for returning damping
    static Matrix eleM;      // class wide matrix for returning mass
    static Vector eleR;      // class wide vector for returning residual
};

#endif
