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

// $Revision: 1.02 $
// $Date: 2021/11/18 00:00:00 $
// Written: Hanlin Dong, Xijun Wang, Tongji University, self@hanlindong.com
//
// Description: This file contains the class definition for DowelType.
// DowelType provides the abstraction of a dowel-type timber joint, including nails, screws, and bolts.
// The envelope curve can be chosen from (1) exponential, (2) Bezier, and (3) piecewise linear.
// The hysteretic law considers force intercept variation, strength degradation,
// unloading, pinching, and reloading stiffness degradation, etc.
// 
// Material definition tcl command:
// uniaxialMaterial DowelType $tag
//     $Fi $Kp $Ru $c $beta $gamma $eta $Dy $alpha_p $alpha_u $alpha_r 
//     <-exponential $K0 $R1 $F0 $Dc $Kd <$Du> <$K0N $R1N $F0N $DcN $KdN $DuN> >
//     <-bezier $D1 $F1 $D2 $F2 $Dc $Fc $Kd <$Du> <$D1N $F1N $D2N $F2N $DcN $FcN $KdN $DuN>> >
//     <-piecewise $D1 $F1 $D2 $F2 $D3 $F3 <$D4 $F4 ...>>

#ifndef DowelType_h
#define DowelType_h

#include <UniaxialMaterial.h>

class DowelType : public UniaxialMaterial
{
  public:
    // Construction method for exponential envelope
    DowelType(int tag, 
              double Fi, double Kp, double Ru, double c,
              double beta, double gamma, double eta,
              double Dy, double Ap, double Au, double Ar,
              double K0_p, double R1_p, double F0_p, double Dc_p, double Kd_p, double Du_p,
              double K0_n, double R1_n, double F0_n, double Dc_n, double Kd_n, double Du_n);
    // Construction method for bezier envelope
    DowelType(int tag, 
              double Fi, double Kp, double Ru, double c,
              double beta, double gamma, double eta,
              double Dy, double Ap, double Au, double Ar,
              double D1_p, double F1_p, double D2_p, double F2_p, 
              double Dc_p, double Fc_p, double Kd_p, double Du_p,
              double D1_n, double F1_n, double D2_n, double F2_n, 
              double Dc_n, double Fc_n, double Kd_n, double Du_n);
    // Construction method for piecewise envelope
    DowelType(int tag, 
              double Fi, double Kp, double Ru, double c,
              double beta, double gamma, double eta,
              double Dy, double Ap, double Au, double Ar,
              int Size, double *Denvs, double *Fenvs);
    // Construction method for copied material
    DowelType();
    ~DowelType();

    const char *getClassType(void) const { return "DowelType"; };

    // Rewrite base class methods
    int setTrialStrain(double strain, double strainRate = 0.0);
    double getStrain(void);
    double getStress(void);
    double getTangent(void);
    double getInitialTangent(void);
    int commitState(void);
    int revertToLastCommit(void);
    int revertToStart(void);
    UniaxialMaterial *getCopy(void);
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel,
        FEM_ObjectBroker &theBroker);
    void Print(OPS_Stream &s, int flag =0);

  protected:

  private:
    // 11 parameters for hysteresis
    double fi;      // reference intercept of pinching line
    double kp;      // pinching line stiffness
    double ru;      // reference unloading stiffness factor
    double c;       // curvature factor
    double beta;    // strength degradation factor
    double gamma;   // energy factor
    double eta;     // intercept degradation factor
    double dyield;  // yield displacement
    double alpha_p; // pinching stiffness degradation factor
    double alpha_u; // unloading stiffness degradation factor
    double alpha_r; // reloading stiffness degradation factor
    
    // common variables for the envelope
    int envType;    // type of envelope curve. 1-exp, 2-bez, 3-pw
    double k0_p;      // initial stiffness
    double k0_n;      // initial stiffness
    double dcap_p;    // envelope cap displacement
    double dcap_n;    // envelope cap displacement
    double fcap_p;    // envelope cap force
    double fcap_n;    // envelope cap force
    double fyield_p;  // yield force
    double fyield_n;  // yield force
    double dult_p;    // ultimate displacement
    double dult_n;    // ultimate displacement
    double kdesc_p;   // descending stiffness
    double kdesc_n;   // descending stiffness

    // other variables to initialize in construction methods
    double dinter_p;  // intersection of pinching and ascending env.
    double dinter_n;  // intersection of pinching and ascending env.
    double eMono_p;   // energy for monotonic loading
    double eMono_n;   // energy for monotonic loading
    
    // extra parameters for exponential envelope
    double k1_p;      // asymptotes stiffness
    double k1_n;      // asymptotes stiffness
    double f0_p;      // asymptotes force intercept
    double f0_n;      // asymptotes force intercept

    // extra parameters for bezier envelope
    double denv1_p;   // first controlling point displacement
    double denv1_n;   // first controlling point displacement
    double fenv1_p;   // first controlling point force
    double fenv1_n;   // first controlling point force
    double denv2_p;   // second controlling point displacement
    double denv2_n;   // second controlling point displacement
    double fenv2_p;   // second controlling point force
    double fenv2_n;   // second controlling point force

    // extra parameters for piecewise envelope
    int envSize;    // size of points
    int envZero;    // zero index of envelope
    double *denvs;  // pointer tp displacement of interp points
    double *fenvs;  // pointer to force of interp points

    // energy related
    bool isPHC;      // if current path is PHC.
    double ePHC_p;   // total energy of primary half-cycles
    double ePHC_n;   // total energy of primary half-cycles
    double eFHC_p;   // total energy of following half-cycles
    double eFHC_n;   // total energy of following half-cycles   

    // pinching controlling points
    double pxs[20], pys[20];     // pinching curve controlling points

    // physical params
    double tStrain;
    double tStress;
    double tTangent;
    double cStrain;
    double cStress;
    double cTangent;

    // loading history params
    int tPath;     // trial loading path
    double tDmin;  // trial minimum displacement in history
    double tFdmin;  // trial force corresponding to the minimum displacement in history    
    double tDmax;  // trial maximum displacement in history
    double tFdmax;  // trial force corresponding to the maximum displacement in history    
    int cPath;     // committed loading path
    double cDmin;  // committed minimum displacement in history
    double cFdmin;  // committed force corresponding to the minimum displacement in history
    double cDmax;  // committed maximum displacement in history
    double cFdmax;  // committed force corresponding to the maximum displacement in history

    // private methods
    double envelope(double disp);             // envelope curve
    double denvelope(double disp);            // envelope tangent
    double envIntersection(double k, double b);  // solve intersection of ascending envelope and a line (positive disp)
    double envWithSlope(double k, bool flag, double xinit);  // solve envelope point with a given slope. Search from xinit.
    // flag is true: right->left. flag is false: left->right
    void resetReversePoints(double disp, double force, bool flag);  // calculate the hysteretic controlling points.
    void getReverseYK(double x, bool flag, double *y, double *k);  // get force and stiffness on hysteretic curve
    double getBezierYK(double x1, double x2, double x3, double x4,
                     double y1, double y2, double y3, double y4,
                     double x, double *other=NULL, bool returnY=true);  //Get Y and K on a Bezier curve given control points an x
};


#endif



