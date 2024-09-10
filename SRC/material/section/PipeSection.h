/* ******************************************************************
***
**    OpenSees - Open System for Earthquake Engineering Simulation **
**          Pacific Earthquake Engineering Research Center **
** **
** **
** (C) Copyright 1999, The Regents of the University of California **
** All Rights Reserved. **
** **
** Commercial use of this program without express permission of the **
** University of California, Berkeley, is strictly prohibited.  See **
** file 'COPYRIGHT'  in main directory for information on usage and **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES. **
** **
** Developed by: **
**   Frank McKenna (fmckenna@ce.berkeley.edu) **
**   Gregory L. Fenves (fenves@ce.berkeley.edu) **
**   Filip C. Filippou (filippou@ce.berkeley.edu) **
** **
** ******************************************************************
*/

// Minjie

#ifndef PipeSection_h
#define PipeSection_h

#include <Matrix.h>
#include <SectionForceDeformation.h>
#include <Vector.h>
#include <classTags.h>
#include <elementAPI.h>

class Channel;
class FEM_ObjectBroker;
class Information;

class PipeSection : public SectionForceDeformation {
    // line 9939
   public:
    PipeSection(int tag, double d, double t, double a, double r);
    ~PipeSection(void);

    int commitState(void);
    int revertToLastCommit(void);
    int revertToStart(void);

    const char *getClassType(void) const;

    int setTrialSectionDeformation(const Vector &);
    const Vector &getSectionDeformation(void);

    const Vector &getStressResultant(void);
    const Matrix &getSectionTangent(void);
    const Matrix &getInitialTangent(void);
    const Matrix &getSectionFlexibility(void);
    const Matrix &getInitialFlexibility(void);

    SectionForceDeformation *getCopy(void);
    const ID &getType(void);
    int getOrder(void) const;

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel,
                 FEM_ObjectBroker &theBroker);

    void Print(OPS_Stream &s, int flag = 0);

    int setParameter(const char **argv, int argc, Parameter &param);
    int updateParameter(int parameterID, Information &info);
    int activateParameter(int parameterID);
    const Vector &getStressResultantSensitivity(int gradIndex,
                                                bool conditional);
    const Matrix &getSectionTangentSensitivity(int gradIndex);
    const Matrix &getInitialTangentSensitivity(int gradIndex);
    const Matrix &getSectionFlexibilitySensitivity(int gradIndex);
    const Matrix &getInitialFlexibilitySensitivity(int gradIndex);

    double DOUT() const { return dout; }
    double WALL() const { return thk; }
    double ALFAV() const { return alphaV; }
    double RHO() const { return rho; }
    double ROUT() const { return rout; }
    double RIN() const { return rin; }
    double AREA() const { return Ax; }
    double IY() const { return Iy; }
    double IZ() const { return Iz; }
    double JX() const { return Jx; }

   protected:
   private:
    double dout;    // SECTD 9508
    double thk;     // WALL 9507
    double alphaV;  // SHEAR 9500
    double rho;     // SECTM 9509

    // DOUT         =  OUTSIDE DIAMETER
    // WALL         =  WALL THICKNESS
    // ALFAV        =  SHAPE FACTOR FOR SHEAR DISTORTION
    // AREA         =  CROSS SECTIONAL AREA
    // XMI          =  SECTION PRINCIPAL MOMENT OF INERTIA

    Vector e;
    Matrix mat;
    ID code;

    double rout;
    double rin;
    double Ax;
    double Iy, Iz;
    double Jx;
    double Ay, Az;

    int parameterID;
};

#endif
