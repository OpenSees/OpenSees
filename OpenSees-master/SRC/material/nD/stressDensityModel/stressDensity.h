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

// Written: Chris McGann
//          July 2018, University of Canterbury
//

#ifndef stressDensity_h
#define stressDensity_h

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <NDMaterial.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>

#ifdef _WIN32
#define sdmuc_ SDMUC

extern "C" {void sdmuc_(
            double strhs[],
            double strsg[],
            double props[],
            double stran[],
            int nmats,
            int nstrp,
            int istep,
            int iiter,
            int ielem,
            double strhs0[],
            double etahs[][3],
            double hdp[][3],
            double oths[]);}

#else

extern "C" {void sdmuc_(
            double strhs[],
            double strsg[],
            double props[],
            double stran[],
            int nmats,
            int nstrp,
            int istep,
            int iiter,
            int ielem,
            double strhs0[],
            double etahs[][3],
            double hdp[][3],
            double oths[]);}
#endif

class stressDensity : public NDMaterial
{
  public:

    // full constructor
	stressDensity(int tag, int classTag, double massDen,
			      // SD model  parameters		
				  double eInit, double constA, double expN, double nu, double a1, double b1,	
                  double a2, double b2, double a3, double b3, double fd, double muNot,
                  double muCyc, double sc, double M, double patm, 
                  // steady state line definition (optional input arguments)
                  double ssl1=0.877, double ssl2=0.877, double ssl3=0.873, double ssl4=0.870, double ssl5=0.860, 
                  double ssl6=0.850, double ssl7=0.833, double hsl=0.895, double pmin=1.0);

    // null constructor
    stressDensity();

    // destructor
    ~stressDensity();

    int commitState(void);
    int revertToLastCommit(void);
    int revertToStart(void);

    NDMaterial *getCopy(void);
    NDMaterial *getCopy(const char *type);
    const char *getType(void) const;
    int getOrder(void) const;

    // get trial strain from the element
    int setTrialStrain(const Vector &strain_from_element);
    // unused strain rate function
    int setTrialStrain(const Vector &v, const Vector &r);

    // send back the strain
    const Vector& getStrain();
    // send back the stress
    const Vector& getStress();
    // send back the tangent
    const Matrix& getTangent();
    const Matrix& getInitialTangent();
    // send mass density to element 
    double getRho(void);

    Response *setResponse(const char **argv, int argc, OPS_Stream &output);
    int getResponse(int responseID, Information &matInformation);

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

    void Print(OPS_Stream &s, int flag =0);

    int setParameter(const char **argv, int argc, Parameter &param);
    int updateParameter(int responseID, Information &eleInformation);

  protected:

    // integers to send to fortran
    static const int nmats = 100;
    static const int nstrp = 100;
    int istep;
    int iiter;
    static const int ielem = 1;

    // arrays to send to fortran
    double strhs[nstrp];
    double strsg[4];
    double props[nmats];
    double stran[4];
    double strhs0[280];
    double etahs[40][3];
    double hdp[80][3];
    double oths[12];

    // member variables for stressDensity
    int theStage;
    int pFlag;
    double pInit;
    double massDensity;
    Vector stressCurrent;
    Vector stressNext;
    Vector strainCurrent;
    Vector strainNext;
    Vector materialParam;
    Matrix initialTangent;
    Matrix currentTangent;

    // member functions
    void initialise();
    void calInitialTangent(void);
    void getCurrentStress(void);

};
#endif 
