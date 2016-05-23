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

// Written: Saumyasuchi Das
//          May 2013, University of Canterbury
// Updated: Chris McGann
//          June 2015, Washington State University

#ifndef StressDensityModel_h
#define StressDensityModel_h

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <NDMaterial.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>

class StressDensityModel : public NDMaterial
{
  public:

    // full constructor
	StressDensityModel(int tag, int classTag, double constDensity,
					   // SD model  parameters		
					   double initialVoidRatio, double constA, double exponentN,
                       double poissonRatio, double constAlpha1, double constBeta1,	
                       double constAlpha2, double constBeta2, double constAlpha3,
                       double constBeta3, double constDegradation, double constMumin,
                       double constMucyclic, double constDilatancyStrain,	
                       double constMumax, double constPatm, 
                       // steady state line void ratio
                       double constsslvoidatP1=0.877, double constsslvoidatP2=0.877, double constsslvoidatP3=0.873,
                       double constsslvoidatP4=0.870, double constsslvoidatP5=0.860, double constsslvoidatP6=0.850,
                       double constsslvoidatP7=0.833, double constsslvoidatP8=0.833, double constsslvoidatP9=0.833,
                       double constsslvoidatP10=0.833,
                       // hydrostatic state line void ratio
                       double consthslvoid=0.895, 
                       // reference pressures
                       double constP1=1.0, double constP2=10.0, double constP3=30.0, double constP4=50.0, 
                       double constP5=100.0, double constP6=200.0, double constP7=400.0, double constP8=400.0, 
                       double constP9=400.0, double constP10=400.0,
                       // offset of the failure surface
                       double constRxx=0.0, double constRyy=0.0, double constRzz=0.0,
                       double constRxy=0.0, double constRyz=0.0, double constRzx=0.0);

    // null constructor
    StressDensityModel();

    // destructor
	~StressDensityModel();

    NDMaterial *getCopy(const char *type);

    virtual int commitState(void);
    int revertToLastCommit(void);
    int revertToStart(void);

    NDMaterial *getCopy(void);
    const char *getType(void) const;
    int getOrder(void) const;

    Response *setResponse (const char **argv, int argc, OPS_Stream &output);
    int getResponse (int responseID, Information &matInformation);

    virtual int sendSelf(int commitTag, Channel &theChannel);  
    virtual int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker); 

    void Print(OPS_Stream &s, int flag =0);

	int setParameter(const char **argv, int argc, Parameter &param);
  	int updateParameter(int responseID, Information &eleInformation);

	// send mass density to element in dynamic analysis
	double getRho(void);

protected:

    double theDensity;
    double modelParameter[16];
    double sslVoidratio[10];
    double hslVoidratio[10];
    double refPressure[10];
    double refOrigin[6];

	int hard_para_int[2];
	int theStage;

	static const int Nsurface = 35;
};

#endif
