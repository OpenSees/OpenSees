/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 2001, The Regents of the University of California    **
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
** Reliability module developed by:                                   **
**   Terje Haukaas (haukaas@ce.berkeley.edu)                          **
**   Armen Der Kiureghian (adk@ce.berkeley.edu)                       **
**                                                                    **
** ****************************************************************** */
                                                                        
// $Revision: 1.1 $
// $Date: 2003-03-06 18:34:14 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/BoucWenMaterial.h,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu) 
//

#ifndef BoucWenMaterial_h
#define BoucWenMaterial_h

#include <UniaxialMaterial.h>
#include <Matrix.h>

class BoucWenMaterial : public UniaxialMaterial
{
  public:
    BoucWenMaterial(int tag, 
					double alpha,
					double ko,
					double n,
					double gamma,
					double beta,
					double Ao,
					double deltaA,
					double deltaNu,
					double deltaEta,
					double tolerance,
					int maxNumIter);
	BoucWenMaterial();
    ~BoucWenMaterial();
    int setTrialStrain(double strain, double strainRate = 0.0); 
    double getStrain(void);          
    double getStress(void);
    double getTangent(void);
    double signum(double);
    int commitState(void);
    int revertToLastCommit(void);    
    int revertToStart(void);        
    UniaxialMaterial *getCopy(void);
    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);    
    void Print(OPS_Stream &s, int flag =0);
    
	// Reliability and sensitivity stuff
    double getInitialTangent        (void);
    int    setParameter             (const char **argv, int argc, Information &info);
    int    updateParameter          (int parameterID, Information &info);
	int    activateParameter        (int parameterID);
	double getStressSensitivity     (int gradNumber, bool conditional);
	double getStrainSensitivity     (int gradNumber);
	double getTangentSensitivity    (int gradNumber);
	double getDampTangentSensitivity(int gradNumber);
	double getRhoSensitivity        (int gradNumber);
	int    commitSensitivity        (double strainGradient, int gradNumber, int numGrads);

  protected:
    
  private:
    // Material parameters
    double alpha;
    double ko;
	double n;
    double gamma;
    double beta;
    double Ao;
    double deltaA;
    double deltaNu;
    double deltaEta;

    // History variables (trial and commited)
    double Tstrain, Cstrain;
	double Tz, Cz;
	double Te, Ce;

	// Ohter variables
	double Tstress, Ttangent;

	double tolerance;
	int maxNumIter;

	// Sensitivit stuff
    int parameterID;
	Matrix *SHVs;


};


#endif

