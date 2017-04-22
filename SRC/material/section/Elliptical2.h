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
                                                                        
#ifndef Elliptical2_h
#define Elliptical2_h

#include <SectionForceDeformation.h>

#include <Matrix.h>
#include <Vector.h>
#include <ID.h>

class Elliptical2 : public SectionForceDeformation
{
  public:
    Elliptical2(int tag, double E1, double E2, double sigY1, double sigY2,
		double Hiso, double Hkin1, double Hkin2,
		int c1 = SECTION_RESPONSE_MZ, int c2 = SECTION_RESPONSE_VY);
    Elliptical2();
    ~Elliptical2();

    const char *getClassType(void) const {return "Elliptical2";};

    int setTrialSectionDeformation(const Vector &v);
    const Matrix &getSectionTangent(void);
    const Matrix &getInitialTangent(void);
    const Vector &getStressResultant(void);
    const Vector &getSectionDeformation(void);

    int commitState(void);
    int revertToLastCommit(void);
    int revertToStart(void);
    
    SectionForceDeformation *getCopy(void);
    const ID &getType(void);
    int getOrder(void) const;
    
    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);

    void Print(OPS_Stream &s, int flag = 0);

    Response *setResponse(const char **argv, int argc, OPS_Stream &output);
    int getResponse(int responseID, Information &eleInformation);

    int setParameter(const char **argv, int argc, Parameter &param);
    int updateParameter(int parameterID, Information &info);
    int activateParameter(int paramID);
    const Vector& getStressResultantSensitivity(int gradIndex,
						bool conditional);
    int commitSensitivity(const Vector &dedh, int gradIndex, int numGrads);

  protected:

  private:
	double E[2];
	double sigY[2];
	double Hiso;
	double Hkin[2];

	double e_n1[2];
	double eP_n[2];
	double eP_n1[2];

	double alpha_n;
	double alpha_n1;

	double dg_n1;

	int code1, code2;
	
	int parameterID;
	Matrix *SHVs;

	static Vector s;
	static Matrix ks;
	static ID code;

	// private functions, as per Matlab implementation by E. Taciroglu
	void isotropicHardeningRule(double H, double alpha, double &k, double &dk) {k = 1.0 + H*alpha; dk = H; return;}
	void yieldFunction(const Vector &stress, const double alpha,
			   double &f, double *dfds, double &dfdalpha);
	void computeConsistentTangent(const double dfds[]);
};


#endif
