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
// $Date: 2008-02-29 19:43:53 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/telm/NewStaticSensitivityIntegrator.h,v $


#ifndef NewStaticSensitivityIntegrator_h
#define NewStaticSensitivityIntegrator_h

#include <SensitivityIntegrator.h>
#include <StaticIntegrator.h>
class AnalysisModel;

class NewStaticSensitivityIntegrator : public SensitivityIntegrator,
									public StaticIntegrator
{
  public:
    NewStaticSensitivityIntegrator(AnalysisModel *theModel, LinearSOE *theLinSOE);
    ~NewStaticSensitivityIntegrator();
    

	// Methods promised by the ordinary integrator
    int newStep(void);    
    int update(const Vector &deltaU);
    int setDeltaLambda(double newDeltaLambda);
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
    void Print(OPS_Stream &s, int flag =0);    

	// Sensitivity related methods
	int formEleResidual(FE_Element *theEle);
	int formSensitivityRHS(int gradNum);
	int formIndependentSensitivityRHS();
	int saveSensitivity(const Vector &v, int gradNum, int numGrads);
	int commitSensitivity(int gradNum, int numGrads);
	////////////////////////////////////////////////////////////////
	////////////////// Added by K Fujimura /////////////////////////
	////////////////////////////////////////////////////////////////
	int updateGradNumber(int passedGradNumber);
	int sensitivityDomainChanged(int NumGrads);
	bool staticSensitivity(void);
	bool NewSensitivity(void);

  protected:
    
  private:
	int gradNumber;
	AnalysisModel *theAnalysisModel;
	LinearSOE *theSOE;
};

#endif

