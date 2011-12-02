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
                                                                        
// $Revision: 1.2 $
// $Date: 2010-02-04 20:12:57 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/telm/NewNewmarkSensitivityIntegrator.h,v $

#ifndef NewNewmarkSensitivityIntegrator_h
#define NewNewmarkSensitivityIntegrator_h

#include <Newmark.h>
#include <SensitivityIntegrator.h>

class FE_Element;
class DOF_Group;
#include <AnalysisModel.h>
#include <Domain.h>

#include <Information.h>
#include <LoadPattern.h>
#include <LoadPatternIter.h>
#include <Node.h>
#include <NodeIter.h>
#include <LinearSOE.h>

//class FE_Element;
//class DOF_Group;
//class Vector;
//class Information;


class NewNewmarkSensitivityIntegrator : public SensitivityIntegrator , public Newmark
{
  public:
    NewNewmarkSensitivityIntegrator();
    NewNewmarkSensitivityIntegrator(int assemblyFlag, double gamma, double beta, bool disp = true);
    NewNewmarkSensitivityIntegrator(int assemblyFlag, double gamma, double beta, double alphaM, double betaKcurrent,
	    double betaKinit, double betaKlastCommit, bool disp = true);
    ~NewNewmarkSensitivityIntegrator();
    
    int setParameter      (char **argv, int argc, Information &info);
    int updateParameter   (int parameterID, Information &info);
	int activateParameter (int parameterID);

	int formEleResidual(FE_Element *theEle);
	int formNodUnbalance(DOF_Group *theDof);

	int formSensitivityRHS(int gradNum);
	int formIndependentSensitivityRHS();
	int saveSensitivity   (const Vector &v, int gradNum, int numGrads);
	int commitSensitivity (int gradNum, int numGrads);  
	////////////////////////////////////////////////////////////////
	////////////////// Added by K Fujimura /////////////////////////
	////////////////////////////////////////////////////////////////
	int updateGradNumber(int passedGradNumber);
	int sensitivityDomainChanged(int NumGrads);
	bool staticSensitivity(void);
	bool NewSensitivity(void);

  protected:
    
  private:

	LinearSOE *theLinSOE;
	AnalysisModel *theModel;
	Domain *theDomain;
	Node *nodePtr;
	LoadPattern *loadPatternPtr;
	FE_Element *elePtr;
	DOF_Group *dofPtr;

	int parameterID;
	int sensitivityFlag;
	int gradNumber;
	Vector *massMatrixMultiplicator;
	Vector *dampingMatrixMultiplicator;
	int assemblyFlag;
	Vector independentRHS;

		int currentNumGrads;
	Vector** USens; // disp Sensitivity Vector
	Vector** UdotSens; 
	Vector** UdotdotSens;
	Vector* USensWork;
	Vector* UdotSensWork;
	Vector* UdotdotSensWork;

	double	a1;
	double	a2;
	double	a3;
	double	a4;
	double	a5;
	double	a6;
	double	a7;
	double	a8;
	double dt;
	double alphaM;
	double betaK;

};

#endif

