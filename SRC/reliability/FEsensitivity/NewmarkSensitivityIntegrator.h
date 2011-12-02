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
// $Date: 2003-03-04 00:46:01 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/FEsensitivity/NewmarkSensitivityIntegrator.h,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#ifndef NewmarkSensitivityIntegrator_h
#define NewmarkSensitivityIntegrator_h

#include <Newmark.h>
#include <SensitivityIntegrator.h>
class FE_Element;
class DOF_Group;
class Vector;
class Information;


class NewmarkSensitivityIntegrator : public SensitivityIntegrator , public Newmark
{
  public:
    NewmarkSensitivityIntegrator();
    NewmarkSensitivityIntegrator(int assemblyFlag, double gamma, double beta, bool disp = true);
    NewmarkSensitivityIntegrator(int assemblyFlag, double gamma, double beta, double alphaM, double betaKcurrent,
	    double betaKinit, double betaKlastCommit, bool disp = true);
    ~NewmarkSensitivityIntegrator();
    
    int setParameter      (char **argv, int argc, Information &info);
    int updateParameter   (int parameterID, Information &info);
	int activateParameter (int parameterID);

	int formEleResidual(FE_Element *theEle);
	int formNodUnbalance(DOF_Group *theDof);

	int formSensitivityRHS(int gradNum);
	int formIndependentSensitivityRHS();
	int saveSensitivity   (const Vector &v, int gradNum, int numGrads);
	int commitSensitivity (int gradNum, int numGrads);  

  protected:
    
  private:

	int parameterID;
	int sensitivityFlag;
	int gradNumber;
	Vector *massMatrixMultiplicator;
	Vector *dampingMatrixMultiplicator;
	int assemblyFlag;
	Vector independentRHS;

};

#endif

