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
                                                                        
// $Revision: 1.0 $
// $Date: 2014-01-16 10:03:47 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/FEsensitivity/PFEMSensitivityIntegrator.h,v $


//
// Written by Minjie Zhu (Oregon State University)
//

#ifndef PFEMSensitivityIntegrator_h
#define PFEMSensitivityIntegrator_h

#include <PFEMIntegrator.h>
#include <SensitivityIntegrator.h>
class FE_Element;
class DOF_Group;
class Vector;
class Information;


class PFEMSensitivityIntegrator : public SensitivityIntegrator , public PFEMIntegrator
{
public:
    PFEMSensitivityIntegrator();
    explicit PFEMSensitivityIntegrator(int assemblyFlag);
    ~PFEMSensitivityIntegrator();
    
    int setParameter      (char **argv, int argc, Information &info);
    int updateParameter   (int parameterID, Information &info);
    int activateParameter (int parameterID);

    int formEleResidual(FE_Element *theEle);
    int formNodUnbalance(DOF_Group *theDof);

    int formSensitivityRHS(int gradNum);
    int formIndependentSensitivityRHS();
    int saveSensitivity   (const Vector &v, int gradNum, int numGrads);
    int commitSensitivity (int gradNum, int numGrads);  
    /////S added by K Fujimura /////
    int updateGradNumber(int passedGradNumber);
    int sensitivityDomainChanged(int NumGrads);
    bool staticSensitivity(void);
    bool NewSensitivity(void);
    /////E added by K Fujimura /////

protected:
    
private:

    int parameterID;
    int sensitivityFlag;
    int gradNumber;
    int assemblyFlag;
    Vector independentRHS;
    Vector dVn;
};

#endif

