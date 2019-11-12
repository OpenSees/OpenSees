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
                                                                        
//----------------------------------------------------------------------------------------------------------------------------
 // Developed by:
 // Adam M. Knaack (adam.knaack@schaefer-inc.com) 
 // Schaefer-Inc, Cincinnati, Ohio, USA
 // Nikola D. Tosic (ntosic@imk.grf.bg.ac.rs)
 // Department for Materials and Structure, Faculty of Civil Engineering, University of Belgrade, Serbia
 // Yahya C. Kurama (ykurama@nd.edu)
 // Department of Civil and Environmental Engineering and Earth Sciences, College of Engineering, University of Notre Dame, Notre Dame, Indiana, USA
 //----------------------------------------------------------------------------------------------------------------------------

 //----------------------------------------------------------------------------------------------------------------------------
 // Created: 2012
 // Last updated: 2019
 //----------------------------------------------------------------------------------------------------------------------------

 //----------------------------------------------------------------------------------------------------------------------------
 // Description: This file contains the source code of CreepAnalysis. 
 // CreepAnalysis is an analysis procedure that calculates time-dependent
 // concrete creep and shrinkage strains.
 //---------------------------------------------------------------------------------------------------------------------------- 
 // Detailed descriptions of the model and its implementation can be found in the following:
 // (1) Knaack, A.M., Kurama, Y.C. 2018. Modeling Time-Dependent Deformations: Application for Reinforced Concrete Beams with 
 //     Recycled Concrete Aggregates. ACI Structural J. 115, 175–190. doi:10.14359/51701153
 // (2) Knaack, A.M., 2013. Sustainable concrete structures using recycled concrete aggregate: short-term and long-term behavior
 //     considering material variability. PhD Dissertation, Civil and Environmental Engineering and Earth Sciences, University of Notre Dame, Notre Dame, Indiana, USA, 680 pp.
 // A manual describing the use of the model and sample files can be found at:
 // http://dx.doi.org/10.17632/z4gxnhchky.2
 //----------------------------------------------------------------------------------------------------------------------------

 //----------------------------------------------------------------------------------------------------------------------------
 // Disclaimer: This software is provided “as is”, without any warranties, expressed or implied. In no event shall the developers be liable for any claim, damages, or liability arising from or in connection with this software.
 //----------------------------------------------------------------------------------------------------------------------------
                                                                        
                                                                        
#ifndef CreepAnalysis_h
#define CreepAnalysis_h

#include <Analysis.h>


// AddingSensitivity:BEGIN //////////////////////////////////
#ifdef _RELIABILITY
#include <SensitivityAlgorithm.h>
#endif
// AddingSensitivity:END ////////////////////////////////////


class ConstraintHandler;
class DOF_Numberer;
class AnalysisModel;
class StaticIntegrator;
class LinearSOE;
class EquiSolnAlgo;
class ConvergenceTest;
class EigenSOE;

class CreepAnalysis: public Analysis
{
  public:
    CreepAnalysis(Domain &theDomain,
		   ConstraintHandler &theHandler,
		   DOF_Numberer &theNumberer,
		   AnalysisModel &theModel,
		   EquiSolnAlgo &theSolnAlgo,		   
		   LinearSOE &theSOE,
		   StaticIntegrator &theIntegrator,
		   ConvergenceTest *theTest = 0);

    ~CreepAnalysis();

    void clearAll(void);	    
    
    int analyze(int numSteps);
    int eigen(int numMode, bool generalized, bool findSmallest);
    int initialize(void);
    int domainChanged(void);

    int setNumberer(DOF_Numberer &theNumberer);
    int setAlgorithm(EquiSolnAlgo &theAlgorithm);
    int setIntegrator(StaticIntegrator &theIntegrator);
    int setLinearSOE(LinearSOE &theSOE);
    int setConvergenceTest(ConvergenceTest &theTest);
    int setEigenSOE(EigenSOE &theSOE);

    EquiSolnAlgo     *getAlgorithm(void);
    StaticIntegrator *getIntegrator(void);
    ConvergenceTest  *getConvergenceTest(void);

    // AddingSensitivity:BEGIN ///////////////////////////////
#ifdef _RELIABILITY
    int setSensitivityAlgorithm(SensitivityAlgorithm *theSensitivityAlgorithm);
#endif
    // AddingSensitivity:END /////////////////////////////////
    
  protected: 
    
  private:
    ConstraintHandler 	*theConstraintHandler;    
    DOF_Numberer 	*theDOF_Numberer;
    AnalysisModel 	*theAnalysisModel;
    EquiSolnAlgo 	*theAlgorithm;
    LinearSOE 		*theSOE;
    EigenSOE 		*theEigenSOE;
    StaticIntegrator    *theIntegrator;
    ConvergenceTest     *theTest;

    int domainStamp;

    // AddingSensitivity:BEGIN ///////////////////////////////
#ifdef _RELIABILITY
    SensitivityAlgorithm *theSensitivityAlgorithm;
#endif
    // AddingSensitivity:END ///////////////////////////////
};

#endif
