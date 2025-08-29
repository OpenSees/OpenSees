//===----------------------------------------------------------------------===//
//
//                                   xara
//                              https://xara.so
//
//===----------------------------------------------------------------------===//
//
// Copyright (c) 2025, Claudio M. Perez
// All rights reserved.  No warranty, explicit or implicit, is provided.
//
// This source code is licensed under the BSD 2-Clause License.
// See LICENSE file or https://opensource.org/licenses/BSD-2-Clause
//
//===----------------------------------------------------------------------===//
//
// BasicAnalysisBuilder is an aggregate class which manages the analysis 
// objects:
//
// - LinearSOE
// - Domain                    *theDomain;
// - ConstraintHandler 	       *theHandler;
// - DOF_Numberer 	           *theNumberer;
// - AnalysisModel 	           *theAnalysisModel;
// - EquiSolnAlgo 	           *theAlgorithm;
// - EigenSOE 		             *theEigenSOE;
// - StaticIntegrator          *theStaticIntegrator;
// - TransientIntegrator       *theTransientIntegrator;
// - ConvergenceTest           *theTest;
//
// The BasicAnalysisBuilder assumes responsibility for
// deleting these objects, but ownership of the SOE may be
// given up.
//
// Written: cmp
//
#ifndef BasicAnalysisBulider_h
#define BasicAnalysisBulider_h

class Domain;
class BasicModelBuilder;
class G3_Table;
class ConstraintHandler;
class DOF_Numberer;
class AnalysisModel;
class EquiSolnAlgo;
class LinearSOE;
class EigenSOE;
class StaticIntegrator;
class TransientIntegrator;
class ConvergenceTest;

class BasicAnalysisBuilder
{
public:
    BasicAnalysisBuilder(BasicModelBuilder&);
    ~BasicAnalysisBuilder();

    enum CurrentAnalysis {
      EMPTY_ANALYSIS,
      STATIC_ANALYSIS, 
      TRANSIENT_ANALYSIS
    };

    enum Perform : int {
      Increment = 1<<0,
      Iterate   = 1<<1,
      Commit    = 1<<2
    };

    void set(ConstraintHandler* obj);
    void set(DOF_Numberer* obj);
    void set(EquiSolnAlgo* obj);
    void set(LinearSOE*  obj, bool free=true);
    void set(StaticIntegrator& obj);
    void set(TransientIntegrator& obj, bool free=true);
    void set(ConvergenceTest* obj);
    void set(EigenSOE& obj);

    LinearSOE* getLinearSOE();

    Domain* getDomain();
    const BasicModelBuilder& getContext() const { return context; }

    int initialize();

    int  newTransientAnalysis();
    int  setStaticAnalysis();
    int  setTransientAnalysis();

    //   Eigen
    void newEigenAnalysis(int typeSolver, double shift);
    int  eigen(int numMode, bool generalized, bool findSmallest);
    int  getNumEigen() {return numEigen;};

    int formUnbalance();

    const EquiSolnAlgo*  getAlgorithm() const;
    StaticIntegrator*    getStaticIntegrator();
    TransientIntegrator* getTransientIntegrator();

    // for getCTestIter command
    ConvergenceTest*     getConvergenceTest();

    int domainChanged();

    // Performing analysis
    int analyze(int num_steps, double size_steps, int flag=Increment|Iterate|Commit);
    int analyzeStatic(int num_steps, int flag);
    
    int analyzeTransient(int numSteps, double dT);
    int analyzeVariable(int numSteps, double dT, double dtMin, double dtMax, int Jd);
private:
    int analyzeStep(double dT);
    int analyzeSubLevel(int level, double dT);

public:
    int analyzeGradient();
    int setGradientType(int flag);

    void wipe();

    
    enum CurrentAnalysis  CurrentAnalysisFlag = EMPTY_ANALYSIS;

private:
    void setLinks(CurrentAnalysis flag = EMPTY_ANALYSIS);
    void fillDefaults(enum CurrentAnalysis flag);

    BasicModelBuilder&         context;
    Domain                    *theDomain;
    ConstraintHandler         *theHandler;
    DOF_Numberer              *theNumberer;
    AnalysisModel             *theAnalysisModel;
    EquiSolnAlgo              *theAlgorithm;
    LinearSOE                 *theSOE;
    EigenSOE                  *theEigenSOE;
    StaticIntegrator          *theStaticIntegrator;
    TransientIntegrator       *theTransientIntegrator;
    ConvergenceTest           *theTest;

    int domainStamp;
    int numEigen = 0;

    int numSubLevels = 0;
    int numSubSteps  = 0;

    bool freeSOE = true;
    bool freeTI  = true;

};

#endif
