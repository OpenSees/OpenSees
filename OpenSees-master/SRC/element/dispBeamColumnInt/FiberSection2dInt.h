// $Source: /usr/local/cvs/OpenSees/SRC/element/dispBeamColumnInt/FiberSection2dInt.h,v $
// $Revision: 1.3 $
// $Date: 2008-04-14 21:22:20 $

// Created: 07/04
// Modified by: LMS 
// Description: This file contains the class implementation of FiberSection2dInt.Based on FiberSection2d.cpp.

                                                                        
#ifndef FiberSection2dInt_h
#define FiberSection2dInt_h

#include <SectionForceDeformation.h>
#include <Vector.h>
#include <Matrix.h>

class UniaxialMaterial;
class Fiber;


class FiberSection2dInt : public SectionForceDeformation
{
  public:
    FiberSection2dInt(); 
    FiberSection2dInt(int tag, 
		      int numFibers, 
		      Fiber **fibers, 
		      int numHFibers, 
		      Fiber **Hfibers, 
		      int NStrip1, 
		      double tavg1, 
		      int NStrip2, 
		      double tavg2, 
		      int NStrip3, 
		      double tavg3); 
    ~FiberSection2dInt();

    int   setTrialSectionDeformation(const Vector &deforms); 
    int   setTrialSectionDeformationB(const Vector &deforms, double L); 
    const Vector &getSectionDeformation (void);
    int   revertToLastCommitB(double L);    
    int   revertToStartB(void);
    int   commitStateB(void);
    void   beta(double e0, double e2, double &sc1, double &tc1, double &tc12, double &beta);



    const Vector &getSigmaY(void);
    const Vector &getTau(void);
    const Vector &getAlpha(void);
    const Vector &getIter(void);
    const Vector &getEX(void);
    const Vector &getEY(void);
    const Vector &getE1(void);
    const Vector &getE2(void);
    const Vector &getSX(void);
    const Vector &getSY(void);
    const Vector &getS1(void);
    const Vector &getS2(void);

    const Vector &getStressResultant(void);
    const Matrix &getSectionTangent(void);
    const Matrix &getInitialTangent(void);

    int   commitState(void);
    int   revertToLastCommit(void);    
    int   revertToStart(void);

    SectionForceDeformation *getCopy(void);
    const ID &getType (void);
    int getOrder (void) const;

    int sendSelf(int cTag, Channel &theChannel);
    int recvSelf(int cTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);
    void Print(OPS_Stream &s, int flag = 0);

    Response *setResponse(const char **argv, int argc, OPS_Stream &s);
    int getResponse(int responseID, Information &info);
    int addFiber(Fiber &theFiber);

    // AddingSensitivity:BEGIN //////////////////////////////////////////
    const Vector & getStressResultantSensitivity(int gradNumber, bool conditional);
    const Vector & getSectionDeformationSensitivity(int gradNumber);
    const Matrix & getSectionTangentSensitivity(int gradNumber);
    int   commitSensitivity(const Vector& sectionDeformationGradient, int gradNumber, int numGrads);
    // AddingSensitivity:END ///////////////////////////////////////////

  protected:
    
  private:
    int numFibers;                   // number of fibers in the section
    UniaxialMaterial **theMaterials1; // array of pointers to materials
    UniaxialMaterial **theMaterials2; // array of pointers to materials
    
    double   *matData;               // data for the materials [yloc and area]
    
    int numHFibers;                   // number of fibers in the section
    UniaxialMaterial **theHMaterials; // array of pointers to materials
    double   *matHData;               // data for the materials [yloc and area]
    
    double   kData[9];               // data for ks matrix 
    double   sData[3];               // data for s vector 
    
    static ID code;

    int NStrip;    
    int NStrip1;
    double tavg1;
    int NStrip2;
    double tavg2;
    int NStrip3;
    double tavg3;

    
    double   sy[100];       //100 = max num of strips        
    double   txy[100];     
    double   alfa[100];
    double   iterOut[100];
    double   iterCommit[100];
    double   exOut[100];
    double   exCommit[100];
    double   eyCommit[100];
    double   e1Commit[100];
    double   e2Commit[100];

    double   alfaCommit[100];
    double   sxCommit[100];
    double   syCommit[100];
    double   s1Commit[100];
    double   s2Commit[100];

    Vector StripCenterLoc;	
    Matrix StripLoc;		
    Vector FiberLoc;		

    double yBar;       // Section centroid
    double ymax;       // Section centroid
    double ymin;       // Section centroid
    Vector e;          // trial section deformations 
    Vector eCommit;    // committed section deformations 
    Vector *s;         // section resisting forces  (axial force, bending moment)
    Matrix *ks;        // section stiffness    
    
    Vector *sigmaY;
    Vector *tau;
    Vector *alpha;
    Vector *alphaCommit;
    Vector *iterFile;
    Vector *exf;
    Vector *e1f;
    Vector *e2f;
    Vector *eyf;
    Vector *sxf;
    Vector *s1f;
    Vector *s2f;
    Vector *syf;
    
    
    // AddingSensitivity:BEGIN //////////////////////////////////////////
    int parameterID;
    // AddingSensitivity:END ///////////////////////////////////////////
};

#endif

