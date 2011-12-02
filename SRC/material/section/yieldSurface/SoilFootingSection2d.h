// Date: 12/08/2004
//
// File name: SoilFootingSection2d.h
// 
// Coded by: Sivapalan Gajan <sgajan@ucdavis.edu>
//
// Description: This file contains the class definition for SoilFootingSection2d.
// SoilFootingSection2d is a footing-soil interface contact material, a subclass 
// of SectionForceDeformation, and should be used with a ZeroLengthSection 
// element.



#ifndef SoilFootingSection2d_h
#define SoilFootingSection2d_h

#include <SectionForceDeformation.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>

#include <fstream>
#include <iostream>

using namespace std;

class Channel;
class FEM_ObjectBroker;
class Information;

class SoilFootingSection2d: public SectionForceDeformation
{
   public:

      // constructors and destructor
      SoilFootingSection2d (int tag, double fs, double vult, double l, double kv, 
                            double kh, double rv, double deltaL);
      SoilFootingSection2d (void);    
      ~SoilFootingSection2d (void);

      // commit methods
      virtual int commitState (void);
      virtual int revertToLastCommit (void);
      virtual int revertToStart (void);
  
      // methods related to section deofrmation
      virtual int setTrialSectionDeformation (const Vector&);
      virtual const Vector &getSectionDeformation (void);
  
      // methods related to section forces/stiffness
      const Vector &getStressResultant (void);
      const Matrix &getSectionTangent (void);
      const Matrix &getSectionFlexibility (void);
  
      // other get methods
      const ID &getType (void);
      int getOrder (void) const;
      virtual SectionForceDeformation *getCopy (void);
  
      // methods used by database/parallel programming
      int sendSelf (int commitTag, Channel &theChannel);
      int recvSelf (int commitTag, Channel &theChannel,
		    FEM_ObjectBroker &theBroker);
  
      virtual void Print (OPS_Stream &s, int flag =0);

      const Matrix& getInitialTangent();
  
//      virtual void getSectionStiffness (Matrix &Ks);
  
   protected:

      // used only by this class
      void initializeInternalVariables (void);
      void initializeBoundingSurface (void);
      void tempFunction (void); 
      int applyLoading (Vector);  



   private:


   
      // basic model parameters
      double FS;
      double Vult, qult;
      double V;
      double L;
      double Kv, Kh, Kt;
      double M;

      double h, hCurr, hPrev, dhCurr, dVtemp;

      // info about section force, disp, and stiffness
      Vector e, eCommit, deModel;
      Vector s, sCommit;
      Matrix ks, ksE;
      double dTheta, dThetaPrev;
      double HPrevCommit, MPrevCommit;


      // other members
      double Rv;
      double dL;
      double ecc, eccCommit;
      double Mult, Melastic, Mmaxpast;

      int noNodes;
      int incr;
      int ini_size;
      int c1, c1T, c2, c2T;
      int c1Commit, c1TCommit, c2Commit, c2TCommit;

      // internal variables
      double **foot;
      double **soilMin, **soilMax;
      double **pressure, **pressMax;
      double thetaPlus, thetaMinus, thetaRange, Mlimit;
      double thetaPlusPrev, thetaMinusPrev;
 
      // bounding surface parameters
      double a, b, ccc, d, eee, f;
      double A, B;
      double Fv, Fh, Fm;
      double alpha, beta, pult;
      double tolerance;
      double soilFree;

      static ID code;
      int isOver, isdV;
      int isElastic;       
      double dTh, dThP;

};

#endif
