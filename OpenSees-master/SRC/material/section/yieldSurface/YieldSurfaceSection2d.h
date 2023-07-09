// @ rkaul@stanford.edu
// @ ggd@stanford.edu

#ifndef YieldSurfaceSection2d_h
#define YieldSurfaceSection2d_h

#include <SectionForceDeformation.h>
#include <YieldSurface_BC.h>
#include <Matrix.h>
#include <Vector.h>

class Channel;
class FEM_ObjectBroker;
class Information;

class YieldSurfaceSection2d: public SectionForceDeformation
{
 public:
  YieldSurfaceSection2d ( int tag, int classtag,
			  YieldSurface_BC *ptrys, 
			  bool use_kr=true);
  YieldSurfaceSection2d (void);    
  ~YieldSurfaceSection2d (void);
  
  virtual int commitState (void);
  virtual int revertToLastCommit (void);
  virtual int revertToStart (void);
  
  virtual int setTrialSectionDeformation (const Vector&);
  virtual const Vector &getSectionDeformation (void);
  
  const Vector &getStressResultant (void);
  const Matrix &getSectionTangent (void);
  const Matrix &getSectionFlexibility (void);
  
  const ID &getType (void);
  int getOrder (void) const;
  
  int sendSelf (int commitTag, Channel &theChannel);
  int recvSelf (int commitTag, Channel &theChannel,
		FEM_ObjectBroker &theBroker);
  
  virtual void Print (OPS_Stream &s, int flag =0);
  
  virtual SectionForceDeformation *getCopy (void)=0;
  
 protected:
  virtual void getSectionStiffness(Matrix &Ks)=0;
  const   bool use_Kr_orig;
  YieldSurface_BC *ys;
  Vector e; // section trial deformations
  Vector s;
  Vector eCommit;
  Vector sCommit;
  Matrix ks;
  
 private:
  //    int algo;
  bool use_Kr, split_step;
  
  static ID code;
  static Vector dele;
  static Vector surfaceForce;
  static Matrix G;
  static Matrix Ktp;
};

#endif
