#ifndef INELASTICYS2DGNL_H
#define INELASTICYS2DGNL_H

#include <UpdatedLagrangianBeam2D.h>
#include <YieldSurface_BC.h>
#include <Vector.h>

#define DISPLAY_YS 2745

/**Inelastic Element - concentrated hinge model, Fi - Fj interaction
   at each ends using yield surfaces
  *@author rkaul
  */

class InelasticYS2DGNL : public UpdatedLagrangianBeam2D
{
 public:
  InelasticYS2DGNL(int tag, 
		   int Nd1, int Nd2,
		   YieldSurface_BC *ysEnd1, YieldSurface_BC *ysEnd2,
		   int rf_algo = -1, // updated
		   bool islinear = false, double rho = 0.0);
  
  ~InelasticYS2DGNL();
  
  virtual const 	Vector &getResistingForce(void);
  virtual const	Matrix &getTangentStiff(void);
  virtual int		commitState(void);
  virtual int		update(void);
  virtual int displaySelf(Renderer &theViewer, int displayMode, float fact);
  void Print(OPS_Stream &s, int flag =0);
  int sendSelf(int commitTag, Channel &theChannel);
  int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
  //    void createView(char *title, WindowManager *theWM, char displaytype);

  void createView(char *title, double scale, int x, int y, int cx, int cy, char displaytype = 'l');

	virtual Response *setResponse(const char **argv, int argc);
	virtual int getResponse(int responseID, Information &eleInformation);

 protected:
  virtual void getLocalStiff(Matrix &K) = 0;
  virtual void getLocalMass(Matrix &M);
  
protected:
  int  computeTrueEleForce(Vector &trial_force);
  void checkSpecialCases(void);
 
private: 
  void forceBalance(Vector &force, int algo);
  void plastifyOneEnd(int end, YieldSurface_BC *ys,  Vector &trial_force,
		      Vector &incrDisp, Matrix &K, Vector &total_force, int algo);
  
  void splitStep(int end_shoot, YieldSurface_BC *ys_shoots, YieldSurface_BC *ys_drifts,
		 Vector &trial_force, Matrix &K, Vector &total_force);
  
  void driftOneEnd(YieldSurface_BC *ys, Vector &trial_force, Vector &surface_force,
		   Matrix &K, Vector &total_force);
  
  void driftBothEnds(Vector &trial_force, Vector &surface_force,
		     Matrix &K, Vector &total_force);
  
  void plastifyBothEnds(Vector &trial_force, Vector &incrDisp,
			Matrix &K, Vector &total_force);
  void checkEndStatus(bool &end1drifts, bool &end2drifts, Vector &trialForce);
  int  plasticPredictor(Vector &trialForce);
  int  elasticCorrector(Vector &trialForce, int algo);
  
 protected:
  static Vector elasticForce;
  static Vector F1, F2, Fs;
  

  YieldSurface_BC *ys1, *ys2;
  char displayType;
  Renderer *pView;
  ColorMap *theMap;

  bool end1Plastify, end2Plastify;
  bool end1Plastify_hist, end2Plastify_hist;

  Matrix end1G, end2G;
  Matrix Stiff;
  int  forceRecoveryAlgo, forceRecoveryAlgo_orig;
  bool end1Damage, end2Damage;
  bool split_step;

  int debug, fdebug, pdebug, ydebug, statusDebug;
  bool init;
  bool updateKt;


  
  const static int INSIDE, OUTSIDE, WITHIN;
  
  static double storage;
};

#endif
