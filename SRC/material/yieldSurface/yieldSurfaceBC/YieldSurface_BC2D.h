// YieldSurfaceBC_2D.h: interface for the HingeForceDeformation class.
//
//////////////////////////////////////////////////////////////////////

#if !defined YIELDSURFACE_BC2D_H
#define YIELDSURFACE_BC2D_H

#include "YieldSurface_BC.h"
#include <UniaxialMaterial.h>

class YieldSurface_BC2D : public YieldSurface_BC
{
public:
    YieldSurface_BC2D(int tag, int classTag, double xmax, double ymax,
                     YS_Evolution &model);

	virtual ~YieldSurface_BC2D();


//    virtual Response *setResponse(char **argv, int argc, Information &info);
//    virtual int 	getResponse(int responseID, Information &info);
    virtual void	Print(OPS_Stream &s, int flag =0)=0;
    virtual int 	sendSelf(int commitTag, Channel &theChannel){return -1;}
    virtual int 	recvSelf(int commitTag, Channel &theChannel,
					FEM_ObjectBroker &theBroker){return -1;}

    virtual void	setTransformation(int xDof, int yDof, int xFact, int yFact);

	virtual void	getCommitGradient(Matrix &G);
	virtual void	getTrialGradient(Matrix &G, Vector &force);
	//virtual const   Vector &getGradient(void);
	//virtual const   Vector &getTrialGradient(void);

	virtual int		update(int flag = 0);
    virtual int 	getState(int stateInfo);
    virtual double  getTrialDrift(Vector &force);
    virtual int		getTrialForceLocation(Vector &force);
    virtual int		getCommitForceLocation();
	//virtual const   Vector &getForce(void);
	//virtual const   Vector &getTrialForce(void);

	//virtual int     setTrialForce(Vector &force);

	// double  getIsotropicFactor(void){ return hModel->getIsotropicFactor();}

    virtual void    addPlasticStiffness(Matrix &K);
//    virtual void	checkState(Vector &trialforce, bool &plastify, bool &shootsthrough);

	virtual double	setToSurface(Vector &force, int algoType, int colorFlag = 0);
	virtual int 	modifySurface(double magPlasticDefo, Vector &Fsurface, Matrix &G, int flag=0);
	//virtual int     trialModifySurface(double magPlasticDefo);
	//virtual double	getElasticForce(Vector &force, Vector &elasticForce);

	virtual int 	commitState(Vector &force);
	virtual int		revertToLastCommit(void);

	virtual YieldSurface_BC *getCopy(void) = 0;
	virtual int		displaySelf(Renderer &theViewer, int displayMode, float fact);

			int		displayCommitForcePoint(Renderer &theViewer, int displayMode, float fact);
	        int 	displayForcePoint(bool toDeformed, double f_x, double f_y, int color);
	        int		displayForcePoint(Vector &force, int color = 4);
//protected:
	virtual Vector&	translationTo(Vector &f_new, Vector &f_dir);
	virtual double	getDrift(double x, double y);
//  For the following 2 methods, x, y already non-dimensionalized
    virtual void 	getGradient(double &gx, double &gy, double x, double y)=0;
    virtual double	getSurfaceDrift(double x, double y)=0;
    virtual void	setExtent()=0;
	virtual const   Vector &getExtent(void);

    virtual int		forceLocation(double drift);
	virtual double	interpolate(double xi, double yi, double xj, double yj);
	virtual void	customizeInterpolate(double &xi, double &yi, double &xj, double &yj);

	double	interpolateClose(double xi, double yi, double xj, double yj);
//  Dimensionalizing taken care at Element System <--> Local System level
//    		void 	toDeformedCoord(double &x, double &y);
//    		void 	toOriginalCoord(double &x, double &y);

protected:
//	UniaxialMaterial  *kpMatX, *kpMatY;

    double xPos, xNeg, yPos, yNeg;          // Extent along X and Y
    double a1, b1, a2, b2, a3, b3, a4, b4;  // y = ax +b -> inner quad
    double offset, increment;
//    double isotropicRatio;          // Ratio of the plastic deformation that
                                    // is isotropic
/////////////////////////////////////////////////////////////////////////////
//  state variables
/////////////////////////////////////////////////////////////////////////////
//    double translateX, translateY;  // Kinematic displacement
//    double isotropicFactor;         // Magnification factor
//    double sumPlasticDeformX, sumPlasticDeformY;

//    double translateX_hist, translateY_hist;
    int    status_hist;
    int state;
//    double isotropicFactor_hist;
//    double sumPlasticDeformX_hist, sumPlasticDeformY_hist;
    double fx_hist, fy_hist, gx_hist, gy_hist;
    double fx_trial, fy_trial, gx_trial, gy_trial;

    static Vector v6;
	static double  error;
	static Vector v2;
	static Vector g2;
	static Vector v4;
	static Vector T2;
	static Vector F2;
public:
//	const  static int dFReturn, RadialReturn, ConstantXReturn, ConstantYReturn;

};

#endif
