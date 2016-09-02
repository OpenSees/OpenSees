// HingeForceDeformation.h: interface for the HingeForceDeformation class.
//
//////////////////////////////////////////////////////////////////////

#if !defined YIELDSURFACE_BC_H
#define YIELDSURFACE_BC_H

#include <Renderer.h>
#include <Vector.h>
#include <Matrix.h>
#include <ID.h>
#include "YS_Evolution.h"

#include <DomainComponent.h>
#include <MovableObject.h>

// min_, max_ provide functions for getting min/max of 2 quantities
#define min_(a,b)  (((a) < (b)) ? (a) : (b))
#define max_(a,b)  (((a) > (b)) ? (a) : (b))

class Information;
class Response;

class YieldSurface_BC : public TaggedObject, public MovableObject
{
public:
    YieldSurface_BC(int tag, int classtag, YS_Evolution &model,
					double capx);
	
	YieldSurface_BC(int tag, int classtag, YS_Evolution &model,
					double capx, double capy);
	
	YieldSurface_BC(int tag, int classtag, YS_Evolution &model,
					double capx, double capy, double capz);
			
	virtual ~YieldSurface_BC();

	virtual void	Print(OPS_Stream &s, int flag =0);
	 void    setEleInfo(int eleTag, int loc);
	 // keep transformation virtual
     virtual void    setTransformation(int xDof, int xFact);
     virtual void    setTransformation(int xDof, int yDof,
                                      int xFact, int yFact);
     virtual void    setTransformation(int xDof, int yDof, int zDof,
                                      int xFact, int yFact, int zFact);

    // in element system
	virtual void	getCommitGradient(Matrix &G) = 0;
	virtual void	getTrialGradient(Matrix &G, Vector &force) = 0;
    virtual double  getTrialDrift(Vector &force) = 0;

//    virtual void  setExtent()=0;
	virtual const   Vector &getExtent(void)=0;
	virtual int		update(int flag = 0);

    // in ys system
    double 	getCap(int dir);
	virtual Vector&	translationTo(Vector &f_new, Vector &f_dir)=0;
	virtual int getState(int stateInfo)=0;

    virtual double getDrift(double x1);
    virtual double getDrift(double x1, double y1);
    virtual double getDrift(double x1, double y1, double z1);

	// needed by evlution model
    virtual double interpolate(double x1, double x2);
    virtual double interpolate(double x1, double y1, double x2, double y2);
    virtual double interpolate(double x1, double y1, double z1, double x2, double y2, double z2);

    virtual int		getTrialForceLocation(Vector &force)=0;
	virtual int		getCommitForceLocation()=0;

    virtual void    addPlasticStiffness(Matrix &K)=0;

	virtual double	setToSurface(Vector &force, int algoType, int flag=0)=0;
	virtual int 	    modifySurface(double magPlasticDefo, Vector &Fsurface, Matrix &G, int flag=0)=0;

	virtual int      commitState(Vector &force);
	virtual int		revertToLastCommit(void)=0;

	virtual YieldSurface_BC *getCopy(void) = 0;
	virtual int	 displaySelf(Renderer &theViewer, int displayMode, float fact);
	virtual void	setView(Renderer *theRenderer);

	// required for debugging only
	virtual int	displayForcePoint(Vector &force, int color = 4);

protected:
    void toLocalSystem  (Vector &eleVector, double &x, bool nonDimensionalize, bool signMult=true);
    void toLocalSystem  (Vector &eleVector, double &x, double &y, bool nonDimensionalize, bool signMult=true);
    void toLocalSystem  (Vector &eleVector, double &x, double &y, double &z, bool nonDimensionalize, bool signMult=true);

	// matrix do not multiply G!
    void toLocalSystem  (Matrix &eleMatrix, double &x, bool nonDimensionalize, bool signMult=true);
    void toLocalSystem  (Matrix &eleMatrix, double &x, double &y, bool nonDimensionalize, bool signMult=true);
    void toLocalSystem  (Matrix &eleMatrix, double &x, double &y, double &z, bool nonDimensionalize, bool signMult=true);

    void toElementSystem(Vector &eleVector, double &x, bool dimensionalize, bool signMult=true);
    void toElementSystem(Vector &eleVector, double &x, double &y, bool dimensionalize, bool signMult=true);
    void toElementSystem(Vector &eleVector, double &x, double &y, double &z, bool dimensionalize, bool signMult=true);

    void toElementSystem(Matrix &eleMatrix, double &x, bool dimensionalize, bool signMult=true);
    void toElementSystem(Matrix &eleMatrix, double &x, double &y, bool dimensionalize, bool signMult=true);
    void toElementSystem(Matrix &eleMatrix, double &x, double &y, double &z, bool dimensionalize, bool signMult=true);

private:
    int checkT(void);

public:
YS_Evolution *hModel;

protected:
	Renderer *theView;
    ID       *T;
    ID       *S;
    double   capX_orig, capY_orig, capZ_orig;
	double   capX, capY, capZ;
	double   capXdim, capYdim, capZdim;
	int      dimension;
	bool     isLoading;

public:
int      ele_Tag, ele_Location;
const  static int dFReturn, RadialReturn, ConstantXReturn, ConstantYReturn;
const  static int NoFP, SurfOnly, StateLoading;
};

/*
enum Fstate { orig, trans };

class ysfp
{
public:
	ysfp(double x, enum Fstate);
	ysfp(double x, double y, enum Fstate);
	ysfp(double x, double y, double z, enum Fstate);

	double getOrig(int dof);
	double getTrans(int dof);
	//setState..
private:
	bool orig, trans;

	Vector *F0;
	Vector *Ft;
};
*/

bool OPS_addYieldSurface_BC(YieldSurface_BC *newComponent);
YieldSurface_BC *OPS_getYieldSurface_BC(int tag);
void OPS_clearAllYieldSurface_BC(void);

#endif
