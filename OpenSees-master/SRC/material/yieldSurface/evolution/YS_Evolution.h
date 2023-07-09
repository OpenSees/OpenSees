// YS_Evolution.h
//
//////////////////////////////////////////////////////////////////////

#if !defined YS_EVOLUTION_H
#define YS_EVOLUTION_H

//#include "YieldSurface_BC.h"

#include <Vector.h>
#include <Matrix.h>
#include <DomainComponent.h>
#include <MovableObject.h>

// min_, max_ provide functions for getting min/max of 2 quantities
#define min_(a,b)  (((a) < (b)) ? (a) : (b))
#define max_(a,b)  (((a) > (b)) ? (a) : (b))
#define sign(a)    (((a) < (0)) ?(-1) : (1))

class Information;
class Response;
class YieldSurface_BC;

class YS_Evolution : public TaggedObject, public MovableObject
{
public:
    YS_Evolution(int tag, int classTag, double iso_ratio, double kin_ratio,
				 int _dimension, double shr_iso=0.5, double shr_kin=0.5);
	virtual ~YS_Evolution();

//  Methods inherited
    virtual void	Print(OPS_Stream &s, int flag =0);
//	friend OPS_Stream &operator<<(OPS_Stream &s, const YS_Evolution &hModel);

	virtual int	update(int flag=0);
	virtual int commitState();
	virtual int	revertToLastCommit(void);
	virtual YS_Evolution *getCopy(void) = 0;

	virtual Response *setResponse(char **argv, int argc, OPS_Stream &output)=0;
	virtual int 	getResponse(int responseID, Information &info)=0;
	virtual int	 displaySelf(Renderer &theViewer, int displayMode, float fact)=0;

// Public methods called by the yield surface
	virtual int 	evolveSurface(YieldSurface_BC *ys, double magPlasticDefo, Vector &G, 
								  Vector &F_Surface, int flag=0)=0;

	// needed by ys->add Kp
	virtual const   Vector &getEquiPlasticStiffness(void)=0;
	  
    virtual double getTrialPlasticStrains(int dof)=0;
    virtual double getCommitPlasticStrains(int dof)=0;
    virtual void   setResidual(double res=1.0);
	void	setInitTranslation(Vector &initTranslate);
    const Vector &getInitTranslation(void);

	double	getCommitTranslation(int dof);
	double	getTrialTranslation(int dof);
	double	getTrialIsotropicFactor(int dof);
	double  getCommitIsotropicFactor(int dof);
	
	void    setDeformable(bool defo);

	void	toDeformedCoord(double &x);
	void	toDeformedCoord(double &x, double &y);
	void	toDeformedCoord(double &x, double &y, double &z);

	 void	toOriginalCoord(double &x);
	 void	toOriginalCoord(double &x, double &y);
	 void	toOriginalCoord(double &x, double &y, double &z);
public:
	bool  freezeEvolution;


	
private:
	void	checkDimension(int dir);
	void    toDeformedCoord(Vector &coord);
	void    toOriginalCoord(Vector &coord);
protected:
	//double  isotropicFactor_hist, isotropicFactor; // Effective magnification
	bool   deformable;
	Vector  isotropicFactor_hist, isotropicFactor;
	Vector  translate_hist,       translate;
	Vector  translate_init;
	double	isotropicRatio_orig,  isotropicRatio, isotropicRatio_shrink;
	double	kinematicRatio_orig,  kinematicRatio, kinematicRatio_shrink;
	int    dimension;
static Vector crd1, crd2, crd3;
};

#endif
