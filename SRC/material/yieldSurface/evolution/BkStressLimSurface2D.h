// BkStressLimSurface2D.h
//
//////////////////////////////////////////////////////////////////////

#if !defined BKSTRESSLIMSURFACE2D_H
#define BKSTRESSLIMSURFACE2D_H

#include "YS_Evolution2D.h"
class YieldSurface_BC;
  
class BkStressLimSurface2D : public YS_Evolution2D
{
public:
    BkStressLimSurface2D(int tag, int classTag, double min_iso_factor, 
				double iso_ratio, double kin_ratio,
				YieldSurface_BC  &lim_surface,
				PlasticHardeningMaterial &kinX,
				PlasticHardeningMaterial &kinY, 
				PlasticHardeningMaterial &isoXPos,
				PlasticHardeningMaterial &isoXNeg,
				PlasticHardeningMaterial &isoYPos,
				PlasticHardeningMaterial &isoYNeg,
				int restype, double res_Fact, double app_Fact, double dir
				);

	virtual ~BkStressLimSurface2D();
	virtual int 	commitState();
	virtual int		revertToLastCommit(void);

	virtual int	 displaySelf(Renderer &theViewer, int displayMode, float fact);
	virtual void	Print(OPS_Stream &s, int flag =0);
	        void    setResidual(double res=1.0);
	virtual YS_Evolution *getCopy(void)=0;
	virtual const   Vector &getEquiPlasticStiffness(void);
	double getTrialPlasticStrains(int dof);
	double getCommitPlasticStrains(int dof);

protected:
	virtual void	setTrialPlasticStrains(double ep, const Vector &f, const Vector &g);
	virtual double getIsoPlasticStiffness(int dir);
	virtual double getKinPlasticStiffness(int dir);
	virtual Vector& getEvolDirection(Vector &f_new);
	
protected:
//	double minIsoFactor;
	PlasticHardeningMaterial  *kinMatX, *kinMatY;
	PlasticHardeningMaterial  *isoMatXPos, *isoMatXNeg;
	PlasticHardeningMaterial  *isoMatYPos, *isoMatYNeg;
	
	YieldSurface_BC   *limSurface;
	bool defPosX, defPosY;
	bool resHardening, resApproach;
	int resAlgo;
	double resFactor, appFactor;
	double direction, direction_orig;
};

#endif
