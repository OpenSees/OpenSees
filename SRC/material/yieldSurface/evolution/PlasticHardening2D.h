// PlasticHardening2D.h
//
//////////////////////////////////////////////////////////////////////

#if !defined PLASTICHARDENING2D_H
#define PLASTICHARDENING2D_H

#include "YS_Evolution2D.h"

class PlasticHardening2D : public YS_Evolution2D
{
public:
    PlasticHardening2D(int tag, int classTag, double min_iso_factor, 
				double iso_ratio, double kin_ratio,
				 PlasticHardeningMaterial &kpx_pos,
				 PlasticHardeningMaterial &kpx_neg,
				 PlasticHardeningMaterial &kpy_pos,
				 PlasticHardeningMaterial &kpy_neg, double dir);

	virtual ~PlasticHardening2D();
	
	virtual int 	commitState();
	virtual int		revertToLastCommit(void);

	virtual void	Print(OPS_Stream &s, int flag =0);
	virtual YS_Evolution *getCopy(void)=0;
	virtual const Vector &getEquiPlasticStiffness(void);
	double getTrialPlasticStrains(int dof);
	double getCommitPlasticStrains(int dof);
	
protected:
	virtual void	setTrialPlasticStrains(double ep, const Vector &f, const Vector &g);
	virtual double	getIsoPlasticStiffness(int dir);
	virtual double	getKinPlasticStiffness(int dir);
	virtual Vector& getEvolDirection(Vector &f_new);

protected:
//	double minIsoFactor;
	PlasticHardeningMaterial  *kpMatXPos, *kpMatYPos;
	PlasticHardeningMaterial  *kpMatXNeg, *kpMatYNeg;
	bool defPosX, defPosY;
	double direction;
//	bool pinchX; 	
};

#endif
