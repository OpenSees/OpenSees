#ifndef PlasticHardeningMaterial_h
#define PlasticHardeningMaterial_h

#include <Material.h>
class Information;
class Response;

class PlasticHardeningMaterial : public Material
{
  public:
    PlasticHardeningMaterial (int tag, int classTag);
    virtual ~PlasticHardeningMaterial();

	// could be delta_plastic, |back-stress|, 
	// distance between force-point on YS and conjugate
	// point on Bounding-Surface
    virtual int  setTrialValue(double xVal, double factor);
	virtual int  setTrialIncrValue(double dxVal);
	        void setResidual(double res=1);

    virtual int  commitState ();
    virtual int  revertToLastCommit (void);
	virtual int  revertToStart (void);

	virtual double getTrialPlasticStiffness()=0;
			 double getTrialValue(void);
    virtual PlasticHardeningMaterial *getCopy (void) = 0;

    virtual Response *setResponse (char **argv, int argc, OPS_Stream &output);
    virtual int getResponse (int responseID, Information &matInformation);
	virtual void Print(OPS_Stream &s, int flag =0);

	virtual int  sendSelf(int commitTag, Channel &theChannel){return -1;}
    virtual int  recvSelf(int commitTag, Channel &theChannel,
					FEM_ObjectBroker &theBroker){return -1;}


  protected:
  double val_hist, val_trial;
  double residual, sFactor;
    
  private:
};


#endif

