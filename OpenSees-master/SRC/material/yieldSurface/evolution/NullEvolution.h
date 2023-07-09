/***************************************************************************
                          NullEvolution.h  -  description
                             -------------------
    begin                : Thu Aug 1 2002
    email                : rkaul@ce-blume215-pent-2.stanford.edu
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *                                                                         *
 ***************************************************************************/

#ifndef NULLEVOLUTION_H
#define NULLEVOLUTION_H

#include "YS_Evolution.h"

/**Useful for declaring inner-surfaces or pinching surfaces or just
plain elastic-perfectly plastic surfaces that do not evolve
  *@author rkaul
  */

class NullEvolution : public YS_Evolution  {
public: 

	NullEvolution(int tag, double isox);
	NullEvolution(int tag, double isox, double isoy);
	NullEvolution(int tag, double isox, double isoy, double isoz);
	~NullEvolution();

  int evolveSurface(YieldSurface_BC *ys, double magPlasticDefo,
                    Vector & G, Vector & F_Surface, int flag);

  const Vector & getEquiPlasticStiffness();

  YS_Evolution* getCopy();

  int getResponse(int responseID, Information & info);
  Response* setResponse(char **argv, int argc, OPS_Stream &output);
  
  int	 displaySelf(Renderer &theViewer, int displayMode, float fact) { return 0;}

  virtual int sendSelf(int commitTag, Channel &theChannel){return -1;}
  virtual int recvSelf(int commitTag, Channel &theChannel,
                         FEM_ObjectBroker &theBroker){return -1;}
  /** No descriptions */
  int revertToLastCommit();
  int commitState(int status);
  double getTrialPlasticStrains(int dof);
  double getCommitPlasticStrains(int dof);

private:
static Vector vec_dim_1;
static Vector vec_dim_2;  
static Vector vec_dim_3;
};

#endif
