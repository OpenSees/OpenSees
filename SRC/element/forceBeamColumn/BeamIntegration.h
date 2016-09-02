/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */

// $Revision: 1.7 $
// $Date: 2006-09-05 22:57:11 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/forceBeamColumn/BeamIntegration.h,v $

#ifndef BeamIntegration_h
#define BeamIntegration_h

#include <OPS_Globals.h>
#include <MovableObject.h>
#include <TaggedObject.h>
#include <ID.h>

class Matrix;
class ElementalLoad;
class Information;

class BeamIntegration : public MovableObject
{
 public:
  BeamIntegration(int classTag);
  virtual ~BeamIntegration();

  virtual void getSectionLocations(int nIP, double L, double *xi) = 0;
  virtual void getSectionWeights(int nIP, double L, double *wt) = 0;

  virtual void addElasticDeformations(ElementalLoad *theLoad,
				      double loadFactor,
				      double L, double *v0) {return;}
  // Return 0 if there is no elastic interior, -1 otherwise
  virtual int addElasticFlexibility(double L, Matrix &fe) {return 0;}

  virtual double getTangentDriftI(double L, double LI, double q2,
				  double q3, bool yAxis = false) {return 0.0;}
  virtual double getTangentDriftJ(double L, double LI, double q2,
				  double q3, bool yAxis = false) {return 0.0;}

  virtual BeamIntegration *getCopy(void) = 0;

  virtual void getLocationsDeriv(int nIP, double L, double dLdh,
				 double *dptsdh);
  virtual void getWeightsDeriv(int nIP, double L, double dLdh,
			       double *dwtsdh);
  // Return 0 if there is no elastic interior, -1 otherwise
  virtual int addElasticFlexDeriv(double L, Matrix &dfedh,
				  double dLdh = 0.0) {return 0;}

  virtual void Print(OPS_Stream &s, int flag = 0) = 0;
};

// a BeamIntegrationRule store BeamIntegration and section tags
class BeamIntegrationRule : public TaggedObject
{
public:
    BeamIntegrationRule(int tag, BeamIntegration* bi, const ID& stags)
	:TaggedObject(tag),theInt(bi),secTags(stags){}
    ~BeamIntegrationRule(){
	if (theInt != 0) delete theInt;
    }

    BeamIntegration* getBeamIntegration(){return theInt;}
    const ID& getSectionTags() const {return secTags;}

    void Print(OPS_Stream &s, int flag){return;}
private:
    BeamIntegration* theInt;
    ID secTags;
};

bool       OPS_addBeamIntegrationRule(BeamIntegrationRule *newComponent);
BeamIntegrationRule *OPS_getBeamIntegrationRule(int tag);
void       OPS_clearAllBeamIntegrationRule(void);

#endif
