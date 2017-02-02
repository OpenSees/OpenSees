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
                                                                        
// Written: fmk 2017

#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 

#include <Vector.h>
#include <Matrix.h>
#include <ID.h>
#include <NDMaterial.h>

class PlaneStressLayeredMaterial : public NDMaterial {
  public : 
    PlaneStressLayeredMaterial();
    PlaneStressLayeredMaterial(int tag, 
			    int iLayers, 
			    double *thickness, 
			    NDMaterial **fibers);

    virtual ~PlaneStressLayeredMaterial();

    const char *getClassType(void) const {return "PlaneStressLayeredMaterial";};

    //mass per unit area
    double getRho() ;

    int setTrialStrain (const Vector &v);
    int setTrialStrain (const Vector &v, const Vector &r);
    int setTrialStrainIncr (const Vector &v);
    int setTrialStrainIncr (const Vector &v, const Vector &r);

    const Matrix &getTangent (void);
    const Matrix &getInitialTangent (void);

    const Vector &getStress (void);
    const Vector &getStrain (void);
        
    int commitState (void);
    int revertToLastCommit (void);
    int revertToStart (void);
    
    NDMaterial *getCopy (void);
    virtual NDMaterial *getCopy(const char *code);

    const char *getType (void) const {return "PlaneStressLayeredMaterial";};

    int getOrder (void) const;

    void Print( OPS_Stream &s, int flag ) ;
    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);    


    Response *setResponse (const char **argv, int argc, 
			   OPS_Stream &s);
    int getResponse (int responseID, Information &matInformation);

    
  private :
    double h ; // total thickness of section
    int nLayers;
    double *wg;

    NDMaterial **theFibers;  //pointers to the materials (fibers)

    Vector strain;
    static Vector stress;
    static Matrix tangent ;
    static ID array ;  

} ; //end of PlaneStressLayeredMaterial declarations





