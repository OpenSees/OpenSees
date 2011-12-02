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
                                                                        
// $Revision: 1.1.1.1 $
// $Date: 2000-09-15 08:23:21 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/ElasticSection3d.h,v $
                                                                        
                                                                        
///////////////////////////////////////////////////////
// File:  ~/Src/element/hinge/ElasticSection3d.h
//
// Written: MHS
// Date: May 2000
//
//
// Purpose:  This header file contains the prototype
// for the ElasticSection class.

#ifndef ElasticSection3d_h
#define ElasticSection3d_h

#include <SectionForceDeformation.h>
#include <Matrix.h>
#include <Vector.h>

class Channel;
class FEM_ObjectBroker;
class Information;

class ElasticSection3d: public SectionForceDeformation
{
  public:
	ElasticSection3d (int tag, double E, double A, double Iz, 
		double Iy, double G, double J);
    ElasticSection3d (int tag, double EA, double EIz, double EIy, double GJ);
	ElasticSection3d (void);
    ~ElasticSection3d (void);

    int commitState (void);
    int revertToLastCommit (void);
    int revertToStart (void);

    int setTrialSectionDeformation (const Vector&);
    const Vector &getSectionDeformation (void);

    const Vector &getStressResultant (void);
    const Matrix &getSectionTangent (void);
    const Matrix &getSectionFlexibility (void);
    
    SectionForceDeformation *getCopy (void);
    const ID &getType (void) const;
    int getOrder (void) const;
    
    int sendSelf (int commitTag, Channel &theChannel);
    int recvSelf (int commitTag, Channel &theChannel,
		  FEM_ObjectBroker &theBroker);
    
    void Print (ostream &s, int flag = 0);

  protected:

  private:
   
    double E, A, Iz, Iy, G, J;

    Matrix k;			// section stiffness matrix
    Matrix f;			// section flexibility matrix
    Vector e;			// section trial deformations
    Vector eCommit;

	static Vector s;			// section resisting forces, static for returns
		
    static ID code;
};

#endif
