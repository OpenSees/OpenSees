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

#ifndef InitialStateAnalysisWrapper_h
#define InitialStateAnalysisWrapper_h

// Written: Chris McGann
//          February 2011

// Description: This file contains the class definition for InitialStateAnalysisWrapper.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <NDMaterial.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>

class InitialStateAnalysisWrapper : public NDMaterial
{
	public:

		// full constructor
		InitialStateAnalysisWrapper(int tag, NDMaterial &mainMat, int ndim);

		// null constructor
		InitialStateAnalysisWrapper();

		// destructor
		~InitialStateAnalysisWrapper();

		int commitState(void);
		int revertToLastCommit(void);
		int revertToStart(void);
		
		// set the strain to be sent to the main material
		int setTrialStrain(const Vector &strain_from_element);

		// get mass density from main material
		double getRho(void);

		// send back strain
		const Vector& getStrain();
		// send back stress
		const Vector& getStress();
		// send back the tangent
		const Matrix& getTangent();
		const Matrix& getInitialTangent();

		NDMaterial *getCopy(const char *type);
		NDMaterial *getCopy(void);
		const char *getType(void) const;
		int getOrder(void) const;

		Response *setResponse (const char **argv, int argc, OPS_Stream &output);
    	int getResponse (int responseID, Information &matInformation);

    	int sendSelf(int commitTag, Channel &theChannel);  
    	int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker); 

    	void Print(OPS_Stream &s, int flag =0);

	int setParameter(const char **argv, int argc, Parameter &param);
	int updateParameter(int responseID, Information &eleInformation);
	
	friend class PyLiq1;
	friend class TzLiq1;
	friend class QzLiq1; // Sumeet

	int getMainClassTag();           // sends class tag of main material object
	
 protected:
	
	NDMaterial *theMainMaterial;     // pointer to main material object

		//int getMainClassTag();           // sends class tag of main material object

		// input variables
		int mDIM;                        // number of dimensions in problem
		
		// member variables
		Vector mEpsilon_o;               // initial strain stored here
		Vector mStrain;                  // strain sent to the main material
};
#endif
