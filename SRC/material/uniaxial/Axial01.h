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
                                                                        
//
// $Date: 2016/05/29 12:19
// 
// Written: Erkan Bicici 
//
// Description: Axial deflection model with respect ot lateral deflection  
//
                                                                        
#ifndef Axial01_h
#define Axial01_h




#include <UniaxialMaterial.h>

class Element;
class Domain;
class DomainComponent;
class Node;

class Axial01 : public UniaxialMaterial
{
  public:
    Axial01(	int tag, 
			int nodeT, int nodeB, 
			double Kax, double Dsh, double Da, 
		Domain *theDomain, Node *theNodeT, Node *theNodeB);




    Axial01();    

    ~Axial01();

    int setTrialStrain(double strain, double strainRate = 0.0); 
    double getStrain(void);          
    double getStress(void);
    double getTangent(void);
	double getInitialTangent(void) { return 9.9e9; };
	double getLateralDrift(void);
	double getVerticalDrift(void);
	
	//void setDomain(Domain *theDomain);
    int commitState(void);
    int revertToLastCommit(void);    
    int revertToStart(void);    

    UniaxialMaterial *getCopy(void);
    
    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);    
    
    void Print(OPS_Stream &s, int flag =0);
    
  protected:
    
  private:
	  int dof;
	  int nodeTop;    // Top node of the column flexural section
	  int nodeBot;    // Bottom node of the column flexural section
	  double Kax;
	  double Dsh;
	  double Da;
	//  Domain *theDomain;
	  Domain *theDomain;
	  Node *theNodeT;
	  Node *theNodeB;

    double trialStrain;	// trial strain
    double trialStress;      // current trial stress
    double trialTangent;     // current trial tangent
    double commitStrain;     // last commited strain
    double commitStress;     // last commited stress
    double commitTangent;    // last committed  tangent
	double dx;
	double dy;
	double failure;
	double FailureStiff;

	double T_min_lt;
	double C_min_lt;
	double T_max_lt;
	double C_max_lt;

	double Clt;
	double Tlt;

	double TMinTangent;
	double CMinTangent;



	void MaxCheck(void);
};


#endif



