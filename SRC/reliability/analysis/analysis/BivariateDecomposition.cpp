/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 2001, The Regents of the University of California    **
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
** Reliability module developed by:                                   **
**   Terje Haukaas (haukaas@ce.berkeley.edu)                          **
**   Armen Der Kiureghian (adk@ce.berkeley.edu)                       **
**                                                                    **
**   Quan Gu (qgu@ucsd.edu)                                           **
**   Joel P. Conte (jpconte@ucsd.edu)                                 **
** ****************************************************************** */
                                                                        
 
//
// Written by  Quan Gu UCSD
//
 
#include <math.h>
#include "BivariateDecomposition.h"
//#include <fstream.h>
 
BivariateDecomposition::BivariateDecomposition( int pNumAxis,/* Vector * pGradient,*/ PrincipalAxis ** pPrincipalAxes, GridPlane ** thePlane, bool pIsTimeVariant)  
//  numAxes is the number of principalAxes+1
{
	
	numAxes = pNumAxis;
	isTimeVariant  = pIsTimeVariant;
	
//	if (pGradient ==0) gradientDP =0;
//	else gradientDP = new Vector(*pGradient);

    thePrincipalAxes = pPrincipalAxes;
	theGridPlanes = thePlane;

	strcpy(type,"BivariateDecomposition");

	int numberOfPlanes = numAxes*(numAxes-1)/2;   

	coefficients = new Matrix * [numberOfPlanes]; 
	int i;
	for(i=0; i<numberOfPlanes; i++)
		coefficients[i]=0;

	axisCoeff = new Vector * [numAxes];
	for(  i=0; i<numAxes; i++)
		axisCoeff[i]=0;



	if ( isTimeVariant){
		coefficients_2 = new Matrix * [numberOfPlanes]; 
		int i;
		for(i=0; i<numberOfPlanes; i++)
			coefficients_2[i]=0;

		axisCoeff_2 = new Vector * [numAxes];
		for( i=0; i<numAxes; i++)
			axisCoeff_2[i]=0;

	}


	tmp1=0;
	linearCorrection =0.0;

	

}

BivariateDecomposition::~BivariateDecomposition()
{
//   if (gradientDP !=0) delete gradientDP;
   for(int i=0; i<numAxes; i++){
		if (coefficients[i] !=0) delete coefficients[i];
   }
   delete [] coefficients;

   if (isTimeVariant){
   
	   for(int i=0; i<numAxes; i++){
			if (coefficients_2[i] !=0) delete coefficients_2[i];
	   }
	   delete [] coefficients_2;
   }

   if (tmp1 !=0 ) delete tmp1;
}

/*

  coefficients pointer matrix: ii


                  pointer:     0    1    2   .....  numberOfGridPlanes-1
                                                         |
                             ____________________________|__
       coefficients:         |    |    |    |    |    |  | |  
                    |--------|- 0 | 1  |  2 | .. |    |    |  
                    |        |____|____|____|____|____|____|  
                    |
                   \|/  corresponding to gridPlane[0]
			 __________________________________________________________________ ---------------> x2: n points
			|            |            |            |            |             |
			|            |            |            |            |             |
			|   c(0,0)   |  c(0,1)    |  c(0,2)    |    ...     |  c(0,n-1)   |
	0       |  x1^(m-1)* |  x1^(m-1)* |  x1^(m-1)* |            |  x1^(m-1)*  |
			|  x2^(n-1)  |  x2^(n-2)  |  x2^(n-3)  |            |  x2^0       |
            |____________|____________|____________|____________|_____________|
			|            |            |            |            |             |
			|            |            |            |            |             |
	1		|   c(1,0)   |  c(1,1)    |  c(1,2)    |    ...     |  c(1,n-1)   |
			|   x1^(m-2)*|  x1^(m-2)* |  x1^(m-2)* |            |  x1^(m-2)*  |
			|   x2^(n-1) |  x2^(n-2)  |  x2^(n-3)  |            |  x2^0       |
 			|____________|____________|____________|____________|_____________|
			|            |            |            |            |             |
			|            |            |            |            |             |
	2		|   c(2,0)   |  c(2,1)    |  c(2,2)    |    ...     |  c(2,n-1)   |
			|   x1^(m-3)*|  x1^(m-3)* |  x1^(m-3)* |            |  x1^(m-3)*  |
			|   x2^(n-1) |  x2^(n-2)  |  x2^(n-3)  |            |  x2^0       |
 			|____________|____________|____________|____________|_____________|
			|            |            |            |            |             |
			|            |            |            |            |             |
			|      ..    |            |            |            |  ..         |
  		    |            |            |            |            |             |
			|            |            |            |            |             |
			|____________|____________|____________|____________|_____________|
			|            |            |            |            |             |
			|            |            |            |            |             |
	m-1		|   c(m-1,0) |  c(m-1,1)  |  c(m-1,2)  |    ...     |  c(m-1,n-1) |
			|   x1^0*    |  x1^0    * |  x1^0*     |            |  x1^0*      |
			|   x2^(n-1) |  x2^(n-2)  |  x2^(n-3)  |            |  x2^0       |
			|____________|____________|____________|____________|_____________|

		 |
         |        0             1           2                          n-1
         |
        \|/ 
        x1: m points



          m points in x2 direction, n points in x1 direction
	
                   m-1  n-1
      gFunValue =  SUM  SUM [ c(i,j) * x1^(m-1-i) * x2^(n-1-j)]
                   i=0  j=0        



*/

int BivariateDecomposition::fitCurve() {  // Lagrange interpolation . refer "Decomposition methods for structural reliability analysis" H. Xu, etc. 2005

  int ii;
	for(ii=0; ii<numAxes; ii++){ // for a special axis
		thePrincipalAxes[ii]->computeShapeFuncCoeff();
	}

//* --- deal  with axes --------
	
	ofstream tt("point_shape_coeff.out");


	

	for(   ii =0; ii<numAxes; ii++){
//		tt << "ii="<<ii<<endln;	
		if (axisCoeff[ii] !=0) delete axisCoeff[ii];

		int m = thePrincipalAxes[ii]->getExperimentalPointRule()->getNumberOfPoints();
		axisCoeff[ii] = new Vector(m);
		axisCoeff[ii]->Zero();


		if (isTimeVariant){ //---
			
			if (axisCoeff_2[ii] !=0) delete axisCoeff_2[ii];
			axisCoeff_2[ii] = new Vector(m);
			axisCoeff_2[ii]->Zero();
		
		
		} //---
		
		for ( int point=0; point < m; point++){
			Vector * theCoeff = thePrincipalAxes[ii]->getShapeFuncCoeff(point);
			axisCoeff[ii]->addVector(1.0, *theCoeff, thePrincipalAxes[ii]->getValueOnAxis(point));

			if (isTimeVariant){
				axisCoeff_2[ii]->addVector(1.0, *theCoeff, thePrincipalAxes[ii]->getValueG2OnAxis(point));			
			}
			tt<<"point num="<<point<<", shape_coeff:"<<endln;
			for(int i =0; i<theCoeff->Size(); i++)
				tt<<(*theCoeff)(i)<<endln;

		} // for each point
	
	} // for ii
	tt.close();  
	
	//* -- debug 
	ofstream file_coeff("Bivariate_coefficient_recorder.out");
	for( ii=0; ii<numAxes; ii++){
		opserr<<"\n\axisCoeff["<<ii<<"]"<<endln;
		file_coeff<<"axisCoeff["<<ii<<"]"<<endln;

		for (int m=0; m<(*axisCoeff[ii]).Size(); m++){
			opserr<<(*axisCoeff[ii])(m)<<endln;
			file_coeff<<(*axisCoeff[ii])(m)<<endln;
		}
	}
//*/
// --- deal with gridPlanes ---	


	int numOfGridPlanes = numAxes*(numAxes-1)/2;

	int m,n;
	for (ii=0; ii<numOfGridPlanes; ii++) {
		PrincipalAxis * axis_1 = theGridPlanes[ii]->getAxisPtr(1); // first axis
		PrincipalAxis * axis_2 = theGridPlanes[ii]->getAxisPtr(2);  // second axis

		 m = axis_1->getExperimentalPointRule()->getNumberOfPoints();
		 n = axis_2->getExperimentalPointRule()->getNumberOfPoints();

		if (coefficients[ii] !=0) delete coefficients[ii];

		coefficients[ii] = new Matrix(m,n);
		coefficients[ii]->Zero();

		if (isTimeVariant) {   // add time variant----
			if (coefficients_2[ii] !=0) delete coefficients_2[ii];
			coefficients_2[ii] = new Matrix(m,n);
			coefficients_2[ii]->Zero();	
		} //-------


		Matrix * valueOnGrid = theGridPlanes[ii]->getGridValuesPtr();
		Matrix * valueG2OnGrid=0;
		
		if(isTimeVariant){
			valueG2OnGrid = theGridPlanes[ii]->getGridValuesG2Ptr();
		}

        for( int pointX = 0; pointX < m; pointX++){

			Vector * shapeFunCoeff_1 = axis_1->getShapeFuncCoeff(pointX);
			for( int pointY = 0; pointY < n; pointY++){
			    
 
				Vector * shapeFunCoeff_2 = axis_2->getShapeFuncCoeff(pointY); 

				for(int i=0; i<m; i++) { //x
					for(int j=0; j<n; j++){ //y

						(*coefficients[ii])(i,j) += (*shapeFunCoeff_1)(i)*(*shapeFunCoeff_2)(j)*(*valueOnGrid)(pointX,pointY); 
						if (isTimeVariant) {   // add time variant----
							(*coefficients_2[ii])(i,j) += (*shapeFunCoeff_1)(i)*(*shapeFunCoeff_2)(j)*(*valueG2OnGrid)(pointX,pointY); 
						}
					
					} //j
				} //i
			} //pointY
		}//pointX
	} //ii =1:numOfGridPlanes
		
  //* --- debug	
	for( ii=0; ii<numOfGridPlanes; ii++){

		Matrix * valueOnGrid = theGridPlanes[ii]->getGridValuesPtr();
		Matrix * valueG2OnGrid=0;
		
		if(isTimeVariant){
			valueG2OnGrid = theGridPlanes[ii]->getGridValuesG2Ptr();
		}


		opserr<<"\n\n\n\ncoefficients["<<ii<<"]"<<endln;
		file_coeff<<"plane +++++++++++"<<ii<<"+++++++++++++++"<<endln;

		file_coeff<<"values on grid ==============================="<<endln;

		int i;
		for (i=0; i<m; i++){
			for( int j=0; j<n; j++){
				opserr<<"values"<<i<<","<<j<<":"<<(*valueOnGrid)(i,j)<<endln;
				file_coeff<<"values"<<i<<","<<j<<":"<<(*valueOnGrid)(i,j)<<endln;
			} //j
		}//i

		file_coeff<<"coeff on grid ==============================="<<endln;



		for (i=0; i<m; i++){
			for( int j=0; j<n; j++){
				opserr<<"coefficients"<<i<<","<<j<<":"<<(*coefficients[ii])(i,j)<<endln;
				file_coeff<<"coefficients"<<i<<","<<j<<":"<<(*coefficients[ii])(i,j)<<endln;
			} //j
		}//i
	}

	file_coeff.close();
//*/


	return 0;

};

double BivariateDecomposition::getFunctionValue(Vector * point){

	int numOfGridPlanes = numAxes*(numAxes-1)/2;
	double gFunValue=0;
	double constant =0;
	

	int m,n;
	int ii;
	for (ii =0; ii< numOfGridPlanes; ii++){

		PrincipalAxis * axis_1 = theGridPlanes[ii]->getAxisPtr(1); // first axis
		PrincipalAxis * axis_2 = theGridPlanes[ii]->getAxisPtr(2);  // second axis


		Vector * vector1 = axis_1->getAxisDirection();
		Vector * vector2 = axis_2->getAxisDirection();

		double x =  (*point)^(*vector1)/(vector1->Norm());
		double y =  (*point)^(*vector2)/(vector2->Norm());

		m = axis_1->getExperimentalPointRule()->getNumberOfPoints();
		n = axis_2->getExperimentalPointRule()->getNumberOfPoints();



		for (int i =0; i<m; i++){
			for (int j =0; j<n; j++){
				gFunValue += (*coefficients[ii])(i,j) * pow(x,m-1-i)*pow(y,n-1-j);	
						
			}
		}

	}// ii =0:numOfGridPlanes
	constant = (*coefficients[numOfGridPlanes-1])(m-1,n-1);
//	opserr<<"constant = "<<constant<<endln;


	double valueOnAxis =0.0;
    for (ii =0; ii< numAxes; ii++){

		Vector * direction = thePrincipalAxes[ii]->getAxisDirection();
		double x = (*point)^(*direction)/direction->Norm();
		int n = (*axisCoeff[ii]).Size();
		for (int j=0; j<n; j++){  // each power
			valueOnAxis += (*axisCoeff[ii])(j) * pow(x,n-j-1);
		
		} //j

	} //for
	gFunValue -= (numAxes-2)*valueOnAxis;


	gFunValue += (numAxes-1)*(numAxes-2)/2.0*constant;  // the constant is zero.


   return gFunValue;
};


double BivariateDecomposition::getFunctionValue2(Vector * point, Vector * dp2prime, Vector * gradG2){

	int numOfGridPlanes = numAxes*(numAxes-1)/2;
	double gFunValue2=0;
	double constant =0;
	
// --- grid planes ---
	int m,n;
	int ii;
	for (ii =0; ii< numOfGridPlanes; ii++){

		PrincipalAxis * axis_1 = theGridPlanes[ii]->getAxisPtr(1); // first axis
		PrincipalAxis * axis_2 = theGridPlanes[ii]->getAxisPtr(2);  // second axis


		Vector * vector1 = axis_1->getAxisDirection();
		Vector * vector2 = axis_2->getAxisDirection();

		double x =  (*point)^(*vector1)/(vector1->Norm());
		double y =  (*point)^(*vector2)/(vector2->Norm());

		m = axis_1->getExperimentalPointRule()->getNumberOfPoints();
		n = axis_2->getExperimentalPointRule()->getNumberOfPoints();



		for (int i =0; i<m; i++){
			for (int j =0; j<n; j++){
				gFunValue2 += (*coefficients_2[ii])(i,j) * pow(x,m-1-i)*pow(y,n-1-j);	
						
			}
		}

	}// ii =0:numOfGridPlanes
	constant = (*coefficients_2[numOfGridPlanes-1])(m-1,n-1);

// --- grid axis ----

	    double valueOnAxis =0.0;
    for (ii =0; ii< numAxes; ii++){

		Vector * direction = thePrincipalAxes[ii]->getAxisDirection();
		double x = (*point)^(*direction)/direction->Norm();
		int n = (*axisCoeff[ii]).Size();
		for (int j=0; j<n; j++){  // each power
			valueOnAxis += (*axisCoeff_2[ii])(j) * pow(x,n-j-1);
		
		} //j

	}
	gFunValue2 -= (numAxes-2)*valueOnAxis;

// --- origin ---
	gFunValue2 += (numAxes-1)*(numAxes-2)/2.0*constant;  // the constant is zero.

// ============ linear contribution ========

	// -- remove the first numAxes th contribution ---

	Vector reducedPt(*point);

	reducedPt.addVector(1.0, *dp2prime, -1.0);

	for ( int i=0; i<numAxes; i++){ // each axis
		Vector theAxis( *(thePrincipalAxes[i]->getAxisDirection()));

		double norm = theAxis.Norm();
		theAxis /=norm;

		double x = reducedPt^theAxis;
			   
		reducedPt.addVector(1.0,theAxis, -x); 
	
	}


	double tmp = reducedPt^(*gradG2);

    gFunValue2 += tmp; 

	if (linearCorrection ==0) {   // linear Correction = deltG * (u2*-u1*) but without nonlinear-ndir contributions
		Vector reducedDP2Prime(*dp2prime);
		for (int ii =0; ii< numAxes; ii++){
				Vector * direction = thePrincipalAxes[ii]->getAxisDirection();
				double x = (*dp2prime)^(*direction)/direction->Norm();
				reducedDP2Prime.addVector(1.0, *direction, -x/direction->Norm());
		}
		
		linearCorrection = reducedDP2Prime^(*gradG2);


	
	} //if (linearCorrection ==0)

	
 /*  if (linearCorrection ==0) { 

		linearCorrection = (*dp2prime)^(*gradG2);

   }*/
	gFunValue2 += linearCorrection; 
	
	
	
	


    return gFunValue2; 
 

};

char * BivariateDecomposition::getType(){
	return type;
};

Vector * BivariateDecomposition::Poly(Vector *x)
{


/*	% Expand recursion formula
for j=1:n
    
    for k=2:j+1
        c(k) = c1(k) - e(j)*c1(k-1);
    end
    
     c1 = c;
end

*/

	int n = x->Size();

	if (tmp1 ==0) {
		tmp1 = new Vector(n+1);
	}
	else if (tmp1->Size() != n+1){
		delete tmp1;
		tmp1 = new Vector(n+1);
	}

		
	tmp1->Zero();
	(*tmp1)(0)=1.0;

	Vector  tmp_old(*tmp1);

	for (int j=0; j<=n-1; j++){
		for (int k =1; k<=j+1; k++){  
			(*tmp1)(k) = tmp_old(k) - (*x)(j)*tmp_old(k-1);
		}
		tmp_old.addVector(0.0, *tmp1,1.0);
	}

	
	return tmp1;
}

/*
void BivariateDecomposition::setGradient(Vector *gradient)
{
	if (gradientDP != 0) delete gradientDP;
	gradientDP = new Vector(*gradient); 
}
*/
void BivariateDecomposition::setPincipalAxesPtr(PrincipalAxis ** pAxis)
{
	thePrincipalAxes = pAxis;
}

void BivariateDecomposition::setGridPlanesPtr(GridPlane **theGPlane)
{
	theGridPlanes = theGPlane;
}




double BivariateDecomposition::debug(Vector * point, Vector * dp2prime, Vector * gradG2){


/////////////////////////////// ------------------- for debug only by Quan ----

	static int iii =1; 
	if (iii==1){
		iii++;
		ofstream file_coeff_2("Bivariate_surface_fitting_debug.out");

     	int numOfGridPlanes = numAxes*(numAxes-1)/2;

		for (int ii=0; ii<numOfGridPlanes; ii++) {

			file_coeff_2 <<"\n the "<<ii<<"  th plane: "<<endln;
			PrincipalAxis * axis_1 = theGridPlanes[ii]->getAxisPtr(1); // first axis
			PrincipalAxis * axis_2 = theGridPlanes[ii]->getAxisPtr(2);  // second axis

			 int m = axis_1->getExperimentalPointRule()->getNumberOfPoints();
			 int n = axis_2->getExperimentalPointRule()->getNumberOfPoints();




			Matrix * valueOnGrid = theGridPlanes[ii]->getGridValuesPtr();
			Matrix * valueG2OnGrid=0;
			
			if(isTimeVariant){
				valueG2OnGrid = theGridPlanes[ii]->getGridValuesG2Ptr();
			}

			double x_begin = axis_1->getExperimentalPointRule()->getBeginOfGrid();
			double x_end =   axis_1->getExperimentalPointRule()->getEndOfGrid();
			double x_interval = (x_end-x_begin)/(m-1);

			double y_begin=  axis_2->getExperimentalPointRule()->getBeginOfGrid();
			double y_end =   axis_2->getExperimentalPointRule()->getEndOfGrid();
			double y_interval = (y_end-y_begin)/(n-1);

			Vector * vector1 = axis_1->getAxisDirection();
			Vector * vector2 = axis_2->getAxisDirection();

			Vector tmpp(50); tmpp.Zero();
			file_coeff_2.precision(16);
			for( int pointX = 0; pointX < m; pointX++){
   

				for( int pointY = 0; pointY < n; pointY++){
					
					double x = x_begin + pointX*x_interval;
					double y = y_begin + pointY*y_interval;
					

					tmpp = x*(*vector1)/vector1->Norm()+y*(*vector2)/vector2->Norm();



					file_coeff_2<<"x:"<<x<<"y:"<<y<<endln;
					file_coeff_2<<"grid value G saved in the valueOnGrid is: "<<(*valueOnGrid)(pointX,pointY)<<endln;
					double gValue = this->getFunctionValue(&tmpp);
					file_coeff_2<<"grid value G computed by the response surface is: "<<gValue<<endln;

					file_coeff_2<<"grid value G2 saved in the valueOnGrid is: "<<(*valueG2OnGrid)(pointX,pointY)<<endln;
					double g2Value = this->getFunctionValue2(&tmpp,dp2prime,gradG2);
					file_coeff_2<<"grid value G2 computed by the response surface is: "<<g2Value<<endln;

	  			 

				} //pointY
			}//pointX
		} //ii =1:numOfGridPlanes
			
 
	//*/
		file_coeff_2.close();

	}// if iii==1     
////////////////////------ end debug -----------------
    return 0; 
 

};
