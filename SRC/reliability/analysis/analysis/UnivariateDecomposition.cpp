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
 

#include "UnivariateDecomposition.h"
//#include <fstream.h>
#include <math.h>

UnivariateDecomposition::UnivariateDecomposition( int pNumAxis, /*Vector * pGradient,*/ PrincipalAxis ** pPrincipalAxes, bool pIsTimeVariant)  //  numAxes is the number of principalAxes+1
{
	
	numAxes = pNumAxis;
	isTimeVariant  = pIsTimeVariant;
	
//	if (pGradient ==0) gradientDP =0;
//	else gradientDP = new Vector(*pGradient);

    thePrincipalAxes = pPrincipalAxes;

	strcpy(type,"UnivariateDecomposition");

	coefficients = new Vector * [numAxes]; 
	
	for(int i=0; i<numAxes; i++)
		coefficients[i]=0;

	if ( isTimeVariant){
		coefficients_2 = new Vector * [numAxes]; 
		for(int i=0; i<numAxes; i++)
			coefficients_2[i]=0;

	}


	tmp=0;
    linearCorrection =0.0; 
	

}

UnivariateDecomposition::~UnivariateDecomposition()
{
   //if (gradientDP !=0) delete gradientDP;
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

   if (tmp !=0 ) delete tmp;
}

/*

  coefficients pointer array: ii
  pointer ii: 0    1    2    3    4    5    6    7    8    9  size = numOfPrincipalAxis+1 = size of thePrincipalAxes          
  
			 _________________________________________________ 
			|    |    |    |    |    |    |    |    |    |    |
			| x  | y  | z  |... |    |    |    |    |    |    |
			|____|____|____|____|____|____|____|____|____|____|
              |    |    |
              |    |    |
              |    |    |
              |    |    ----> Vector(n) save coefficients of this dim.  n = Num of points in this dim.   
              |    |    
              |    |           0    1 ...               n-1
			  |   \|/        _______________________________
              |  1st_eigen   |    |    |    |    |    |    |  
              |              | c0 |c1  | c2 | .. |    |    |  
              |              |____|____|____|____|____|____|  
              |            
             \|/           polynormal = c0*z^(n-1) + c1*z^(n-2) + ....+c(n-1) 
             DP direction

*/

int UnivariateDecomposition::fitCurve() {  // Lagrange interpolation . refer "Decomposition methods for structural reliability analysis" H. Xu, etc. 2005


	

	ofstream tt("point_shape_coeff_uni.out");


	
	int ii;
	for(ii =0; ii<numAxes; ii++){
		tt << "ii="<<ii<<endln;	
		if (coefficients[ii] !=0) delete coefficients[ii];

		int m = thePrincipalAxes[ii]->getExperimentalPointRule()->getNumberOfPoints();
		coefficients[ii] = new Vector(m);
		coefficients[ii]->Zero();


		if (isTimeVariant){ //---
			
			if (coefficients_2[ii] !=0) delete coefficients_2[ii];
			coefficients_2[ii] = new Vector(m);
			coefficients_2[ii]->Zero();
		
		
		} //---
		
		for ( int point=0; point < m; point++){
			Vector * theCoeff = thePrincipalAxes[ii]->getShapeFuncCoeff(point);
			coefficients[ii]->addVector(1.0, *theCoeff, thePrincipalAxes[ii]->getValueOnAxis(point));

			if (isTimeVariant){
				coefficients_2[ii]->addVector(1.0, *theCoeff, thePrincipalAxes[ii]->getValueG2OnAxis(point));			
			}
			tt<<"point num="<<point<<", shape_coeff:"<<endln;
			for(int i =0; i<theCoeff->Size(); i++)
				tt<<(*theCoeff)(i)<<endln;

		} // for each point
	
	} // for ii
	tt.close();
	
 


//* -- debug 
	ofstream file_coeff("univariate_coefficient_recorder.out");
	for( ii=0; ii<numAxes; ii++){
		opserr<<"\n\ncoefficients["<<ii<<"]"<<endln;
		file_coeff<<"coefficients["<<ii<<"]"<<endln;

		for (int m=0; m<(*coefficients[ii]).Size(); m++){
			opserr<<(*coefficients[ii])(m)<<endln;
			file_coeff<<(*coefficients[ii])(m)<<endln;
		}
	}
//

	if (isTimeVariant){  // add time variant
	
		for( ii=0; ii<numAxes; ii++){
		opserr<<"\n\ncoefficients_2["<<ii<<"]"<<endln;
		file_coeff<<"coefficients_2["<<ii<<"]"<<endln;

		for (int m=0; m<(*coefficients[ii]).Size(); m++){
			opserr<<(*coefficients_2[ii])(m)<<endln;
			file_coeff<<(*coefficients_2[ii])(m)<<endln;
		}
		} // time variant
	
	
	}

	file_coeff.close();



	return 0;

};

double UnivariateDecomposition::getFunctionValue(Vector * point){

	int size = point->Size();
	
	double gFunValue=0;
	for (int i=0; i<numAxes; i++){ // each axis

		int n = coefficients[i]->Size();
		double norm = thePrincipalAxes[i]->getAxisDirection()->Norm();	
		double x = (*point)^(*thePrincipalAxes[i]->getAxisDirection());
		       x /= norm;

		for (int j=0; j<n; j++){  // each power
			gFunValue += (*coefficients[i])(j) * pow(x,n-j-1);
		
		}
	
	}

    return gFunValue;

};


double UnivariateDecomposition::getFunctionValue2(Vector * point, Vector * dp2prime, Vector * gradG2){

	int size = point->Size();
	
	double gFunValue2=0;
	double constant =0.0;
	int n =0;
	int i;
	for (i=0; i<numAxes; i++){ // each axis

		n = coefficients_2[i]->Size();
		double norm = thePrincipalAxes[i]->getAxisDirection()->Norm();	
		double x = (*point)^(*thePrincipalAxes[i]->getAxisDirection());
		       x /= norm;

		for (int j=0; j<n; j++){  // each power
			gFunValue2 += (*coefficients_2[i])(j) * pow(x,n-j-1);
		
		}
	
	}

	constant = (*coefficients_2[numAxes-1])(n-1);
	
// --------- linear contribution ----

	// -- remove the first numAxes th contribution ---

	Vector reducedPt(*point);

	reducedPt.addVector(1.0, *dp2prime, -1.0);

	for (  i=0; i<numAxes; i++){ // each axis
		Vector theAxis( *(thePrincipalAxes[i]->getAxisDirection()));

		double norm = theAxis.Norm();
		theAxis /=norm;

//		double x = (*point)^theAxis;  // ---------old ***********************
		
		double x = reducedPt^theAxis;
		
		reducedPt.addVector(1.0,theAxis, -x); 
	
	}


	double tmp = reducedPt^(*gradG2);

    gFunValue2 += tmp; 

	gFunValue2 -= (numAxes-1)*constant;

	//*********	if (debug == true ) {   // linear Correction = deltG * (u2*-u1*) but without nonlinear-ndir contributions 
			
/*		//	=== debug ===
		//Vector reducedDP2Prime = Vector(*dp2prime);
		ofstream file_coeff_2("Univariate_surface_fitting_linear_debug.out");
		file_coeff_2.precision(16);
       file_coeff_2<<"================\n==================\n================"<<endln;
       opserr<<"================\n==================\n================"<<endln;

       int len=(*point).Size();

	   file_coeff_2<<"point=[";
	   for (  i=0; i<len; i++) 
	      file_coeff_2<<(*point)(i)<<", "; 
	   file_coeff_2<<"];"<<endln;

	   
	   file_coeff_2<<"dp2prime=[";
	   for ( i=0; i<len; i++) 
	      file_coeff_2<<(*dp2prime)(i)<<", "; 
	   file_coeff_2<<"];"<<endln;

	   file_coeff_2<<"Axis_1=[";
	   for ( i=0; i<len; i++) 
	      file_coeff_2<<(*(thePrincipalAxes[0]->getAxisDirection()))(i)<<", "; 
	   file_coeff_2<<"];"<<endln;

	   file_coeff_2<<"Axis_2=[";
	   for ( i=0; i<len; i++) 
	      file_coeff_2<<(*(thePrincipalAxes[1]->getAxisDirection()))(i)<<", "; 
	   file_coeff_2<<"];"<<endln;

	   file_coeff_2<<"Axis_3=[";
	   for ( i=0; i<len; i++) 
	      file_coeff_2<<(*(thePrincipalAxes[2]->getAxisDirection()))(i)<<", "; 
	   file_coeff_2<<"];"<<endln;

	   file_coeff_2<<"gradG2=[";
	   for ( i=0; i<len; i++) 
	      file_coeff_2<<(*gradG2)(i)<<", "; 
	   file_coeff_2<<"];"<<endln;

	   file_coeff_2<<"reducedPt=[";
	   for ( i=0; i<len; i++) 
	      file_coeff_2<<reducedPt(i)<<", "; 
	   file_coeff_2<<"];"<<endln;

	   file_coeff_2<<"tmp = "<<tmp<<";"<<endln;




	   opserr<<"*dp2prime=["<<*dp2prime<<endln;
	   opserr<<"Axis_1=["<<*(thePrincipalAxes[0]->getAxisDirection())<<"];"<<endln;
	   opserr<<"Axis_2=["<<*(thePrincipalAxes[1]->getAxisDirection())<<"];"<<endln;
	   opserr<<"Axis_3=["<<*(thePrincipalAxes[2]->getAxisDirection())<<"];"<<endln;
	   opserr<<"reducedPt=["<<reducedPt<<"];"endln;
	   opserr<<"*gradG2=["<<*gradG2<<"];"endln;

       //opserr<<"================\n==================\n================"<<endln;
	   //file_coeff_2<<"================\n==================\n================"<<endln;
	   

		Vector reducedDP2Prime(*dp2prime);
		for (int ii =0; ii< numAxes; ii++){
				Vector * direction = thePrincipalAxes[ii]->getAxisDirection();
				double x = (*dp2prime)^(*direction)/direction->Norm();
				reducedDP2Prime.addVector(1.0, *direction, -x/direction->Norm());
		}
		
		double linearCorrection_1 = reducedDP2Prime^(*gradG2);

	   file_coeff_2<<"linearCorrection_1 = "<<linearCorrection_1<<";"<<endln; //debug=====

	
//	} // ******************if (debug == true)

	*/
/*   if (linearCorrection ==0) {    // ---------old **********************

		linearCorrection = (*dp2prime)^(*gradG2);



   }  */


		
	if (linearCorrection ==0) { 
		Vector reducedDP2Prime(*dp2prime);
		for (int ii =0; ii< numAxes; ii++){
			Vector * direction = thePrincipalAxes[ii]->getAxisDirection();
			double x = (*dp2prime)^(*direction)/direction->Norm();
			reducedDP2Prime.addVector(1.0, *direction, -x/direction->Norm());
		}
			
		linearCorrection = reducedDP2Prime^(*gradG2);
	}
    //     file_coeff_2.close();  //debug ====
   
   gFunValue2 += linearCorrection; 

/////////////////////////////// ------------------- for debug only by Quan ----
/*
	
	static iii =1;
	if (iii==1){
		iii++;
		ofstream file_coeff_2("Univariate_surface_fitting_debug.out");
		file_coeff_2.precision(16);

		for (int ii=0; ii<numAxes; ii++) {

			file_coeff_2 <<"\n the "<<ii<<"  th Axis: "<<endln;


		int m = thePrincipalAxes[ii]->getExperimentalPointRule()->getNumberOfPoints();
			 

//			Vector * ValueOnAxis =  thePrincipalAxes[ii]->getValuesOnAxis();
//			Vector * ValueG2OnAxis =  0;
			
			

			double x_begin = thePrincipalAxes[ii]->getExperimentalPointRule()->getBeginOfGrid();
			double x_end =   thePrincipalAxes[ii]->getExperimentalPointRule()->getEndOfGrid();
			double x_interval = (x_end-x_begin)/(m-1);

 

		Vector * vector1 = thePrincipalAxes[ii]->getAxisDirection();
	 

			Vector tmpp(50); tmpp.Zero();
			
			for( int pointX = 0; pointX < m; pointX++){
   
  
					
					double x = x_begin + pointX*x_interval;
 

					tmpp = x*(*vector1)/vector1->Norm();



					file_coeff_2<<"x:"<<x<<endln;
					file_coeff_2<<"grid value G saved in the valueOnGrid is: "<<thePrincipalAxes[ii]->getValueOnAxis(pointX)<<endln;

					double gValue = this->getFunctionValue(&tmpp);
					file_coeff_2<<"grid value G computed by the response surface is: "<<gValue<<endln;



					file_coeff_2<<"grid value G2 saved in the valueOnGrid is: "<<thePrincipalAxes[ii]->getValueG2OnAxis(pointX)<<endln;
					double g2Value = this->getFunctionValue2(&tmpp,dp2prime,gradG2);
					file_coeff_2<<"grid value G2 computed by the response surface is: "<<g2Value<<endln;

	  			 

			 
			}//pointX
		} //ii =1:numAxes
			
 

		file_coeff_2.close();
		exit(-1);

	}// if iii==1     
	//*/	


    return gFunValue2;

};

char * UnivariateDecomposition::getType(){
	return type;
};


void UnivariateDecomposition::setPincipalAxesPtr(PrincipalAxis ** pAxis)
{
	thePrincipalAxes = pAxis;
}
