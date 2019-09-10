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
 

#include <Hessian.h>
//#include <fstream>
//#include <iomanip>
#include <iostream>
using std::ifstream;
using std::ofstream;
using std::ios;
#include <math.h>
//using ios::append;
//using std::setprecision;
//using std::setiosflags;

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////


Hessian::Hessian(int pSize,ReliabilityDomain *passedReliabilityDomain, 
		 ProbabilityTransformation *passedProbabilityTransformation,
		 FunctionEvaluator *passedGFunEvaluator,
		 GradientEvaluator *passedGradGEvaluator, double tol)
{
  theProbabilityTransformation = passedProbabilityTransformation;
  theReliabilityDomain = passedReliabilityDomain;
  theGFunEvaluator = passedGFunEvaluator;
  theGradGEvaluator = passedGradGEvaluator;
  
  perturbTol=tol;
  standSens =0;
  designPointUSpace =0;
  
  if (pSize ==0) {
    opserr<<" Hessian: size reset to number of RVs."<<endln;
    sizeOfHessian = theReliabilityDomain->getNumberOfRandomVariables();
  }
  else sizeOfHessian = pSize;
  theHessian = new Matrix(sizeOfHessian,sizeOfHessian); // not good since hessian is symmetric. may change later..
  theHessian->Zero();
  theReducedHessian = new Matrix(sizeOfHessian-1,sizeOfHessian-1) ;   // in U space
  theReducedHessian->Zero();
  
  theHessianInPhysicalSpace = new Matrix(sizeOfHessian,sizeOfHessian); // not good since hessian is symmetric. may change later..
  theHessianInPhysicalSpace->Zero();
  normOfGradientInStandardNormalSpace =0.0;
  
  
  
}

Hessian::~Hessian()
{
  if (theHessian !=0) delete theHessian;
  if (theReducedHessian !=0) delete theReducedHessian;
  if (theHessianInPhysicalSpace !=0) delete theHessianInPhysicalSpace;
  if (designPointUSpace !=0) delete designPointUSpace;
}

const Matrix  &Hessian::getHessianApproximation(){
  return *theHessian; // U space
};


int 
Hessian::updateHessianApproximation(const Vector &u_old,
				    double g_old,
				    const Vector &gradG_old,
				    double stepSize,
				    const Vector &searchDirection,
				    double g_new,
				    const Vector &gradG_new)
{
  opserr << "Hessian::updateHessianApproximation() - QUAN - Not yet implemented\n";
  return -1;
}


int  Hessian::setHessianToIdentity(int size){

	if (size != sizeOfHessian) {
		opserr<<"Fatal: Hessian::setHessianToIdentity, size does not match!"<<endln;
		exit(-1);
	}

	theHessian->Zero();

	for(int i=0;i<size; i++)
		(*theHessian)(i,i)=1.0;
	theHessianInPhysicalSpace->addMatrix(0.0, *theHessian,1.0);

	return 0;
};

// Perturb in standard normal space by full perturbation..
int Hessian::formHessianByFDM(int numOfLimitStateFunction, Vector * theDesignPoint)
{
    Matrix jacobian_x_u(sizeOfHessian,sizeOfHessian);
	Matrix jacobian_u_x(sizeOfHessian,sizeOfHessian);
	//* --debug purpose --
			int ii=0;
			ifstream input_a_1( "hessian_restart_.out", ios::in); 
			input_a_1.precision(16);
			if (theHessian ==0) theHessian = new Matrix(sizeOfHessian, sizeOfHessian);
			//while (input_a_1.endof()){
			double tmp;
			int endoffile =-1;
			while((input_a_1.good())&&(ii<sizeOfHessian)&&(!input_a_1.eof())&&(endoffile ==-1)){ // ii is column
				
				int jj=0; 
				while ((jj<sizeOfHessian) && (endoffile ==-1)){
					if (!input_a_1.eof()){
						input_a_1 >> tmp;
						(*theHessian)(jj,ii)=tmp;
					}
					else{
						//opserr<<"hessian_restart_.out is wrong"<<endln;
						ii--; 
						endoffile =0;
					
					}
					jj++;
				}	
				ii++;	

			}// while 
							
			     
			

            input_a_1.close();
		//ii=100;

                 
 // -- debug --  */

   opserr<<"-- perturbation = "<<perturbTol<<endln;
//	ofstream debug("debiga.out"  );

	Vector x = (*theDesignPoint);

	/*
	// Transform starting point into standard normal space
	double result = theProbabilityTransformation->set_x(x);
	if (result < 0) {
		opserr << "SearchWithStepSizeAndStepDirection::doTheActualSearch() - " << endln
			<< " could not set x in the xu-transformation." << endln;
		return -1;
	}


	result = theProbabilityTransformation->transform_x_to_u();
	if (result < 0) {
		opserr << "SHessian::formHessianByFDM - " << endln
			<< " could not transform from x to u." << endln;
		return -1;

	}
	Vector u = theProbabilityTransformation->get_u();
	*/
	Vector u;
	int result = theProbabilityTransformation->transform_x_to_u(u);
	if (result < 0) {
		opserr << "SHessian::formHessianByFDM - " << endln
			<< " could not transform from x to u." << endln;
		return -1;
	}


	if (designPointUSpace ==0) 
		designPointUSpace = new Vector(u.Size());

	designPointUSpace->addVector(0.0, u,1.0);


	// -- compute and save Jacobi dx/du
	/*	result = theProbabilityTransformation->transform_u_to_x_andComputeJacobian();
	if (result < 0) {
		opserr << "Hessian::formHessianByFDM - " << endln
			<< " could not transform from u to x." << endln;
		return -1;
	}
	*/
	

	// --- compute standard sensitivity --
	//	jacobian_x_u.addMatrix(0.0,theProbabilityTransformation->getJacobian_x_u(),1.0);
	theProbabilityTransformation->getJacobian_x_to_u(jacobian_x_u);	

	
	//jacobian_u_x.addMatrix(0.0,theProbabilityTransformation->getJacobian_u_x(),1.0);
	theProbabilityTransformation->getJacobian_u_to_x(x, jacobian_u_x);	

	if (theGFunEvaluator->setVariables() < 0) {
	  opserr << "ERROR Hessian -- error setting variables in namespace" << endln;
	  return -1;
	}

	result = theGFunEvaluator->runAnalysis();
	if (result < 0) {
		opserr << "Hessian - " << endln
			<< " could not run analysis to evaluate limit-state function. " << endln;
		return -1;
	}
	double gFunctionValue = theGFunEvaluator->evaluateExpression();


	// Gradient in original space
	result = theGradGEvaluator->computeGradient(gFunctionValue);
	if (result < 0) {
		opserr << "Hessian - " << endln
			<< " could not compute gradients of the limit-state function. " << endln;
		return -1;
	}
	Vector gradientOfgFunction = theGradGEvaluator->getGradient();

	
/*	debug.precision(16);
	
	debug <<"dx_du=[\n";
	for(int i =0; i<x.Size(); i++){
		for(int j =0; j<x.Size(); j++){
			debug<<jacobian_x_u(i,j);
		}
		debug <<endln;
	}	
	debug <<"];"<<endln;

	
	debug <<"u_stand=[\n";
	for(i =0; i<x.Size(); i++){
		debug<<u(i)<<endln;
	}
	debug<<"]';"<<endln;

	debug <<"x_stand=[\n";
	for(i =0; i<x.Size(); i++){
		debug<<x(i)<<endln;
	}
	debug<<"]';"<<endln;


	debug <<"grad_stand=[\n";
	for(i =0; i<x.Size(); i++){
		debug<<gradientOfgFunction(i)<<endln;
	}
    debug<<"]';"<<endln;
 
*/

	Vector gradientInStandardNormalSpace = jacobian_x_u ^ gradientOfgFunction;


	if (standSens ==0) 
		standSens = new Vector(gradientInStandardNormalSpace);
	else 
	    standSens->addVector(0.0, gradientInStandardNormalSpace, 1.0); // --- save for further use.

    this->normOfGradientInStandardNormalSpace = standSens->Norm();

 // ---- ------------ -----
	
	ofstream out_a_2( "normOfGradientInStandardNormalSpace__.out");
	out_a_2 << normOfGradientInStandardNormalSpace;
    out_a_2.close();


// ==================================================================================================
	// --- perform the perturbations, since hessian is symmetric, may not need compute full of Matrix ----
 // ==================================================================================================

	double u_save;
	
	Vector column(sizeOfHessian);
	for (/*int ii=0*/; ii<sizeOfHessian; ii++) {
		gradientOfgFunction.Zero();
        opserr<<"perturbation to get hessian, ii:"<<ii<<endln;
		u_save = u(ii);

		if (fabs(u(ii))<1.0e-10) { 
			opserr<<"Warning: u("<<ii<<") is "<<u(ii)<<"!!!!"<<endln;
			u(ii)=perturbTol/100.0;
		}
		else u(ii) = u(ii)*(1.0+perturbTol);

		

		/*
		theProbabilityTransformation->set_u(u);
		result = theProbabilityTransformation->transform_u_to_x_andComputeJacobian();
			if (result < 0) {
				opserr << "Hessian::formHessianByFDM - " << endln
					<< " could not transform from u to x." << endln;
				return -1;
			}

		
		x = theProbabilityTransformation->get_x();
		*/
		result = theProbabilityTransformation->transform_u_to_x(u, x);
		if (result < 0) {
		  opserr << "Hessian::formHessianByFDM - " << endln
			 << " could not transform from u to x." << endln;
		  return -1;
		}

		if (theGFunEvaluator->setVariables() < 0) {
		  opserr << "ERROR Hessian -- error setting variables in namespace" << endln;
		  return -1;
		}

		result = theGFunEvaluator->runAnalysis();
				if (result < 0) {
					opserr << "Hessian - " << endln
						<< " could not run analysis to evaluate limit-state function. " << endln;
					return -1;
				}
		gFunctionValue = theGFunEvaluator->evaluateExpression();

		// Gradient in original space
		result = theGradGEvaluator->computeGradient(gFunctionValue);
				if (result < 0) {
					opserr << "Hessian - " << endln
						<< " could not compute gradients of the limit-state function. " << endln;
					return -1;
				}
			gradientOfgFunction = theGradGEvaluator->getGradient();
			// Gradient in standard normal space

		gradientInStandardNormalSpace = jacobian_x_u ^ gradientOfgFunction;
	
	
		
		column = (gradientInStandardNormalSpace-(*standSens))/(u(ii)-u_save);
		
        for (int i=0; i<sizeOfHessian;i++)
			(*theHessian)(i,ii)=column(i);  


		// --- recover u ---
		u(ii)=u_save;
		

	//* --debug purpose --
			
			ofstream out_a_1( "hessian_restart_.out", ios::app);
			out_a_1.precision(16);
			if (out_a_1.good()){
				for (int jj=0; jj<sizeOfHessian;jj++){
					out_a_1 << column(jj)<<endln;;
				}	
			}	
			else {
				opserr<<"Fatal: can not alloc memory for hessian_restart_.out"<<endln;
			 }
		
            out_a_1.close();
                 
     // -- debug --  */



	}  //ii loop

	// == symmetric hessian ===
	int i;
	for(i =0; i<sizeOfHessian; i++){
		for( int j=0; j<i; j++)
			(*theHessian)(i,j) = ((*theHessian)(i,j)+(*theHessian)(j,i))/2.0;
	}

	for( i =0; i<sizeOfHessian; i++){
		for( int j=i+1; j<sizeOfHessian; j++)
			(*theHessian)(i,j) = (*theHessian)(j,i);
	}


	(*theHessianInPhysicalSpace) = jacobian_u_x^(*theHessian);
    
	(*theHessianInPhysicalSpace)  =  (*theHessianInPhysicalSpace)* jacobian_u_x;
	


//	} // ----else debug purpose 




	/* ===== try to fix the singularity by discontinuity of sensitivity ===

	for ( ii=0; ii<sizeOfHessian; ii++ ) {
		int kk=0;
		int singular =-1;
		if (diagnose(ii)==-1) {
			while ((kk<10)&& (singular ==-1)){
				kk++;
				if(refineHessian(8,ii)==-1)
				{
					(*designPointUSpace)(ii) +=perturbTol;
				}
				else 
					singular =0;
			
			
			} //while 
			if (singular==-1){
				opserr<<"Fatal: can not fix hessian in column:"<<ii<<endln;
				opserr<<"You must fix it by hand before use hessian !"<<endln;
			
			}
		} //if

	}

*/
// ====== end fix of singularity  ====================================
	

	return 0;

} 

int Hessian::formHessianBySNOPT()
{
	opserr<<"Hessian::formHessianBySNOPT not done yet."<<endln;
	return -1;
}

const Matrix &
Hessian::getHessianInPhysicalSpace()
{
	
  opserr<<"Hessian in physical space is not accurate at all! and will be useless."<<endln;
  return * theHessianInPhysicalSpace;
}


int    
Hessian::sendSelf(int commitTag, Channel &theChannel)
{
  return 0;
}

int
Hessian::recvSelf(int commitTag, Channel &theChannel, 
			     FEM_ObjectBroker &theBroker)
{
  return 0;
}

const Matrix &
Hessian::getRotationMatrix(Vector alpha)
{

int dim = alpha.Size();

static Matrix R(dim,dim);

R.Zero();

for(int ii=0; ii<dim; ii++)
    R(ii,ii)=1.0;

for(int i=0; i<dim; i++)
    R(dim-1,i)=alpha(i);

Vector r0j(dim);
double tmp;

for (int j=dim-2;j>=0; j--){
 
	for(int ii=0; ii<dim; ii++){
		r0j(ii) = R(j,ii);
	}


	for (int i=dim-1; i>j; i--){
		for(int k=0; k<dim; k++){
			tmp=0;
			for(int kk=0; kk<dim; kk++)
			    tmp += R(i,kk)*r0j(kk);

            R(j,k)=R(j,k)-tmp*R(i,k);
		}
	}
	

	tmp=0;
	for(int kk=0; kk<dim; kk++)
		tmp += R(j,kk)*R(j,kk);
	tmp = pow(tmp,0.5);

	if (tmp > 1.0e-12){
		for(int k=0;k<dim;k++)
            R(j,k) /=tmp;
	}
	else {
		for(int k=0;k<dim-1;k++)
              R(j,k)=0.0;
        R(j,dim-1)=1.0;
		
	}
}
//opserr<<"R is: "<<R<<endln;
return R;
}


// --- A=(R*H*R')/|dG|, but remove last column and row --- Conte's notes middle p29
int Hessian::formReducedHessian(Vector *pDesignPt)
{
	this->formHessianByFDM(1,pDesignPt);
	Vector alpha = (*designPointUSpace)/(*designPointUSpace).Norm();
	Matrix tmp = this->getRotationMatrix(alpha);
	Matrix rot_transp = tmp;

	for(int i=0; i< tmp.noCols(); i++)
		for (int j=0; j< tmp.noRows(); j++)
			rot_transp(i,j) = tmp(j,i);
	
    Matrix A = *theHessian;




/*	if (normOfGradientInStandardNormalSpace ==0){
		ifstream out_a_3( "normOfGradientInStandardNormalSpace__.out", ios::in);
		if (out_a_3.good()){
			out_a_3 >> normOfGradientInStandardNormalSpace; //input_a_1 >> tmp;
		}
		else {
			opserr<<"can not get normOfGradientInStandardNormalSpace from normOfGradientInStandardNormalSpace__.out";
			exit(-1);
		}

		out_a_3.close();
	}


*/




//	opserr<<"normOfGradientInStandardNormalSpace"<<normOfGradientInStandardNormalSpace<<endln;
	tmp.addMatrixTripleProduct(0.0,rot_transp,A,1.0/normOfGradientInStandardNormalSpace);

	for(int z=0; z<sizeOfHessian-1; z++)
		for(int j=0; j<sizeOfHessian-1; j++)
			(*theReducedHessian)(z,j)=tmp(z,j);
	
	return 0;
}

const Matrix &
Hessian::getReducedHessian()
{
	return *theReducedHessian;  // U space
} 

int Hessian::formReducedHessian(Vector *Designpoint_X, Matrix *pHessian)
{
	// ---------------------------------------------------------
	//this->formHessianByFDM(1,pDesignPt);
	Vector x = (*Designpoint_X);


	// Transform starting point into standard normal space
	/*
	double result = theProbabilityTransformation->set_x(x);
	if (result < 0) {
		opserr << "SearchWithStepSizeAndStepDirection::doTheActualSearch() - " << endln
			<< " could not set x in the xu-transformation." << endln;
		return -1;
	}


	result = theProbabilityTransformation->transform_x_to_u();
	result = theProbabilityTransformation->transform_u_to_x_andComputeJacobian();

	if (result < 0) {
		opserr << "SHessian::formHessianByFDM - " << endln
			<< " could not transform from x to u." << endln;
		return -1;
	}
	Vector u = theProbabilityTransformation->get_u();
	*/
	Vector u;
	int result = theProbabilityTransformation->transform_x_to_u(u);
	if (result < 0) {
	  opserr << "SHessian::formHessianByFDM - " << endln
		 << " could not transform from x to u." << endln;
	  return -1;
	}

	if (designPointUSpace ==0) 
		designPointUSpace = new Vector(u.Size());

	designPointUSpace->addVector(0.0, u,1.0);

	theHessian->addMatrix(0.0, *pHessian,1.0); 


/* --- debug 1. ---
	
	for(int i=0; i<sizeOfHessian; i++){
	
	for(int j=0; j<sizeOfHessian; j++){
			
            opserr<<"(*theHessian)("<<i<<","<<j<<")"<<(*theHessian)(i,j)<<endln;
	}
	}
//  ---*/



//	Matrix jacobian_u_x = theProbabilityTransformation->getJacobian_u_x();
	Matrix jacobian_u_x(u.Size(), x.Size());

	result = theProbabilityTransformation->getJacobian_u_to_x(u, jacobian_u_x);

	(*theHessianInPhysicalSpace) = jacobian_u_x^(*theHessian);
    
	(*theHessianInPhysicalSpace)  =  (*theHessianInPhysicalSpace)* jacobian_u_x;


	// ---------------------------------------------------------

/*/ --- debug 2. ---
	for(i=0; i<sizeOfHessian; i++){
	for( int j=0; j<sizeOfHessian; j++){
			
            opserr<<"(jacobian_u_x)("<<i<<","<<j<<")"<<(jacobian_u_x)(i,j)<<endln;
	}
	}
*/
	Vector alpha = (*designPointUSpace)/(*designPointUSpace).Norm();
	Matrix tmp = this->getRotationMatrix(alpha);
	Matrix rot_transp = tmp;


/*/ --- debug 3.
    for(i=0; i<sizeOfHessian; i++){
	for( int j=0; j<sizeOfHessian; j++){
			
            opserr<<"(rot_transp)("<<i<<","<<j<<")"<<(rot_transp)(i,j)<<endln;
	}
	}

*/


	for(int i=0; i< tmp.noCols(); i++)
		for (int j=0; j< tmp.noRows(); j++)
			rot_transp(i,j) = tmp(j,i);
	
    Matrix A = *theHessian;


	if (normOfGradientInStandardNormalSpace ==0){
		ifstream out_a_3( "normOfGradientInStandardNormalSpace__.out", ios::in);
		if (out_a_3.good()){
			out_a_3 >> normOfGradientInStandardNormalSpace; //input_a_1 >> tmp;
		}
		else {
			opserr<<"can not get normOfGradientInStandardNormalSpace from normOfGradientInStandardNormalSpace__.out";
			exit(-1);
		}

		out_a_3.close();
	}




//	opserr<<"normOfGradientInStandardNormalSpace"<<normOfGradientInStandardNormalSpace<<endln;
	tmp.addMatrixTripleProduct(0.0,rot_transp,A,1.0/normOfGradientInStandardNormalSpace);


  
	for(int k=0; k<sizeOfHessian-1; k++)
		for(int j=0; j<sizeOfHessian-1; j++){
			(*theReducedHessian)(k,j)=tmp(k,j);
           // opserr<<"(*theReducedHessian)("<<i<<","<<j<<")"<<(*theReducedHessian)(i,j)<<endln;
		}
	
	return 0;
}

int Hessian::diagnose(int i)
{
	/*
sizeOfHessian=100;
     for i =1:sizeOfHessian
         j=1; endmark=0;
         while (j<=sizeOfHessian) & (endmark ==0)
             %index(i) = mean( hessian(:,i) ./ hessian(i,:)');
             if abs(hessian(i,j))>1.0e-6
                 sum =  abs(hessian(j,i)/hessian(i,j));
                 if sum > 20
                     endmark =-1;
                     [i j sum]
                     
                 end
             end
             j=j+1;
         end
          
         diff(i)=endmark;
        
         
     end
    diff'
     

	*/

	for(int j=0; j<sizeOfHessian; j++){
		if (fabs((*theHessian)(i,j))>1.0e-6) {
			double a = (*theHessian)(j,i)/(*theHessian)(i,j);
		//	if ((fabs(a)>5) || ((a<-1.e-12) &&(a>-1.0))){
			if (fabs(a)>20){
			   return -1;
			}
		} 
	}
		
	

	return 0;

}

int Hessian::refineHessian(int time, int colOfHessian)
{
	Vector column(sizeOfHessian);
	
	double tol = perturbTol;
	int ii = colOfHessian;
	int success =-1;
	int kk=0;
	static Vector grad_1(sizeOfHessian);
	static Vector grad_2(sizeOfHessian);
	static Vector grad_3(sizeOfHessian);
	
	Vector u = *designPointUSpace;

	/*
	theProbabilityTransformation->set_u(u);
	int result = theProbabilityTransformation->transform_u_to_x_andComputeJacobian();
	if (result < 0) {
		opserr << "Hessian::formHessianByFDM - " << endln
			<< " could not transform from u to x." << endln;
		return -1;
	}
	Matrix jacobian_x_u = theProbabilityTransformation->getJacobian_x_u();
	*/
	Matrix jacobian_x_u;
	int result = theProbabilityTransformation->getJacobian_x_to_u(jacobian_x_u);
	if (result < 0) {
		opserr << "Hessian::formHessianByFDM - " << endln
			<< " could not transform from u to x." << endln;
		return -1;
	}

	double u_1, u_2, u_3;

    u_1 =u(ii);
	grad_1.addVector(0.0, *standSens,1.0);
	
	//u_2 = u(ii)*tol;
	//grad_2=g(u_2);

	Vector gradientOfgFunction(sizeOfHessian );

    double u_save = u(ii);

	if (fabs(u(ii))<1.0e-10) { 
		opserr<<"Warning: u("<<ii<<") is "<<u(ii)<<"!!!!"<<endln;
	    u(ii)=-tol/100.0;
	}
	else u(ii) = u(ii)*(1.0-tol);
	
	u_2=u(ii);
		
    opserr<<"re-perturbation to get hessian, ii:"<<ii<<". u_2:"<<u_2<<endln;
	

    /*
	theProbabilityTransformation->set_u(u);
	result = theProbabilityTransformation->transform_u_to_x_andComputeJacobian();
		if (result < 0) {
			opserr << "Hessian::formHessianByFDM - " << endln
				<< " could not transform from u to x." << endln;
			return -1;
		}

		
	Vector x = theProbabilityTransformation->get_x();
    */
    Vector x;
    theProbabilityTransformation->transform_u_to_x(u, x);
    if (result < 0) {
      opserr << "Hessian::formHessianByFDM - " << endln
	     << " could not transform from u to x." << endln;
      return -1;
    }

	if (theGFunEvaluator->setVariables() < 0) {
	  opserr << "ERROR Hessian -- error setting variables in namespace" << endln;
	  return -1;
	}
		result = theGFunEvaluator->runAnalysis();
				if (result < 0) {
					opserr << "Hessian - " << endln
						<< " could not run analysis to evaluate limit-state function. " << endln;
					return -1;
				}
		double	gFunctionValue = theGFunEvaluator->evaluateExpression();

        // Gradient in original space
		result = theGradGEvaluator->computeGradient(gFunctionValue);
				if (result < 0) {
					opserr << "Hessian - " << endln
						<< " could not compute gradients of the limit-state function. " << endln;
					return -1;
				}
			gradientOfgFunction = theGradGEvaluator->getGradient();
			// Gradient in standard normal space

	
	
	 grad_2 = jacobian_x_u ^ gradientOfgFunction;
	
	
	column = (grad_2-grad_1)/(u_2-u_1);
		
    for (int i=0; i<sizeOfHessian;i++)
		(*theHessian)(i,ii)=column(i);  

    u(ii)=u_save;

	if (diagnose(ii)==0) success =0;

	while( (success<0) && (kk<time)){
	    kk++;
		tol /=2.0;
		opserr<<"re-perturbation to get hessian,kk:"<<kk<<".  tol:"<<tol<<endln;
        //==u_3=u(ii)*(1-tol);
		//==grad_3=g(u_3);
		//==column=(grad_3-grad_1)/u_3-u_1;

		if (fabs(u(ii))<1.0e-10) { 
			opserr<<"Warning: u("<<ii<<") is "<<u(ii)<<"!!!!"<<endln;
			u(ii)=-tol/100.0;
		}
		else u(ii) = u(ii)*(1.0-tol);
		
		u_3=u(ii);
		
		opserr<<"re-perturbation to get hessian,u_3:" <<u_3<<endln;
	
	/*
		theProbabilityTransformation->set_u(u);
		result = theProbabilityTransformation->transform_u_to_x_andComputeJacobian();
			if (result < 0) {
				opserr << "Hessian::formHessianByFDM - " << endln
					<< " could not transform from u to x." << endln;
				return -1;
			}

			
		Vector x = theProbabilityTransformation->get_x();
	*/
		Vector x; 
		result = theProbabilityTransformation->transform_u_to_x(u,x);
		if (result < 0) {
		  opserr << "Hessian::formHessianByFDM - " << endln
			 << " could not transform from u to x." << endln;
		  return -1;
		}

	if (theGFunEvaluator->setVariables() < 0) {
	  opserr << "ERROR Hessian -- error setting variables in namespace" << endln;
	  return -1;
	}

			result = theGFunEvaluator->runAnalysis();
					if (result < 0) {
						opserr << "Hessian - " << endln
							<< " could not run analysis to evaluate limit-state function. " << endln;
						return -1;
					}
			gFunctionValue = theGFunEvaluator->evaluateExpression();

            // Gradient in original space
			result = theGradGEvaluator->computeGradient(gFunctionValue);
					if (result < 0) {
						opserr << "Hessian - " << endln
							<< " could not compute gradients of the limit-state function. " << endln;
						return -1;
					}
				gradientOfgFunction = theGradGEvaluator->getGradient();
				// Gradient in standard normal space

		 
		grad_3 = jacobian_x_u ^ gradientOfgFunction;
	
	    u(ii)=u_save;

	    column = (grad_3-grad_1)/(u_3-u_1);
		
		for (int i=0; i<sizeOfHessian;i++)
			(*theHessian)(i,ii)=column(i);  


		if (diagnose(ii)==0)
			success =0;
		else{
			column=(grad_3-grad_2)/u_3-u_2;
		
			if (diagnose(ii)==0)
				success =0;
			else{
				u_2=u_3;
				grad_2=grad_3;
			}
		}
	} //while 
	
	return success;
}

double Hessian::getNormOfGradientU()
{
	return this->normOfGradientInStandardNormalSpace;
}
