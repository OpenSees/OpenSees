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
 

#include <DP_RSM_Sim.h>
#include <MatrixOperations.h>
#include <ostream>
#include <math.h>


#ifndef _WIN32
void itoa(int, char *, int base);
#endif

//////////////////////////////////////////////////////////////////////
// command: runDP_RSM_SimTimeInvariantAnalysis -designPt dp.out  -output results.out  -ndir $n <-experimentalPointRule Uniform -gridInfo {-1  minY  maxY nPts 0  minY  maxY nPts0 1 minX1  maxX1 nPts1 2 minX2  maxX2 nPts2 ...}> 
//  -saveHessian hession.out <-surfaceDesign UnivariateDecomposition -simulation ImportanceSampling -tarCOV 0.1 -numSimulation 100000>

//////////////////////////////////////////////////////////////////////

DP_RSM_Sim::DP_RSM_Sim(ReliabilityDomain *passedReliabilityDomain,
					FunctionEvaluator *passedGFunEvaluator,
					ProbabilityTransformation *passedProbabilityTransformation,
					char *passedOutputFileName,
					GradientEvaluator * passedGradGEvaluator, Vector * pDesignPt, int numAxis, 
					char * typeExpPtRule,Tcl_Interp *passedTclInterp, 
					Matrix * passedHessian, char * passedHessianFile, char * typeSurfaceDesign, 
					char * typeRespSurfaceSimulation, Vector * gridInfo,
					RandomNumberGenerator * pRandomNumberGenerator,
					double pTargetCOV,
					int pNumberOfSimulations)
					
					:ReliabilityAnalysis()
{
    theReliabilityDomain = passedReliabilityDomain;
 	theProbabilityTransformation = passedProbabilityTransformation;
	theGFunEvaluator = passedGFunEvaluator;	
    theGradGEvaluator = passedGradGEvaluator;
	theRandomNumberGenerator =  pRandomNumberGenerator;


	
	HessianFileName = 0;
	if (passedHessianFile !=0) {
		HessianFileName = new char[30];
		strcpy(HessianFileName, passedHessianFile);
	
	}

	HessianMatrix =0;
	numOfPrincipalAxes = numAxis;

	strcpy(outputFileName,passedOutputFileName);
	
	rotation =0;

	if (pDesignPt ==0) {opserr<<"no designpt."<<endln; exit(-1);};
	theDesignPtXSpace = new Vector (*pDesignPt);

	/*
	theProbabilityTransformation->set_x(*theDesignPtXSpace);
	theProbabilityTransformation->transform_x_to_u();
	int numRV = theReliabilityDomain->getNumberOfRandomVariables();
	theDesignPoint = new Vector(numRV);
	theDesignPoint->addVector(0.0,theProbabilityTransformation->get_u(),1.0);
	*/

	int numRV = theReliabilityDomain->getNumberOfRandomVariables();
	theDesignPoint = new Vector(numRV);	
	theProbabilityTransformation->transform_x_to_u(*theDesignPoint);
	
	if (strcmp(typeExpPtRule, "Uniform") ==0 ) {
		theExpPtRule= new UniformExperimentalPointRule1D();
	}
	else {
	    opserr<<"ExperimentalPointRule not implemented yet"; exit(-1);
	}



	int numOfGridPlanes = (numOfPrincipalAxes+1)*numOfPrincipalAxes/2;
	theGridPlanes = new GridPlane * [numOfGridPlanes]; // uplimit 20 for now.refer arrayOfTaggedObject.. theComponents = new TaggedObject *[size]
	for(int ii = 0;ii< numOfGridPlanes; ii++) theGridPlanes[ii]=0;


	

	//theHessian = 0;

	if (passedHessian !=0){
		//HessianMatrix = new Matrix( *passedHessian);
		//theHessian = new Hessian(numRV,theReliabilityDomain,
            //theProbabilityTransformation,theGFunEvaluator,theGradGEvaluator, 1.0e-5);
	
		//theHessian->formReducedHessian(theDesignPtXSpace, passedHessian); 

	}

	this->theInterp=passedTclInterp;

 
// ----
	if (strcmp(typeSurfaceDesign, "UnivariateDecomposition")==0){
		theSurfaceDesign = new UnivariateDecomposition(numOfPrincipalAxes+1, 0, false); // note here numOfPrincipalAxes+1 since dp direction is 1st axis
		// need to set pointer to pointer arry of principalAxes later
	}
	else if (strcmp(typeSurfaceDesign, "BivariateDecomposition")==0){
		//opserr<<"typeSurfaceDesign not implemented yet "; 
		theSurfaceDesign = new BivariateDecomposition(numOfPrincipalAxes+1, 0, 0 , false); // note here numOfPrincipalAxes+1 since dp direction is 1st axis

		// need to set (1) pointer to pointer arry of principalAxes (2) pointer to pointer arry of gridPlanes later
	}
	else {
		opserr<<"unknown typeSurfaceDesign "; exit(-1);
	}
// -----
	theTargetCOV = pTargetCOV;
    theNumberOfSimulations = pNumberOfSimulations;	
	
	if (strcmp(typeRespSurfaceSimulation, "ImportanceSampling")==0){
		theRespSurfaceSimulation = new RespSurfaceSimulation(
			            theReliabilityDomain,
						theDesignPoint,
						outputFileName,
						theSurfaceDesign,
						theRandomNumberGenerator,
                        theTargetCOV,
						theNumberOfSimulations);
	}
	else {
		opserr<<"unknown sampling "; exit(-1);
	}	
// -----	
	thePrincipalAxes = new PrincipalAxis * [numOfPrincipalAxes+1]; // dp direction + all principalAxis
	for(int ii = 0;ii< numOfPrincipalAxes+1; ii++) thePrincipalAxes[ii]=0;


	this->setGridInfo( gridInfo );  // create all principalAxis[i], in which set experimentalPointRule with Grid info. default: [-1.0:11:1.0]
	


  
}

DP_RSM_Sim::~DP_RSM_Sim()
{
   // need to delete every memory alloced by 'new'.... ??????
}

int DP_RSM_Sim::analyze()
{
  //1  form Hessian;
  //2  get  eigenvectors & eigenvalues
  //3  order eigenVectors
  //4. grid values computation
  //5. fit curve
  //6. simulation


	// --1.1. create grid planes, set grid info ----

	int numOfGridPlanes = (numOfPrincipalAxes+1)*numOfPrincipalAxes/2;
 
    			
	ofstream output( outputFileName, ios::out );



	//1.  --- compute and form Hessian by FDM; compute A, Rotation R, A=R*H*R'/|dG|, but one dim less.

	int numRV = theReliabilityDomain->getNumberOfRandomVariables();
  	
	//if (theHessian ==0) {
	//	theHessian = new Hessian(numRV,theReliabilityDomain,theProbabilityTransformation,
    //                    theGFunEvaluator,theGradGEvaluator, 1.0e-3);
	//	theHessian->formReducedHessian(theDesignPtXSpace);
	//}

	if (HessianFileName !=0){  //recorder Hessian

		Matrix tmpHessian(numRV, numRV);
		//tmpHessian.addMatrix(0.0, theHessian->getHessianApproximation(), 1.0);

		ofstream output(HessianFileName);
		
		int ii=0;
		int jj=0;
	
		for(ii=0; ii<numRV; ii++){ 
			for (jj=0; jj<numRV;jj++){
				output << tmpHessian(jj,ii)<<endln ;
			}	
		
		}
		output.close();
		
	}





    if (HessianMatrix !=0) delete HessianMatrix;
	if (rotation !=0) delete rotation;
	 
	//HessianMatrix = new Matrix(theHessian->getReducedHessian()); //A matrix, one dim less the nrv
    HessianMatrix = new Matrix(1,1);
	
	Vector alpha =  (*theDesignPoint)/(*theDesignPoint).Norm();


	//rotation =  new Matrix(theHessian->getRotationMatrix(alpha)); 

    MatrixOperations * theMatrixOperations = new MatrixOperations(*HessianMatrix);
	
	// 2 & 3. -- eigen analysis, insert PrincipalAxis into PrincipalAxes array according to magnitute


  
/*  
  
  ===Before order:
  pointer:    0    1    2    3    4    5    6    7    8            
  #eigenVect  1    2    3    4    5    6    7    8    9 =numOfPrincipalAxes
			 __________________________________________________ 
			|    |    |    |    |    |    |    |    |    |    |
			|    |    |    |    |    |    |    |    |    |   \|
			|____|____|____|____|____|____|____|____|____|____|\    
                                                                _\| 
Hessian size = numRV-1                                           not used





  ===After order change to ===>

  pointer:    0    1    2    3    4    5    6    7    8    9                 size = numOfPrincipalAxes+1          
  #eigenVect  0    1    2    3    4    5    6    7    8    9 =numOfPrincipalAxes
			 _________________________________________________ 
			|    |    |    |    |    |    |    |    |    |    |
			|    |    |    |    |    |    |    |    |    |    |
			|____|____|____|____|____|____|____|____|____|____|



*/
	if (numOfPrincipalAxes*2>=numRV-1){
		theMatrixOperations->performEigenAnalysis(1,numRV-1);

		// save into principalplanes
		double eigenValue;
		Vector eigenVector(numRV);

		Vector tmp(numRV-1);
		Vector tmp1(numRV);
		tmp1(numRV-1)=0.0;

		int i;
		for(i=0; i<numOfPrincipalAxes; i++){
			eigenValue = theMatrixOperations->getEigenvalue(i+1);
			tmp.addVector(0.0, theMatrixOperations->getEigenvector(i+1),1.0);
			for (int j =0; j<numRV-1; j++)
				tmp1(j)=tmp(j);

			eigenVector = (*rotation) ^ tmp1;
			
			thePrincipalAxes[i]->setCurvature(eigenValue);
			thePrincipalAxes[i]->setAxisDirection(&eigenVector);

		}
/*			
  pointer:    0    1    2    3    4    5    6    7    8               
  #eigenVect  1    2    3    4    5    6    7    8    9  =numRV-1                    numRV = 10
			 _____________________________________________                  
			|    |    |    |    |    |    |    |    |    |                 
			|    |    |    |    |    |    |    |    |    |     
			|____|____|____|____|____|____|____|____|____|         
                                   | 
                                  \|/
                         numOfPrincipalAxes =5
  */

		int insertPlace = 0;			
		for( i=numRV-2; i>0; i--){

			eigenValue = theMatrixOperations->getEigenvalue(i+1);
			// -- comp abs of eigenValue


			while(insertPlace < numOfPrincipalAxes){
				if ( fabs(thePrincipalAxes[insertPlace]->getCurvature())>=fabs(eigenValue)){
					insertPlace++;
				}
				else 
					break;

			}

			if (insertPlace != numOfPrincipalAxes){
				for( int k=numOfPrincipalAxes-1; k>insertPlace; k--){
					thePrincipalAxes[k]->copyValues( thePrincipalAxes[k-1]);
				}

				tmp.addVector(0.0, theMatrixOperations->getEigenvector(i+1),1.0);
				for (int j =0; j<numRV-1; j++)
					tmp1(j)=tmp(j);
				tmp1(numRV-1)=0.0;

				eigenVector = (*rotation) ^ tmp1;
			
				thePrincipalAxes[insertPlace]->setCurvature(eigenValue);
				thePrincipalAxes[insertPlace]->setAxisDirection(&eigenVector);
				thePrincipalAxes[insertPlace]->cleanValuesOnAxis();
				insertPlace++;
			}
			else if (insertPlace == numOfPrincipalAxes){ i=0; } // to exit for loop
		}  //for

	}
	else{

/*			
  pointer:    0    1    2    3    4    5    6    7    8               
  #eigenVect  1    2    3    4    5    6    7    8    9  =numRV-1                    numRV = 10
			 _____________________________________________                  
			|    |    |    |    |    |    |    |    |    |                 
			|    |    |    |    |    |    |    |    |    |     
			|____|____|____|____|____|____|____|____|____|         
                         | 
                        \|/
                 numOfPrincipalAxes =3
  */


		theMatrixOperations->performEigenAnalysis(1,numOfPrincipalAxes);
		
		// save into principalplanes
		double eigenValue;
		Vector eigenVector(numRV);

		Vector tmp(numRV-1);
		Vector tmp1(numRV);
		tmp1(numRV-1)=0.0;

		int i;
		for(i=0; i<numOfPrincipalAxes; i++){
			eigenValue = theMatrixOperations->getEigenvalue(i+1);
			tmp.addVector(0.0, theMatrixOperations->getEigenvector(i+1),1.0);
			for (int j =0; j<numRV-1; j++)
				tmp1(j)=tmp(j);
			tmp1(numRV-1)=0.0;

			eigenVector = (*rotation) ^ tmp1;

			thePrincipalAxes[i]->setCurvature(eigenValue);
			thePrincipalAxes[i]->setAxisDirection(&eigenVector);

		}

		theMatrixOperations->performEigenAnalysis(numRV-numOfPrincipalAxes,numRV-1);

		int insertPlace = 0;
		for( i=numRV-2; i>=numRV-numOfPrincipalAxes-1; i--){

			eigenValue = theMatrixOperations->getEigenvalue(i+1);
			// -- comp abs of eigenValue
			
			
			while(insertPlace < numOfPrincipalAxes){
				if ( fabs(thePrincipalAxes[insertPlace]->getCurvature())>=fabs(eigenValue)){
					insertPlace++;
				}
				else 
					break;

			} //while
			

			if (insertPlace != numOfPrincipalAxes){
				for( int k=numOfPrincipalAxes-1; k>insertPlace; k--){
					thePrincipalAxes[k]->copyValues( thePrincipalAxes[k-1]);
				}

				tmp.addVector(0.0, theMatrixOperations->getEigenvector(i+1),1.0);
				for (int j =0; j<numRV-1; j++)
					tmp1(j)=tmp(j);
				tmp1(numRV-1)=0.0;

				eigenVector = (*rotation) ^ tmp1;
							
				thePrincipalAxes[insertPlace]->setCurvature(eigenValue);
				thePrincipalAxes[insertPlace]->setAxisDirection(&eigenVector);
				thePrincipalAxes[insertPlace]->cleanValuesOnAxis();
				
				insertPlace++;
				// note gridInfo does not change..
		
			}
			else {i=numRV-numOfPrincipalAxes-1; }  // stop for loop
	

		}  //for
	} //else if (numOfGridPlanes*2>=numRV-1)

	
	
	


/*  
  
  ===Before order:
  pointer:    0    1    2    3    4    5    6    7    8            
  #eigenVect  1    2    3    4    5    6    7    8    9 =numOfPrincipalAxes
			 __________________________________________________ 
			|    |    |    |    |    |    |    |    |    |    |
			|    |    |    |    |    |    |    |    |    |   \|
			|____|____|____|____|____|____|____|____|____|____|\    
                                                                _\| 
Hessian size = numRV-1                                           not used





  ===After order change to ===>

  pointer:    0    1    2    3    4    5    6    7    8    9                    size = numOfPrincipalAxes+1          
  #eigenVect  0    1    2    3    4    5    6    7    8    9 =numOfPrincipalAxes
			 _________________________________________________ 
			|    |    |    |    |    |    |    |    |    |    |
			|    |    |    |    |    |    |    |    |    |    |
			|____|____|____|____|____|____|____|____|____|____|



*/

	int k;
	
	for (k =numOfPrincipalAxes; k>0; k--){
		thePrincipalAxes[k]->copyValues( thePrincipalAxes[k-1]);
	
	} //for k

	thePrincipalAxes[0]->setCurvature(1.0e20);  // infinite for dp direction
    thePrincipalAxes[0]->setAxisDirection(&alpha);
	thePrincipalAxes[0]->cleanValuesOnAxis();


	ofstream file_coeff("axis_recorder.out", ios::out);



	for ( k =0; k<=numOfPrincipalAxes; k++){

		double a = thePrincipalAxes[k]->getCurvature();
		opserr<<"the "<<k<<"th curvature:"<<a<<endln;
		file_coeff<<"the "<<k<<"th curvature:"<<a<<endln;

		
		Vector * aa = thePrincipalAxes[k]->getAxisDirection();
		opserr<<"the "<<k<<"th direction:\n";
		file_coeff<<"the "<<k<<"th direction:\n";

		for( int m =0; m< aa->Size(); m++){
			opserr<<(*aa)(m)<<endln;
			file_coeff<<(*aa)(m)<<endln;
	
		}

	} //for k




// ===============================================================
//4  get grid values


	// 4.1 for each axis, compute the values
 
	char cmd[50]="remove sensitivityAlgorithm";
	Tcl_Eval(theInterp, cmd);


		// --- create gridPlanes ----

	int ii;
	for (ii=0 ; ii<numOfGridPlanes; ii++){ 

		//	int numOfGrid = thePrincipalAxes[ii]->getExperimentalPointRule()->getNumberOfPoints();
			int i = this->getNumOfAxis(1,ii+1);
			int j = this->getNumOfAxis(2,ii+1);

			theGridPlanes[ii]=new GridPlane(ii+1, theDesignPoint, 
									   thePrincipalAxes[i], thePrincipalAxes[j], rotation,
									   theProbabilityTransformation,
									   theGFunEvaluator);
	} //for ii



	double gFunValue=0;

	for(  ii=0; ii<numOfPrincipalAxes+1; ii++) {  
		opserr<<"========"<<ii<<"th axis: =======\n";
		file_coeff<<"========"<<ii<<"th axis: =======\n";

		if (ii==0){  // dp direction
			int numOfGrid = thePrincipalAxes[ii]->getExperimentalPointRule()->getNumberOfPoints();  // 1st direction: dp, 2nd: eigenvector

			for(int j=0; j<numOfGrid; j++){
                
				opserr<<"j is: "<<j<<" \n"; 
				file_coeff<<"j is: "<<j<<" \n"; 
				
				//Vector coord = theGridPlanes[ii]->getPointCoordinate(j,0);
				double coord = thePrincipalAxes[ii]->getExperimentalPointRule()->getPointCoordinate(j);
				gFunValue = theGridPlanes[ii]->getValueOnGrid(coord, 0.0);

				// --- deal with FE divergence cases, temporary ---

				if(theGridPlanes[ii]->isFEConverged()){
					thePrincipalAxes[ii]->setValueOnAxis(j,gFunValue);
				}
				else{  // divergence
				
					if ((j==0) ||(j==1)){
						opserr<<"Fatal: the range of simulation in design point direction is inside the failure domain, you may decrease the grid bound.."<<endln;
						exit(-1);
					}
					
					else{ // interpolation
						
						
						
						double x = coord;
						//Vector coord0 = theGridPlanes[ii]->getPointCoordinate(j-2,0);
						double x0 = thePrincipalAxes[ii]->getExperimentalPointRule()->getPointCoordinate(j-2);
							//theGridPlanes[ii]->getPointCoordinate(j-2,0)(0);
						double g0 = thePrincipalAxes[ii]->getValueOnAxis(j-2);

						//Vector coord1 = theGridPlanes[ii]->getPointCoordinate(j-1,0);
						double x1 = thePrincipalAxes[ii]->getExperimentalPointRule()->getPointCoordinate(j-1);
							//theGridPlanes[ii]->getPointCoordinate(j-1,0)(0);
						double g1 = thePrincipalAxes[ii]->getValueOnAxis(j-1);

						gFunValue = (g1-g0)/(x1-x0)*(x-x1)+g1;

						if (gFunValue>=0) gFunValue = -1.0e-14;

						thePrincipalAxes[ii]->setValueOnAxis(j,gFunValue);
					}
				
				
				
				}

				opserr<<"jth value  is: "<<gFunValue<<" \n";
				file_coeff<<"jth value  is: "<<gFunValue<<" \n";
			
			}

		} // if (ii==0)  dp direction

		else{  // ii =1,2,3..i ...n,  the ith eigenvector direction. 
			int numOfGrid = thePrincipalAxes[ii]->getExperimentalPointRule()->getNumberOfPoints(); 
			
			// find zero point. 

			int startPt = thePrincipalAxes[ii]->getExperimentalPointRule()->getPointClosestToOrigin();


			int j;
			for(j=startPt; j<numOfGrid; j++){
				opserr<<"j is: "<<j<<" \n";
				file_coeff<<"j is: "<<j<<" \n";

				double coord =  thePrincipalAxes[ii]->getExperimentalPointRule()->getPointCoordinate(j); //theGridPlanes[ii-1]->getPointCoordinate(0,j);
				gFunValue = theGridPlanes[ii-1]->getValueOnGrid(0, coord);
//				thePrincipalAxes[ii]->setValueOnAxis(j,gFunValue);

				// --- deal with FE divergence cases, temporary ---

				if(theGridPlanes[ii-1]->isFEConverged()){
					thePrincipalAxes[ii]->setValueOnAxis(j,gFunValue);
				}
				else{  // divergence
				
					if ((j==startPt) ||(j==startPt+1)){
						opserr<<"Fatal: the range of simulation in design point direction is inside the failure domain, you may change the grid bound.."<<endln;
						exit(-1);
					}
					
					else{ // interpolation
						
						double y = coord;
						
						//Vector coord0 = theGridPlanes[ii-1]->getPointCoordinate(0,j-2); 
						double y0 = thePrincipalAxes[ii]->getExperimentalPointRule()->getPointCoordinate(j-2);
						double g0 = thePrincipalAxes[ii]->getValueOnAxis(j-2);

						//Vector coord1 = theGridPlanes[ii-1]->getPointCoordinate(0,j-1); 
						double y1 = thePrincipalAxes[ii]->getExperimentalPointRule()->getPointCoordinate(j-1);
						double g1 = thePrincipalAxes[ii]->getValueOnAxis(j-1);

						gFunValue = (g1-g0)/(y1-y0)*(y-y1)+g1;

						if (gFunValue>=0) gFunValue = -1.0e-14;

						thePrincipalAxes[ii]->setValueOnAxis(j,gFunValue);
					} //else
				} //else
				opserr<<"jth value  is: "<<gFunValue<<" \n";
				file_coeff<<"jth value  is: "<<gFunValue<<" \n";
			
			} //for (int j=startPt; j<numOfGrid; j++)

			for( j=startPt-1; j>=0; j--){
				opserr<<"j is: "<<j<<" \n";
				file_coeff<<"j is: "<<j<<" \n";

				double coord =  thePrincipalAxes[ii]->getExperimentalPointRule()->getPointCoordinate(j); //theGridPlanes[ii-1]->getPointCoordinate(0,j);
				gFunValue = theGridPlanes[ii-1]->getValueOnGrid(0, coord);


				// --- deal with FE divergence cases, temporary ---

				if(theGridPlanes[ii-1]->isFEConverged()){
					thePrincipalAxes[ii]->setValueOnAxis(j,gFunValue);
				}
				else{  // divergence
				
					if (j==startPt-1) {
						opserr<<"Fatal: the range of simulation in design point direction is inside the failure domain, you may change the grid bound.."<<endln;
						exit(-1);
					}
					
					else{ // interpolation
						
						double y = coord;
						
						//Vector coord0 = theGridPlanes[ii-1]->getPointCoordinate(0,j-2); 
						double y0 = thePrincipalAxes[ii]->getExperimentalPointRule()->getPointCoordinate(j+2);
						double g0 = thePrincipalAxes[ii]->getValueOnAxis(j+2);

						//Vector coord1 = theGridPlanes[ii-1]->getPointCoordinate(0,j-1); 
						double y1 = thePrincipalAxes[ii]->getExperimentalPointRule()->getPointCoordinate(j+1);
						double g1 = thePrincipalAxes[ii]->getValueOnAxis(j+1);

						gFunValue = (g1-g0)/(y1-y0)*(y-y1)+g1;

						if (gFunValue>=0) gFunValue = -1.0e-14;

						thePrincipalAxes[ii]->setValueOnAxis(j,gFunValue);
					} //else
				} //else
				opserr<<"jth value  is: "<<gFunValue<<" \n";
				file_coeff<<"jth value  is: "<<gFunValue<<" \n";
			
			} //for (int j=startPt; j>=0; j--)


		} //else // ii =1,2,3..i ...n,  the ith eigenvector direction.
	


	} // for  each axis..



 

	if (strcmp(theSurfaceDesign->getType(), "UnivariateDecomposition")==0) {
		theSurfaceDesign->setPincipalAxesPtr(thePrincipalAxes);

	}

	else if (strcmp(theSurfaceDesign->getType(), "BivariateDecomposition")==0) {
		theSurfaceDesign->setPincipalAxesPtr(thePrincipalAxes);
		theSurfaceDesign->setGridPlanesPtr(theGridPlanes);

//=================================================================



	   // 4.2 values on grid planes

		char fileName[30];
		strcpy(fileName, outputFileName);
		char seps[]   = ".";
		char * token;

    	ofstream output_v( "gridValues.out", ios::out );

		for(  ii=0; ii<numOfGridPlanes; ii++) {  //for each plane
		
			output_v <<"plane is ii="<<ii<<"\n\n======";

			token = strtok( outputFileName, seps );
			strcpy(fileName, token);
			strcat(fileName, "-vis");
    
			token = strtok( NULL, seps );
			
			char extension[5]=".";
			if (token !=NULL)    
				strcat(extension, token);
			else
				strcat(extension, "out");
   			char theStr[3];

			//itoa(ii+1,theStr,10);
			sprintf(theStr,"%d",ii+1);
			strcat(fileName, theStr);
			strcat(fileName, extension);
    			
			ofstream output_vis( fileName, ios::out );



			double x, y;
			int numOfGridXPt = theGridPlanes[ii]->getNumOfGridXPt();
			int numOfGridYPt = theGridPlanes[ii]->getNumOfGridYPt();

			PrincipalAxis * axis_1 = theGridPlanes[ii]->getAxisPtr(1);
			PrincipalAxis * axis_2 = theGridPlanes[ii]->getAxisPtr(2);

			ExperimentalPointRule1D * rule_1 = axis_1->getExperimentalPointRule();
			ExperimentalPointRule1D * rule_2 = axis_2->getExperimentalPointRule();
			

			
			int jj;
			for(jj=0; jj< numOfGridXPt; jj++){
				for(int kk=0; kk< numOfGridYPt; kk++){

					x = rule_1->getPointCoordinate(jj);
					y = rule_2->getPointCoordinate(kk);

					if (fabs(x)<1.0e-14) {
						gFunValue = axis_2->getValueOnAxis(kk);
						
					}
					else if (fabs(y)<1.0e-14) {
						gFunValue = axis_1->getValueOnAxis(jj);
					}
					else{
						gFunValue = theGridPlanes[ii]->getValueOnGrid(x,y);

						if (! theGridPlanes[ii]->isFEConverged())
						       gFunValue =-1.0e20;  // will correct it later..
					}
						// write x,y valueOfG
						
					theGridPlanes[ii]->setGridValue(jj,kk,gFunValue);
						
					output_v <<"jj="<<jj<<", kk="<<kk<<",x="<<x<<", y="<<y<<", value = "<< gFunValue<<endln;			
						
					opserr<<"-----------------Before correction---------------\n";
					opserr<<"Plain:"<<ii<<": "<<x<<","<<y<<","<<gFunValue<<"\n";	
						


				}  // for kk
			}  //for jj

			// --- write for visulaze --
			for( jj=0; jj< numOfGridXPt; jj++){
				for(int kk=0; kk< numOfGridYPt; kk++){
					double x = rule_1->getPointCoordinate(jj);
					double y = rule_2->getPointCoordinate(kk);
					double z = theGridPlanes[ii]->getSavedValueOnGrid(jj,kk);
					output_vis<<x<<"   "<<y<<"   "<<z<<"\n";
					opserr<<"Before correction\n";
					opserr<<"Plain:"<<ii<<": "<<x<<","<<y<<","<<z<<"\n";	
				
				}
			} //jj


       //  ======== correction for the nonconvergence cases =======
/*	m=3; 
	n=6;
	startPt=1;
	for sum =2:m+n
	   if(sum-1<=m)
		   for i = startPt: min(sum-1,n) 
			   [ i, sum-i]
		   end
	   else
			startPt = startPt+1;
			for i = startPt: min(sum-1,n)
			   [ i, sum-i]
			end 
		end
	end

 */
		int startPtX = axis_1->getExperimentalPointRule()->getPointClosestToOrigin(); // note this startPtX begin from 0
        int startPtY = axis_2->getExperimentalPointRule()->getPointClosestToOrigin();// note this startPtY begin from 0
		
		
		double newValue =0.0;

		// --- 3rd (-1,-1) ---
		int m = startPtX;
		int n = startPtY;

        opserr<<"======================================="<<endln;
		int  startPt=1;
		int sum;
		for (sum =2; sum <=m+n; sum++){
			if(sum-1<=m){
				int min = sum-1;
				if (sum-1>n) min = n;

				for(int i_y = startPt; i_y<= min; i_y++){ 

				   opserr<<"checking: ["<<  -(sum-i_y)<<","<<-i_y<<"]"<<endln;
				   if (fabs(theGridPlanes[ii]->getSavedValueOnGrid(startPtX-(sum-i_y), startPtY-i_y))>1.0e19) {

					   newValue = FEDivergenceCorrection(ii, startPtX, startPtY, -(sum-i_y), -i_y, m,n, 3); // (startPtX, startPtY, incr_i_x, incr_i_y, maxincr_i_x, maxincr_i_y, first)
					   theGridPlanes[ii]->setGridValue(startPtX-(sum-i_y), startPtY-i_y, newValue);
					   opserr<<"correcting: ["<<  -(sum-i_y)<<","<<-i_y<<"] to"<<newValue<<endln;
					   output_v <<"correction::::: jj="<<startPtX-(sum-i_y)<<", kk="<<startPtY-i_y<<", value = "<< newValue<<endln;
				   }

				   //opserr<<"correcting: ["<<  -(sum-i_y)<<","<<-i_y<<"]"<<endln;
				}
			}//if
			else{
				startPt = startPt+1;
				int min = sum-1;
				if (sum-1>n) min = n;
				
				for(int i_y = startPt; i_y<= min; i_y++){ 
				   opserr<<"checking: ["<<  -(sum-i_y)<<","<<-i_y<<"]"<<endln;
				   if (fabs(theGridPlanes[ii]->getSavedValueOnGrid(startPtX-(sum-i_y), startPtY-i_y))>1.0e19) {

					   newValue = FEDivergenceCorrection(ii, startPtX, startPtY, -(sum-i_y), -i_y, m,n, 3); // (startPtX, startPtY, incr_i_x, incr_i_y, maxincr_i_x, maxincr_i_y, first)
					   theGridPlanes[ii]->setGridValue(startPtX-(sum-i_y), startPtY-i_y, newValue);
					   opserr<<"correcting: ["<<  -(sum-i_y)<<","<<-i_y<<"] to"<<newValue<<endln;
					   output_v <<"correction::::: jj="<<startPtX-(sum-i_y)<<", kk="<<startPtY-i_y<<", value = "<< newValue<<endln;
				   }


				   //opserr<<"correcting: ["<<  -(sum-i_y)<<","<<-i_y<<"]"<<endln;
				}
			}//else
			
		}// for sum

		// --- 4th (1,-1)  ---

		m = numOfGridXPt - startPtX-1;
		n = startPtY;

        opserr<<"======================================="<<endln;
		 startPt=1;
		for ( sum =2; sum <=m+n; sum++){
			if(sum-1<=m){
				int min = sum-1;
				if (sum-1>n) min = n;

				for(int i_y = startPt; i_y<= min; i_y++){ 

				   opserr<<"checking: ["<<  (sum-i_y)<<","<<-i_y<<"]"<<endln;
				   if (fabs(theGridPlanes[ii]->getSavedValueOnGrid(startPtX+(sum-i_y), startPtY-i_y))>1.0e19) {

					   newValue = FEDivergenceCorrection(ii, startPtX, startPtY, (sum-i_y), -i_y, m,n, 4); // (startPtX, startPtY, incr_i_x, incr_i_y, maxincr_i_x, maxincr_i_y, first)
					   theGridPlanes[ii]->setGridValue(startPtX+(sum-i_y), startPtY-i_y, newValue);
					   opserr<<"correcting: ["<<  (sum-i_y)<<","<<-i_y<<"] to"<<newValue<<endln;
					   output_v <<"correction::::: jj="<<startPtX+(sum-i_y)<<", kk="<<startPtY-i_y<<", value = "<< newValue<<endln;
				   }


				   //opserr<<"correcting: ["<<  (sum-i_y)<<","<<-i_y<<"]"<<endln;
				}
			}//if
			else{
				startPt = startPt+1;
				int min = sum-1;
				if (sum-1>n) min = n;
				
				for(int i_y = startPt; i_y<= min; i_y++){ 

				   opserr<<"checking: ["<<  (sum-i_y)<<","<<-i_y<<"]"<<endln;
				   if (fabs(theGridPlanes[ii]->getSavedValueOnGrid(startPtX+(sum-i_y), startPtY-i_y))>1.0e19) {

					   newValue = FEDivergenceCorrection(ii, startPtX, startPtY, (sum-i_y), -i_y, m,n, 4); // (startPtX, startPtY, incr_i_x, incr_i_y, maxincr_i_x, maxincr_i_y, first)
					   theGridPlanes[ii]->setGridValue(startPtX+(sum-i_y), startPtY-i_y, newValue);
					   opserr<<"correcting: ["<<  (sum-i_y)<<","<<-i_y<<"] to"<<newValue<<endln;
					   output_v <<"correction::::: jj="<<startPtX+(sum-i_y)<<", kk="<<startPtY-i_y<<", value = "<< newValue<<endln;
				   }

				   // opserr<<"correcting: ["<<  (sum-i_y)<<","<<-i_y<<"]"<<endln;
				}
			}//else
			
		}// for sum





		// --- 2nd (-1,1)  m,n>0---
		m = startPtX;
		n = numOfGridYPt - startPtY-1;

        opserr<<"======================================="<<endln;
		 startPt=1;
		for ( sum =2; sum <=m+n; sum++){
			if(sum-1<=m){
				int min = sum-1;
				if (sum-1>n) min = n;

				for(int i_y = startPt; i_y<= min; i_y++){ 
				   opserr<<"checking: ["<<  -(sum-i_y)<<","<<i_y<<"]"<<endln;
				   if (fabs(theGridPlanes[ii]->getSavedValueOnGrid(startPtX-(sum-i_y), startPtY+i_y))>1.0e19) {

					   newValue = FEDivergenceCorrection(ii, startPtX, startPtY, -(sum-i_y), i_y, m,n, 2); // (startPtX, startPtY, incr_i_x, incr_i_y, maxincr_i_x, maxincr_i_y, first)
					   theGridPlanes[ii]->setGridValue(startPtX-(sum-i_y), startPtY+i_y, newValue);
					   opserr<<"correcting: ["<<  -(sum-i_y)<<","<<i_y<<"] to"<<newValue<<endln;
					   output_v <<"correction::::: jj="<<startPtX-(sum-i_y)<<", kk="<<startPtY+i_y<<", value = "<< newValue<<endln;
				   }
				   //opserr<<"correcting: ["<<  -(sum-i_y)<<","<<i_y<<"]"<<endln;
				}
			}//if
			else{
				startPt = startPt+1;
				int min = sum-1;
				if (sum-1>n) min = n;
				
				for(int i_y = startPt; i_y<= min; i_y++){ 
				   opserr<<"checking: ["<<  -(sum-i_y)<<","<<i_y<<"]"<<endln;
				   if (fabs(theGridPlanes[ii]->getSavedValueOnGrid(startPtX-(sum-i_y), startPtY+i_y))>1.0e19) {

					   newValue = FEDivergenceCorrection(ii, startPtX, startPtY, -(sum-i_y), i_y, m,n, 2); // (startPtX, startPtY, incr_i_x, incr_i_y, maxincr_i_x, maxincr_i_y, first)
					   theGridPlanes[ii]->setGridValue(startPtX-(sum-i_y), startPtY+i_y, newValue);
					   opserr<<"correcting: ["<<  -(sum-i_y)<<","<<i_y<<"] to"<<newValue<<endln;
					   output_v <<"correction::::: jj="<<startPtX-(sum-i_y)<<", kk="<<startPtY+i_y<<", value = "<< newValue<<endln;
				   }


				   //opserr<<"correcting: ["<<  -(sum-i_y)<<","<<i_y<<"]"<<endln;
				}
			}//else
			
		}// for sum

		// --- 1st(1,1)  m,n>0 ---
		m = numOfGridXPt - startPtX-1;
		n = numOfGridYPt - startPtY-1;

        opserr<<"======================================="<<endln;
		 startPt=1;
		for ( sum =2; sum <=m+n; sum++){
			if(sum-1<=m){
				int min = sum-1;
				if (sum-1>n) min = n;

				for(int i_y = startPt; i_y<= min; i_y++){ 
					   
				   opserr<<"checking: ["<<  (sum-i_y)<<","<<i_y<<"]"<<endln;
				   double debug = theGridPlanes[ii]->getSavedValueOnGrid(startPtX+sum-i_y, startPtY+i_y); // debug purpose
				   if (fabs(theGridPlanes[ii]->getSavedValueOnGrid(startPtX+sum-i_y, startPtY+i_y))>1.0e19) {

					   newValue = FEDivergenceCorrection(ii, startPtX, startPtY, sum-i_y, i_y, m,n, 1); // (startPtX, startPtY, incr_i_x, incr_i_y, maxincr_i_x, maxincr_i_y, first)
					   theGridPlanes[ii]->setGridValue(startPtX+sum-i_y, startPtY+i_y, newValue);
					   opserr<<"correcting: ["<<  (sum-i_y)<<","<<i_y<<"] to"<<newValue<<endln;
					   output_v <<"correction::::: jj="<<startPtX+(sum-i_y)<<", kk="<<startPtY+i_y<<", value = "<< newValue<<endln;
				   }
				}
			}//if
			else{
				startPt = startPt+1;
				int min = sum-1;
				if (sum-1>n) min = n;
				
				for(int i_y = startPt; i_y<= min; i_y++){ 
				   opserr<<"checking: ["<<  (sum-i_y)<<","<<i_y<<"]"<<endln;
				   if (fabs(theGridPlanes[ii]->getSavedValueOnGrid(startPtX+sum-i_y, startPtY+i_y))>1.0e19) {

					   newValue = FEDivergenceCorrection(ii, startPtX, startPtY, sum-i_y, i_y, m,n, 1); // (startPtX, startPtY, incr_i_x, incr_i_y, maxincr_i_x, maxincr_i_y, first)
					   theGridPlanes[ii]->setGridValue(startPtX+sum-i_y, startPtY+i_y, newValue);
					   opserr<<"correcting: ["<<  (sum-i_y)<<","<<i_y<<"] to"<<newValue<<endln;
					   output_v <<"correction::::: jj="<<startPtX+(sum-i_y)<<", kk="<<startPtY+i_y<<", value = "<< newValue<<endln;
				   }
				}
			}//else
			
		}// for sum


		// =======================================

			output_vis<<"==================================="<<endln;
			output_vis<<"===========After correction ==========="<<endln;

			// --- write for visulaze --
			for( jj=0; jj< numOfGridXPt; jj++){
				for(int kk=0; kk< numOfGridYPt; kk++){
					double x = rule_1->getPointCoordinate(jj);
					double y = rule_2->getPointCoordinate(kk);
					double z = theGridPlanes[ii]->getSavedValueOnGrid(jj,kk);
					output_vis<<x<<"   "<<y<<"   "<<z<<"\n";
					opserr<<"After correction\n";
					opserr<<"Plain:"<<ii<<": "<<x<<","<<y<<","<<z<<"\n";	
				
				}
			} //jj



			output_vis.close();
			output_v.close();
		}  // for ii



	} // else if "BivariateDecomposition"



	if (theMatrixOperations !=0) delete theMatrixOperations;

	//5. fit curve
	theSurfaceDesign->fitCurve();

	//6. simulation
	theRespSurfaceSimulation->runSimulationAnalysis();

	output<<"failure prob:"<< theRespSurfaceSimulation->getFailureProbability()<<endln; 
	output<<"COV:"<< theRespSurfaceSimulation->getCov()<<endln; 

	

	file_coeff.close();	
	
	return 0;

  }

void DP_RSM_Sim::setGridInfo(Vector *gridInfo)
{





/*
	

  pointer:    0    1    2    3    4    5    6    7    8    9  size = numOfPrincipalAxes+1          
  #eigenVect  0    1    2    3    4    5    6    7    8    9 =numOfPrincipalAxes
			 _________________________________________________ 
			|    |    |    |    |    |    |    |    |    |    |
			|    |    |    |    |    |    |    |    |    |    |
			|____|____|____|____|____|____|____|____|____|____|

*/
	// ---1. create principalAxes

	for (int i =0; i< numOfPrincipalAxes+1; i++){  // numOfPrincipalAxes+1 because dp is accounted as a direction
		if (thePrincipalAxes[i] !=0) delete thePrincipalAxes[i]; 
		thePrincipalAxes[i] = new PrincipalAxis(i+1, this->theExpPtRule);
	}

	
	if (gridInfo==0) return ;

   
	// --2, set grid info ----
	

	int numInfo = (int)gridInfo->Size()/4;
	
	int iNum=0;
  
	opserr<<"grid:"<<(*gridInfo)<<endln;

    if (fabs((*gridInfo)(iNum) +1)<1.0e-10) {  // -1: all axis and dp direction
		    iNum++;
			double beginY=(*gridInfo)(iNum++);
			double endY=(*gridInfo)(iNum++);
			int numOfGridY =(*gridInfo)(iNum++);
			for (int i=0; i<numOfPrincipalAxes+1;i++){
				(thePrincipalAxes[i]->getExperimentalPointRule())->setInfo(numOfGridY, beginY,endY);
			}
		
	
	}
	if (iNum>=(int)gridInfo->Size()) return;

 
	for ( ; iNum<(int)gridInfo->Size();){
		    int num=(*gridInfo)(iNum++);
			double beginX=(*gridInfo)(iNum++);
			double endX=(*gridInfo)(iNum++);
			int numOfGridX =(*gridInfo)(iNum++);
	
			if (num >numOfPrincipalAxes){
			     opserr<<"warning:  DP_RSM_Sim::setGridInfo(), number of information more than number of principalPlane."<<endln;
			}
			else {
				(thePrincipalAxes[num]->getExperimentalPointRule())->setInfo(numOfGridX, beginX,endX);			
			}
	
	} //for
	
return ;

}

int DP_RSM_Sim::getNumOfAxis(int Axis, int numGridPlane)
{
	int i =0; 
	int k=numGridPlane; 
	int n = numOfPrincipalAxes;
	
	while ( k>n){
	
		k -=n;
		n--;
		i++;

	}
 

	if (Axis ==1) return i;

	else if (Axis ==2){
		 return k+i;
	}


	else return -1;

}


/*  mapping between num of grid plane and principal axial. (0) means dp direction, 1,2 .. means the ith eigenVector.
   number the grid planes in a way like:  (numOfPrincipalAxes=4)

  dp      1st      2nd      3rd      4th
  0       1        2        3         4 


  return pointer index    numGridPlane      pointer of GridPlane[i]
  i     j                 numGridPlane        numGridPlane-1
  0     1                     1                    0
  0     2                     2                    1
  0     3                     3                    2
  0     4                     4                    3
  1     2                     5                    4
  1     3                     6                    5
  1     4                     7                    6
  2     3                     8                    7 
  2     4                     9                    8
  3     4                     10                   9
  
  


*/

/*

PrincipalAxes[i]: i is pointer index

  pointer id i      axis
     0               dp
     1              1st eigenVector
     2              2nd
     3              3rd
     4              4th




*/

int DP_RSM_Sim::getNumOfGridPlane(int i, int j)
{
	int k =j;
	int n = numOfPrincipalAxes;
	while (i>0){
		i--;
		n--;
		k +=n;

	
	}

	return k;

}


//FEDivergenceCorrection(ii, startPtX, startPtY, sum-i_y, i_y, m,n, 1); // (startPtX, startPtY, incr_i_x, incr_i_y, maxincr_i_x, maxincr_i_y, first)
double DP_RSM_Sim::FEDivergenceCorrection(int numGridPlane, int startPtX, int startPtY, int incr_i_x, int incr_i_y, int m, int n, int type)
{

	int i_x = startPtX + incr_i_x;
	int i_y = startPtY + incr_i_y;


	int ii = numGridPlane;

	PrincipalAxis * axis_1 = theGridPlanes[ii]->getAxisPtr(1);
	PrincipalAxis * axis_2 = theGridPlanes[ii]->getAxisPtr(2);

	ExperimentalPointRule1D * rule_1 = axis_1->getExperimentalPointRule();
	ExperimentalPointRule1D * rule_2 = axis_2->getExperimentalPointRule();



	double x_0, x_1, x, x_2;
	double y_0, y_1, y, y_2;

	x = rule_1->getPointCoordinate(i_x);
	y = rule_2->getPointCoordinate(i_y);

	double v_x_0 = 0;
	double v_x_1 = 0;
	
	double v = 0;
	double v_x_2 = 0;

	double v_y_0 = 0;
	double v_y_1 = 0;

	double v_y_2 = 0;
	

	switch( type ) 
	{
		case 1:

			x_0 = rule_1->getPointCoordinate(i_x-2);
			x_1 = rule_1->getPointCoordinate(i_x-1);
	//		x = rule_1->getPointCoordinate(i_x);

			y_0 = rule_2->getPointCoordinate(i_y-2);
			y_1 = rule_2->getPointCoordinate(i_y-1);
	//		y = rule_2->getPointCoordinate(i_y);

			v_x_0 = theGridPlanes[ii]->getSavedValueOnGrid(i_x-2, i_y);
			v_x_1 = theGridPlanes[ii]->getSavedValueOnGrid(i_x-1, i_y);
	
			v_y_0 = theGridPlanes[ii]->getSavedValueOnGrid(i_x, i_y-2);
			v_y_1 = theGridPlanes[ii]->getSavedValueOnGrid(i_x, i_y-1);

			if (incr_i_x<m){
				v_x_2 = theGridPlanes[ii]->getSavedValueOnGrid(i_x+1, i_y);
				if (fabs(v_x_2)<1.0e19) { // this is a 'hole' and need to be 'filled'
					x_2 = rule_1->getPointCoordinate(i_x+1);
					v = (v_x_2-v_x_1)/(x_2-x_1)*(x-x_1)+v_x_1;
					return v;
				} 
			}
			if (incr_i_y<n){
				v_y_2 = theGridPlanes[ii]->getSavedValueOnGrid(i_x, i_y+1);
				if (fabs(v_y_2)<1.0e19) { // this is a 'hole' and need to be 'filled'
					y_2 = rule_2->getPointCoordinate(i_y+1);
					v = (v_y_2-v_y_1)/(y_2-y_1)*(y-y_1)+v_y_1;
					return v;
				} 
			}

			// --- exterplolation ---

			if ((fabs(v_x_0)<1.0e19) &&(fabs(v_y_0)<1.0e19))
				v = ((v_x_1-v_x_0)/(x_1-x_0)*(x-x_0)+v_x_0+(v_y_1-v_y_0)/(y_1-y_0)*(y-y_0)+v_y_0)/2.0;
			else if (fabs(v_x_0)<1.0e19)
				v = (v_x_1-v_x_0)/(x_1-x_0)*(x-x_0)+v_x_0;
			else if (fabs(v_y_0)<1.0e19)
				v= (v_y_1-v_y_0)/(y_1-y_0)*(y-y_0)+v_y_0;
			else {
				opserr<<"Fatal: the range of points in axis "<< ii<<" is too big."<<endln;
				exit(-1);
			}
			
			
			break;
		case 2 :   // (-1,1)

			x_0 = rule_1->getPointCoordinate(i_x+2);
			x_1 = rule_1->getPointCoordinate(i_x+1);

			y_0 = rule_2->getPointCoordinate(i_y-2);
			y_1 = rule_2->getPointCoordinate(i_y-1);


			v_x_0 = theGridPlanes[ii]->getSavedValueOnGrid(i_x+2, i_y);
			v_x_1 = theGridPlanes[ii]->getSavedValueOnGrid(i_x+1, i_y);
	
			v_y_0 = theGridPlanes[ii]->getSavedValueOnGrid(i_x, i_y-2);
			v_y_1 = theGridPlanes[ii]->getSavedValueOnGrid(i_x, i_y-1);

			if (abs(incr_i_x)<m){
				v_x_2 = theGridPlanes[ii]->getSavedValueOnGrid(i_x-1, i_y);
				if (fabs(v_x_2)<1.0e19) { // this is a 'hole' and need to be 'filled'
					x_2 = rule_1->getPointCoordinate(i_x-1);
					v = (v_x_2-v_x_1)/(x_2-x_1)*(x-x_1)+v_x_1;
					return v;
				} 
			}
			if (abs(incr_i_y)<n){
				v_y_2 = theGridPlanes[ii]->getSavedValueOnGrid(i_x, i_y+1);
				if (fabs(v_y_2)<1.0e19) { // this is a 'hole' and need to be 'filled'
					y_2 = rule_2->getPointCoordinate(i_y+1);
					v = (v_y_2-v_y_1)/(y_2-y_1)*(y-y_1)+v_y_1;
					return v;
				} 
			}

			// --- exterplolation ---
			if ((fabs(v_x_0)<1.0e19) &&(fabs(v_y_0)<1.0e19))
				v = ((v_x_1-v_x_0)/(x_1-x_0)*(x-x_0)+v_x_0+(v_y_1-v_y_0)/(y_1-y_0)*(y-y_0)+v_y_0)/2.0;
			else if (fabs(v_x_0)<1.0e19)
				v = (v_x_1-v_x_0)/(x_1-x_0)*(x-x_0)+v_x_0;
			else if (fabs(v_y_0)<1.0e19)
				v= (v_y_1-v_y_0)/(y_1-y_0)*(y-y_0)+v_y_0;
			else {
				opserr<<"Fatal: the range of points in axis "<< ii<<" is too big."<<endln;
				exit(-1);
			}
					
			
			break;

		case 3 : // (-1,-1)
			
			x_0 = rule_1->getPointCoordinate(i_x+2);
			x_1 = rule_1->getPointCoordinate(i_x+1);

			y_0 = rule_2->getPointCoordinate(i_y+2);
			y_1 = rule_2->getPointCoordinate(i_y+1);


			v_x_0 = theGridPlanes[ii]->getSavedValueOnGrid(i_x+2, i_y);
			v_x_1 = theGridPlanes[ii]->getSavedValueOnGrid(i_x+1, i_y);
	
			v_y_0 = theGridPlanes[ii]->getSavedValueOnGrid(i_x, i_y+2);
			v_y_1 = theGridPlanes[ii]->getSavedValueOnGrid(i_x, i_y+1);

			if (abs(incr_i_x)<m){
				v_x_2 = theGridPlanes[ii]->getSavedValueOnGrid(i_x-1, i_y);
				if (fabs(v_x_2)<1.0e19) { // this is a 'hole' and need to be 'filled'
					x_2 = rule_1->getPointCoordinate(i_x-1);
					v = (v_x_2-v_x_1)/(x_2-x_1)*(x-x_1)+v_x_1;
					return v;
				} 
			}
			if (abs(incr_i_y)<n){
				v_y_2 = theGridPlanes[ii]->getSavedValueOnGrid(i_x, i_y-1);
				if (fabs(v_y_2)<1.0e19) { // this is a 'hole' and need to be 'filled'
					y_2 = rule_2->getPointCoordinate(i_y-1);
					v = (v_y_2-v_y_1)/(y_2-y_1)*(y-y_1)+v_y_1;
					return v;
				} 
			}

			// --- exterplolation ---

			if ((fabs(v_x_0)<1.0e19) &&(fabs(v_y_0)<1.0e19))
				v = ((v_x_1-v_x_0)/(x_1-x_0)*(x-x_0)+v_x_0+(v_y_1-v_y_0)/(y_1-y_0)*(y-y_0)+v_y_0)/2.0;
			else if (fabs(v_x_0)<1.0e19)
				v = (v_x_1-v_x_0)/(x_1-x_0)*(x-x_0)+v_x_0;
			else if (fabs(v_y_0)<1.0e19)
				v= (v_y_1-v_y_0)/(y_1-y_0)*(y-y_0)+v_y_0;
			else {
				opserr<<"Fatal: the range of points in axis "<< ii<<" is too big."<<endln;
				exit(-1);
			}
			
			break;
		
		case 4:   //(1,-1)

			x_0 = rule_1->getPointCoordinate(i_x-2);
			x_1 = rule_1->getPointCoordinate(i_x-1);
	//		x = rule_1->getPointCoordinate(i_x);

			y_0 = rule_2->getPointCoordinate(i_y+2);
			y_1 = rule_2->getPointCoordinate(i_y+1);
	//		y = rule_2->getPointCoordinate(i_y);

			v_x_0 = theGridPlanes[ii]->getSavedValueOnGrid(i_x-2, i_y);
			v_x_1 = theGridPlanes[ii]->getSavedValueOnGrid(i_x-1, i_y);
	
			v_y_0 = theGridPlanes[ii]->getSavedValueOnGrid(i_x, i_y+2);
			v_y_1 = theGridPlanes[ii]->getSavedValueOnGrid(i_x, i_y+1);

			if (abs(incr_i_x)<m){
				v_x_2 = theGridPlanes[ii]->getSavedValueOnGrid(i_x+1, i_y);
				if (fabs(v_x_2)<1.0e19) { // this is a 'hole' and need to be 'filled'
					x_2 = rule_1->getPointCoordinate(i_x+1);
					v = (v_x_2-v_x_1)/(x_2-x_1)*(x-x_1)+v_x_1;
					return v;
				} 
			}
			if (abs(incr_i_y)<n){
				v_y_2 = theGridPlanes[ii]->getSavedValueOnGrid(i_x, i_y-1);
				if (fabs(v_y_2)<1.0e19) { // this is a 'hole' and need to be 'filled'
					y_2 = rule_2->getPointCoordinate(i_y-1);
					v = (v_y_2-v_y_1)/(y_2-y_1)*(y-y_1)+v_y_1;
					return v;
				} 
			}

			// --- exterplolation ---

			if ((fabs(v_x_0)<1.0e19) &&(fabs(v_y_0)<1.0e19))
				v = ((v_x_1-v_x_0)/(x_1-x_0)*(x-x_0)+v_x_0+(v_y_1-v_y_0)/(y_1-y_0)*(y-y_0)+v_y_0)/2.0;
			else if (fabs(v_x_0)<1.0e19)
				v = (v_x_1-v_x_0)/(x_1-x_0)*(x-x_0)+v_x_0;
			else if (fabs(v_y_0)<1.0e19)
				v= (v_y_1-v_y_0)/(y_1-y_0)*(y-y_0)+v_y_0;
			else {
				opserr<<"Fatal: the range of points in axis "<< ii<<" is too big."<<endln;
				exit(-1);
			}
						
		
			break;
			
		}// switch
		



	return v;
}
