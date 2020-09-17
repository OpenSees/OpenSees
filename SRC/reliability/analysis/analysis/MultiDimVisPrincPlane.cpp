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
 

#include <MultiDimVisPrincPlane.h>
#include <Matrix.h>
#include <Vector.h>
#include <MatrixOperations.h>
#include <ostream>
#include <math.h>

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////


#ifndef _WIN32
void itoa(int n, char *s, int base);
#endif

MultiDimVisPrincPlane::MultiDimVisPrincPlane(
                    ReliabilityDomain *passedReliabilityDomain,
					FunctionEvaluator *passedGFunEvaluator,
					ProbabilityTransformation *passedProbabilityTransformation,
					char *passedOutputFileName,
					GradientEvaluator * passedGradGEvaluator, Vector * pDesignPt, int numPPlane, int pType, Vector *pVector,Tcl_Interp *passedTclInterp, 
					Matrix * passedHessian, char * passedHessianFile, int pAnalysisType, double pLittleDt
					): ReliabilityAnalysis()
{

	analysisType =pAnalysisType; // time invariant   1: timeVariant
	valuesOfAxis =0;
	valuesG2OfAxis = 0;
	
	HessianFileName = 0;
	if (passedHessianFile !=0) {
		HessianFileName = new char[30];
		strcpy(HessianFileName, passedHessianFile);
	
	}


	this->theInterp=passedTclInterp;
	theReliabilityDomain = passedReliabilityDomain;
	theGFunEvaluator = passedGFunEvaluator;
	theProbabilityTransformation = passedProbabilityTransformation;
    strcpy(outputFileName, passedOutputFileName);
    theGradGEvaluator = passedGradGEvaluator;
	if (pDesignPt ==0) {opserr<<"no designpt."<<endln; exit(-1);};

	theDesignPtXSpace = new Vector (*pDesignPt);
	
	//FMK
	/*
	theProbabilityTransformation->set_x(*theDesignPtXSpace);
	theProbabilityTransformation->transform_x_to_u();
	theDesignPoint->addVector(0.0,theProbabilityTransformation->get_u(),1.0);
	*/

	int numRV = theReliabilityDomain->getNumberOfRandomVariables();
	theDesignPoint = new Vector(numRV);
	theProbabilityTransformation->transform_x_to_u(*theDesignPoint);	

	thePrincipalPlanes = new PrincipalPlane * [20]; // uplimit 20 for now.refer arrayOfTaggedObject.. theComponents = new TaggedObject *[size]
	
	for(int i = 0;i< 20; i++) thePrincipalPlanes[i]=0;
	
    numOfPrinPlane = numPPlane;
	
	if (numOfPrinPlane>= numRV){
		opserr<<"warning: numOfPrinPlane can not be larger than"<<numRV-1<<", and it is reset to "<<numRV-1<<endln;
		numOfPrinPlane = numRV-1;
	}

	if (numPPlane>20) {
		numPPlane=20;
		opserr<<"set numPrincPlane to 20."<<endln;
	}

	type = pType;

	this->setGridInfo(pVector, numOfPrinPlane);
	HessianMatrix =0;
	rotation =0;
	//theHessian=0;
	littleDt = pLittleDt;

	if (passedHessian !=0){
		//HessianMatrix = new Matrix( *passedHessian);
		//theHessian = new Hessian(numRV,theReliabilityDomain,theProbabilityTransformation,
        //                         theGFunEvaluator,theGradGEvaluator, 1.0e-5);
	
		//theHessian->formReducedHessian(theDesignPtXSpace, passedHessian); 

	}


}

MultiDimVisPrincPlane::~MultiDimVisPrincPlane()
{
	for(int i = 0;i< 20; i++)
		if (thePrincipalPlanes[i] !=0){
			delete thePrincipalPlanes[i];
		}
	delete thePrincipalPlanes;
    if (HessianMatrix !=0) delete HessianMatrix;
	if (rotation !=0) delete rotation;
	//if (theHessian !=0) delete theHessian;
	if (HessianFileName !=0) delete HessianFileName;
	



}
int MultiDimVisPrincPlane::analyze()
{

// --- type =0:   Limit state surface -----
// --- type =1:   Limit function -------
	if (type ==0){
		opserr<<"not implemented yet."<<endln;
		exit(-1);
	
	}
    //	double tolFun =1.0e-6;

// --1. create principalplanes, set grid info ----

	
	if (thePrincipalPlanes[0]==0) {
		for(int i =0; i<numOfPrinPlane; i++){

			thePrincipalPlanes[i]=new PrincipalPlane(i+1, theDesignPoint, 
								   0, 0, theProbabilityTransformation,
								   theGFunEvaluator, 0.0);
		//	thePrincipalPlanes[i]->setDesignPoint( theDesignPoint);

		}  /// 
	}


	ofstream debug("debug.out"  );

     
	//2.  --- compute and form Hessian by FDM; compute A, Rotation R, A=R*H*R'/|dG|, but one dim less.

	int numRV = theReliabilityDomain->getNumberOfRandomVariables();
  	
	//if (theHessian ==0){
		//theHessian = new Hessian(numRV,	theReliabilityDomain,theProbabilityTransformation,
        //                         theGFunEvaluator,theGradGEvaluator, 1.0e-5);
		//theHessian->formReducedHessian(theDesignPtXSpace);
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


	
	for(int i =0; i<numRV; i++)
		for(int j=0; j<numRV; j++)
			debug<<(*rotation)(i,j)<<endln;
	

	

    MatrixOperations * theMatrixOperations = new MatrixOperations(*HessianMatrix);
	
	// 3. -- eigen analysis, insert into thePrincipalPlane array according to magnitute

	if (numOfPrinPlane*2>=numRV-1){
		theMatrixOperations->performEigenAnalysis(1,numRV-1);

		// save into principalplanes
		double eigenValue;
		Vector eigenVector(numRV);
		Vector tmp(numRV-1);
		Vector tmp1(numRV);
		tmp1(numRV-1)=0.0;
		
		int i;
		for(i=0; i<numOfPrinPlane; i++){
			eigenValue = theMatrixOperations->getEigenvalue(i+1);
			tmp.addVector(0.0, theMatrixOperations->getEigenvector(i+1),1.0);
			for (int j =0; j<numRV-1; j++)
				tmp1(j)=tmp(j);

			eigenVector = (*rotation) ^ tmp1;
	//		opserr<<"====numOfPPlane:" << numOfPrinPlane << "       \n eigenVector:" <<eigenVector<<" =====\n";

			thePrincipalPlanes[i]->setCurvature(eigenValue);
			thePrincipalPlanes[i]->setEigenVector(&eigenVector);

		}
			
		int insertPlace = 0;			
		for( i=numRV-2; i>0; i--){

			eigenValue = theMatrixOperations->getEigenvalue(i+1);
			// -- comp abs of eigenValue


			while(insertPlace < numOfPrinPlane){
				if ( fabs(thePrincipalPlanes[insertPlace]->getCurvature())>=fabs(eigenValue)){
					insertPlace++;
				}
				else 
					break;

			}

			if (insertPlace != numOfPrinPlane){
				for( int k=numOfPrinPlane-1; k>insertPlace; k--){
					thePrincipalPlanes[k]->copyValues( thePrincipalPlanes[k-1]);
				}

				tmp.addVector(0.0, theMatrixOperations->getEigenvector(i+1),1.0);
				for (int j =0; j<numRV-1; j++)
					tmp1(j)=tmp(j);
				tmp1(numRV-1)=0.0;

				eigenVector = (*rotation) ^ tmp1;
			//    opserr<<"====inserted numOfPPlane:"<<numOfPrinPlane<<"       \n eigenVector:"<< eigenVector<<" =====\n";
				//thePrincipalPlanes[insertPlace]->setNumOfPlane(insertPlace+1);
				thePrincipalPlanes[insertPlace]->setCurvature(eigenValue);
				thePrincipalPlanes[insertPlace]->setEigenVector(&eigenVector);
				thePrincipalPlanes[insertPlace]->cleanGridMatrix();
				insertPlace++;
			}
			else if (insertPlace == numOfPrinPlane){ i=0; } // to exit for loop
		}  //for

	}
	else{
		theMatrixOperations->performEigenAnalysis(1,numOfPrinPlane);
		
		// save into principalplanes
		double eigenValue;
		Vector eigenVector(numRV);
		Vector tmp(numRV-1);
		Vector tmp1(numRV);
		tmp1(numRV-1)=0.0;
		
		int i;
		for(i=0; i<numOfPrinPlane; i++){
			eigenValue = theMatrixOperations->getEigenvalue(i+1);
			tmp.addVector(0.0, theMatrixOperations->getEigenvector(i+1),1.0);
			for (int j =0; j<numRV-1; j++)
				tmp1(j)=tmp(j);

			eigenVector = (*rotation) ^ tmp1;

			thePrincipalPlanes[i]->setCurvature(eigenValue);
			thePrincipalPlanes[i]->setEigenVector(&eigenVector);

		}

		theMatrixOperations->performEigenAnalysis(numRV-numOfPrinPlane,numRV-1);

		int insertPlace = 0;
		for( i=numRV-2; i>=numRV-numOfPrinPlane-1; i--){

			eigenValue = theMatrixOperations->getEigenvalue(i+1);
			// -- comp abs of eigenValue
			
			
			while(insertPlace < numOfPrinPlane){
				if ( fabs(thePrincipalPlanes[insertPlace]->getCurvature())>=fabs(eigenValue)){
					insertPlace++;
				}
				else 
					break;

			} //while
			

			if (insertPlace != numOfPrinPlane){
				for( int k=numOfPrinPlane-1; k>insertPlace; k--){
					thePrincipalPlanes[k]->copyValues( thePrincipalPlanes[k-1]);
				}

				tmp.addVector(0.0, theMatrixOperations->getEigenvector(i+1),1.0);
				for (int j =0; j<numRV-1; j++)
					tmp1(j)=tmp(j);
				tmp1(numRV-1)=0.0;

				eigenVector = (*rotation) ^ tmp1;

				//thePrincipalPlanes[insertPlace]->setNumOfPlane(insertPlace+1);
				thePrincipalPlanes[insertPlace]->setCurvature(eigenValue);
				thePrincipalPlanes[insertPlace]->setEigenVector(&eigenVector);
				thePrincipalPlanes[insertPlace]->cleanGridMatrix();
				insertPlace++;
				// note gridInfo does not change..
		
			}
			else {i=numRV-numOfPrinPlane-1; }  // stop for loop
	

		}  //for
	} //else if (numOfPrinPlane*2>=numRV-1)


	ofstream file_coeff("axis_recorder.out", ios::out);



	for (int k =0; k<numOfPrinPlane; k++){

		double a = thePrincipalPlanes[k]->getCurvature();
		opserr<<"the "<<k+1<<"th curvature:"<<a<<endln;
		file_coeff<<"the "<<k+1<<"th curvature:"<<a<<endln;

		
		Vector * aa = thePrincipalPlanes[k]->getEigenVectorPtr();
		opserr<<"the "<<k+1<<"th direction:\n";
		file_coeff<<"the "<<k+1<<"th direction:\n";

		for( int m =0; m< aa->Size(); m++){
			opserr<<(*aa)(m)<<endln;
			file_coeff<<(*aa)(m)<<endln;
	
		}

	} //for k

	file_coeff.close();

	// 4. --- compute the grid values save to file



/*	opserr<<"=============checking purpose:==============\n";
	opserr.setPrecision(16);

	opserr<<"1.)Hessian:"<< theHessian->getHessianApproximation()<<endln;
    

	opserr<<"2.) principalPlane          eigenValue             eigenVector "<<endln;	
	for(int  ii=0; ii<numOfPrinPlane; ii++) {
		opserr<< ii+1<<"        "<< thePrincipalPlanes[ii]->getCurvature()<<"       "<<*(thePrincipalPlanes[ii]->getEigenVectorPtr())<<endln;
	}

*/
	
	double valueOfG =0;
	double valueOfG2 =0;


	char fileName[30];
	strcpy(fileName, outputFileName);
	char seps[]   = ".";
    char * token;

    

	char cmd[50]="remove sensitivityAlgorithm";
	Tcl_Eval(theInterp, cmd);

	if (analysisType ==0){ // time invariant
		for( int   ii=0; ii<numOfPrinPlane; ii++) {

			// == using zero finding to find the zero point ==
			token = strtok( outputFileName, seps );
			strcpy(fileName, token);
    
			token = strtok( NULL, seps );
			

			char extension[5]=".";
			if (token !=NULL)    
				strcat(extension, token);
			else
				strcat(extension, "out");
   			char theStr[3];

				// 4.1 -- open file to write

				//itoa(ii+1,theStr,10);
				sprintf(theStr,"%d",ii+1);
				strcat(fileName, theStr);
				strcat(fileName, extension);
    			
				ofstream output( fileName, ios::out );


				// 4.2 -- compute x,y write into file


			


				double x, y;

				int numOfGridXPt = thePrincipalPlanes[ii]->getNumOfGridXPt();
				double beginX = thePrincipalPlanes[ii]->getBeginOfGridX();
				double endX= thePrincipalPlanes[ii]->getEndOfGridX();
				double stepX = (endX - beginX) / (numOfGridXPt-1);
			

			

				int numOfGridYPt = thePrincipalPlanes[ii]->getNumOfGridYPt();
				double beginY = thePrincipalPlanes[ii]->getBeginOfGridY();
				double endY= thePrincipalPlanes[ii]->getEndOfGridY();
				double stepY = (endY - beginY) / (numOfGridYPt-1);


				// -- grid on axils --

				if (valuesOfAxis ==0)
				{
					valuesOfAxis = new Vector(numOfGridYPt);
					for(int kk=0; kk< numOfGridYPt; kk++){
						x=0;
						y=beginY + kk*stepY;
						if (fabs(y)<1.0e-14) {
							valueOfG =0;
						}
						else{
							valueOfG = thePrincipalPlanes[ii]->getValueOnGrid(x,y);
						}
						(*valuesOfAxis)(kk) = valueOfG;

					}
				}

				for(int jj=0; jj< numOfGridXPt; jj++){
					for(int kk=0; kk< numOfGridYPt; kk++){
					//	output<<"ii="<<ii<<",jj="<<jj<<"kk="<<kk<<"\n";
						x= beginX + jj*stepX;
						y= beginY + kk*stepY;
						if (fabs(x)<1.0e-14) {
							valueOfG =(*valuesOfAxis)(kk);
							thePrincipalPlanes[ii]->setGridValue(jj,kk,valueOfG);
						}
						else{
							valueOfG = thePrincipalPlanes[ii]->getValueOnGrid(x,y);
							thePrincipalPlanes[ii]->setGridValue(jj,kk,valueOfG);
						}
						// write x,y valueOfG

						output<<x<<"   "<<y<<"   "<<valueOfG<<"\n";
						opserr<<x<<","<<y<<","<<valueOfG<<"\n";			
					}
				}

				// 4.3 --close file --
				output.close();
		} // for each principalPlanes
	} // if analysisType==0
	else{  // time variant
		
		
		//double littleDt =1.0e-3;
		char fileName2[40];
		for( int ii=0; ii<numOfPrinPlane; ii++) {

			// == using zero finding to find the zero point ==
			token = strtok( outputFileName, seps );
			strcpy(fileName, token);
    

			strcpy(fileName2,fileName);
			strcat(fileName2, "_dt");

			token = strtok( NULL, seps );
			
			
			char extension[5]=".";

			if (token !=NULL) 
				strcat(extension, token);
			else
				strcat(extension, "out");
   			char theStr[3];

				// 4.1 -- open file to write

				//itoa(ii+1,theStr,10);
				sprintf(theStr,"%d",ii+1);
				strcat(fileName, theStr);
				strcat(fileName, extension);

				strcat(fileName2, theStr);
				strcat(fileName2, extension);
    			
				ofstream output( fileName, ios::out );
				ofstream output2( fileName2, ios::out );



				// 4.2 -- compute x,y write into file


			


				double x, y;

				int numOfGridXPt = thePrincipalPlanes[ii]->getNumOfGridXPt();
				double beginX = thePrincipalPlanes[ii]->getBeginOfGridX();
				double endX= thePrincipalPlanes[ii]->getEndOfGridX();
				double stepX = (endX - beginX) / (numOfGridXPt-1);
			

			

				int numOfGridYPt = thePrincipalPlanes[ii]->getNumOfGridYPt();
				double beginY = thePrincipalPlanes[ii]->getBeginOfGridY();
				double endY= thePrincipalPlanes[ii]->getEndOfGridY();
				double stepY = (endY - beginY) / (numOfGridYPt-1);


				// -- grid on axils that is Y axis --
				
				
				if (valuesG2OfAxis ==0) { 
					valuesG2OfAxis = new Vector(numOfGridYPt);
				}

				if (valuesOfAxis ==0)
				{
					valuesOfAxis = new Vector(numOfGridYPt);
					

					for(int kk=0; kk< numOfGridYPt; kk++){
						x=0;
						y=beginY + kk*stepY;
						if (fabs(y)<1.0e-14) {
							
							valueOfG =0;
							valueOfG2 = thePrincipalPlanes[ii]->getValueG2OnGrid(x,y, valueOfG, littleDt);  
						}
						else{
							valueOfG = thePrincipalPlanes[ii]->getValueOnGrid(x,y);
							valueOfG2 = thePrincipalPlanes[ii]->getValueG2OnGrid(x,y, valueOfG, littleDt);
						}
						(*valuesOfAxis)(kk) = valueOfG;
						(*valuesG2OfAxis)(kk) = valueOfG2;  

					}
				}

				for(int jj=0; jj< numOfGridXPt; jj++){
					for(int kk=0; kk< numOfGridYPt; kk++){
					//	output<<"ii="<<ii<<",jj="<<jj<<"kk="<<kk<<"\n";
						x= beginX + jj*stepX;
						y= beginY + kk*stepY;
				
						



						if (fabs(x)<1.0e-14) {
							valueOfG =(*valuesOfAxis)(kk);
							valueOfG2 =(*valuesG2OfAxis)(kk);
							thePrincipalPlanes[ii]->setGridValue(jj,kk,valueOfG);
							thePrincipalPlanes[ii]->setGridValueG2(jj,kk,valueOfG2); 
						}
						else{
							valueOfG = thePrincipalPlanes[ii]->getValueOnGrid(x,y);
							valueOfG2 = thePrincipalPlanes[ii]->getValueG2OnGrid(x,y,valueOfG, littleDt);
							thePrincipalPlanes[ii]->setGridValue(jj,kk,valueOfG);
							thePrincipalPlanes[ii]->setGridValueG2(jj,kk,valueOfG2);
						}
						// write x,y valueOfG

						output<<x<<"   "<<y<<"   "<<valueOfG<<"\n";
						output2<<x<<"   "<<y<<"   "<<valueOfG2<<"\n";

						opserr<<x<<","<<y<<","<<valueOfG<<"\n";	
						opserr<<x<<","<<y<<","<<valueOfG2<<"\n";							
					}
				}

				// 4.3 --close file --
				output.close();
				output2.close();

		} // for each principalPlanes	
	
	
	
	}
	if (theMatrixOperations !=0) { delete theMatrixOperations; theMatrixOperations =0;}
	if (valuesOfAxis !=0) { delete valuesOfAxis;  valuesOfAxis =0;}

	debug.close();
	return 0;
}



int MultiDimVisPrincPlane::setGridInfo(Vector *gridInfo, int numOfPrinPlane )
{

	if (gridInfo==0) return 0;

	if (numOfPrinPlane>20) {
		opserr<<"warning: MultiDimVisPrincPlane, numOfPrincipalPlane limit =20, automatic reset to 20"<<endln; 
		numOfPrinPlane=20;
	}
	if (valuesOfAxis !=0) { delete valuesOfAxis;  valuesOfAxis =0;}
    
	// --1. create principalplanes, set grid info ----
	for(int i =0; i<numOfPrinPlane; i++){

		thePrincipalPlanes[i]=new PrincipalPlane(i+1, theDesignPoint, 
							   0, 0, theProbabilityTransformation,
							   theGFunEvaluator, 0.0);
		// thePrincipalPlanes[i]->setDesignPoint( theDesignPoint);

	}  /// 

	int numInfo = (int)gridInfo->Size()/4;

	


	
	int iNum=0;
  
	opserr<<"grid:"<<(*gridInfo)<<endln;
	if (fabs((*gridInfo)(iNum) -0)<1.0e-10) {
		    iNum++;
			double beginY=(*gridInfo)(iNum++);
			double endY=(*gridInfo)(iNum++);
			int numOfGridY =(*gridInfo)(iNum++);
			for (int i=0; i<numOfPrinPlane;i++){
				thePrincipalPlanes[i]->setGridYInfo(numOfGridY, beginY,endY);
			}
		
	
	}

	for ( ; iNum<(int)gridInfo->Size();){
		    int num=(*gridInfo)(iNum++);
			double beginX=(*gridInfo)(iNum++);
			double endX=(*gridInfo)(iNum++);
			int numOfGridX =(*gridInfo)(iNum++);
	
			if (num >numOfPrinPlane){
			     opserr<<"warning:  MultiDimVisPrincPlane::setGridInfo(), number of information more than number of principalPlane."<<endln;
			}
			else {
				thePrincipalPlanes[num-1]->setGridXInfo(numOfGridX, beginX,endX);
			}
			
		
	}
	
return 0;
}

