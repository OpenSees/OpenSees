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

// $Revision: 1.2 $
// $Date: 2008-05-13 16:30:27 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/telm/ThresholdIncInitialPointBuilder.cpp,v $
#include <math.h>
#include <ThresholdIncInitialPointBuilder.h>
ThresholdIncInitialPointBuilder::ThresholdIncInitialPointBuilder(ReliabilityDomain *passedReliabilityDomain,
						FunctionEvaluator* passedGFunEvaluator,
						FindDesignPointAlgorithm* passedFindDesignPointAlgorithm,
						//NewSearchWithStepSizeAndStepDirection* passedFindDesignPointAlgorithm,
						int passedmaxDivide,
						double passedeps,
						bool passedstart_mirror,
						int passedMaxLineSearch,
						double passedRatio, 
						bool passedprint)
:InitialPointBuilder(passedReliabilityDomain, passedGFunEvaluator)
//:MirrorImageInitialPointBuilder
//(passedReliabilityDomain,passedGFunEvaluator,passedMaxLineSearch,
// passedRatio,passedprint)
{
    theFindDesignPointAlgorithm=passedFindDesignPointAlgorithm;
	eps=passedeps;
	maxDivide=passedmaxDivide;
	start_mirror=passedstart_mirror;
	prevDesign=0;
	xtemp=0;
	xcompleted=0;
	if(print){
		output.open("ThresholdIncInitialPointBuilder.txt", ios::out);
		output << "\n";
		output << "ThresholdIncInitialPointBuilder::ThresholdIncInitialPointBuilder\n";
		output << "\n";
		output << "eps  "<<eps<<"\n";
		output << "maxDivide "<<maxDivide<<"\n";
		output.flush();
	}
} 
ThresholdIncInitialPointBuilder::~ThresholdIncInitialPointBuilder()
{
	if(prevDesign!=0){ delete prevDesign; prevDesign=0;}
	if(xtemp!=0){ delete xtemp; xtemp=0;}
	if(xcompleted!=0){ delete xcompleted; xcompleted=0;}
} 
void ThresholdIncInitialPointBuilder::clear()
{
//	this->clearMirrorImageInitialPointBuilder();
	if(prevDesign!=0){ delete prevDesign; prevDesign=0;}
	if(xtemp!=0){ delete xtemp; xtemp=0;}
	if(xcompleted!=0){ delete xcompleted; xcompleted=0;}
} 
Vector ThresholdIncInitialPointBuilder::buildInitialPoint(int nstep)
{

	static NormalRV aStdNormRV(1,0.0,1.0);

	if(xinitial !=0){delete xinitial;xinitial=0;}
	xinitial=new Vector(*xmean);
	if(nstep<0){
		if(print){
			output << "\n";
			output << "ThresholdIncInitialPointBuilder::buildInitialPoint\n";
			output << "nstep  "<<nstep<<"\n";
			output.flush();
		}
//		AnalysisStep=-nstep;
//		(*xinitial)=this->buildMirrorInitialPoint();
//		numAna=numAnaMirror;
//		if(print){
//			output << "\n";
//			output << "AnalysisStep  "<<AnalysisStep<<"\n";
//			output << "numAna  "<<numAna<<"\n";
//			output.flush();
//		}
		return (*xinitial);
	}else{
		numAna=0;
		if(print){
			output << "\n";
			output << "ThresholdIncInitialPointBuilder::buildInitialPoint\n";
			output << "nstep  "<<nstep<<"\n";
			output.flush();
		}
		int AnalysisStep=nstep;
		theGFunEvaluator->setNsteps(AnalysisStep);
		if(xtemp !=0){delete xtemp;xtemp=0;}
		xtemp=new Vector(*xmean);
		if(xcompleted !=0){delete xcompleted;xcompleted=0;}
		xcompleted=new Vector(*prevDesign);

		double tmpcompletedThreshold=completedThreshold;
		double drespremain=threshold-tmpcompletedThreshold;
		double ampdresp=1.0;
		double dresptrial;
		double thresholdtrial;
		double ampratio,ampf,ampresp;
		double lower=1.0-eps;
		double upper=1.0+eps;
		double ddd,xresp;
		int res1,res2;
		int ibreak=0;
		int iloop=0;
		int idivide=0;
		if(print){
			output << "\n";
			output << "before loop\n";
			output << "threshold  "<<threshold<<"\n";
			output << "completedThreshold  "<<completedThreshold<<"\n";
			output << "tmpcompletedThreshold  "<<tmpcompletedThreshold<<"\n";
			output << "drespremain  "<<drespremain<<"\n";
			output << "ampdresp  "<<ampdresp<<"\n";
			output << "lower  "<<lower<<"\n";
			output << "upper  "<<upper<<"\n";
			output.flush();
		}
		while(ibreak==0){
			iloop++;
            dresptrial=drespremain*ampdresp;
            thresholdtrial=tmpcompletedThreshold+dresptrial;
			ampf=thresholdtrial/tmpcompletedThreshold;
			(*xtemp)=(*xcompleted)*ampf;
			theGFunEvaluator->inactivateSensitivty();
			theGFunEvaluator->setThreshold(0.0);
			res1=theGFunEvaluator->runAnalysis();
			numAna++;
			if ( res1<0 ) {
				opserr << "GFuncError OutCrossingAnalysis::MakeInitialPoint" << res1 << "\n";
				exit(-1);
			}
			xresp=-theGFunEvaluator->evaluateExpression();
			ampresp=xresp/tmpcompletedThreshold;
			if(print){
			output << "\n";
			output << "loop "<<iloop << "\n";
			output << "dresptrial  "<<dresptrial<<"\n";
			output << "thresholdtrial  "<<thresholdtrial<<"\n";
			output << "ampf  "<<ampf<<"\n";
			output << "xresp  "<<xresp<<"\n";
			output << "ampresp  "<<ampresp<<"\n";
			output.flush();
			}
			if(idivide<=maxDivide&&ampresp<1.0){
				ampdresp*=0.5;
				idivide++;
				if(print){
				output << "\n";
				output << "idivide<=maxDivide&&ampresp<1.0\n";
				output << "ampdresp  "<<ampdresp<<"\n";
				output << "idivide  "<<idivide<<"\n";
				output.flush();
				}
			}else{
				ampratio=ampresp/ampf;
				if(print){
				output << "\n";
				output << "ampratio  "<<ampratio<<"\n";
				output.flush();
				}
				if((ampratio>=lower&&ampratio<=upper)||idivide>=maxDivide){
					ddd=fabs(drespremain-dresptrial);
					if(print){
					output << "\n";
					output << "idivide  "<<idivide<<"\n";
					output << "drespremain  "<<drespremain<<"\n";
					output << "dresptrial  "<<dresptrial<<"\n";
					output << "ddd  "<<ddd<<"\n";
					output.flush();
					}
					if(ddd/drespremain<=0.001){
						(*xinitial)=(*xtemp);
						ibreak=1;
						if(print){
						output << "\n";
						output << "ddd/drespremain  "<<ddd/drespremain<<"\n";
						output << "return\n";
						output.flush();
						}
					}else{
						if(print){
						output << "\n";
						output << "try to find des point\n";
						output << "threshold  "<<xresp<<"\n";
						output.flush();
						}
//						theOutCrossingResults->clear(lsf,0, xresp, 0.0, 0, numRV);
						theGFunEvaluator->setThreshold(xresp);
						theFindDesignPointAlgorithm->set_x(*xtemp);
						int iresult=theFindDesignPointAlgorithm->findDesignPoint();

						// save results //
//						theOutCrossingResults->setnumAna(0,theFindDesignPointAlgorithm->getNumberOfEvaluations());
//						theOutCrossingResults->setnumLinSearch(0,numAna);
//						theOutCrossingResults->setnumAnaIncSens(0,theFindDesignPointAlgorithm->getNumberOfSensAna());
//						theOutCrossingResults->setxDesPoints(0,theFindDesignPointAlgorithm->get_x());
//						theOutCrossingResults->setuDesPoints(0,theFindDesignPointAlgorithm->get_u());
//						theOutCrossingResults->setDesAlpha(0,theFindDesignPointAlgorithm->get_alpha());
//						double betaTmp =(theFindDesignPointAlgorithm->get_alpha())^(theFindDesignPointAlgorithm->get_u());
//						double pfTmp = 1.0 - aStdNormRV->getCDFvalue(betaTmp);
//						theOutCrossingResults->setbeta(0,betaTmp);
//						theOutCrossingResults->setpf(0,pfTmp);
//						theOutCrossingResults->setnu(0,0.0);
//						theOutCrossingResults->settime(0,dt*(float)AnalysisStep);
//						theOutCrossingResults->setnumSteps(0,AnalysisStep);
//						theOutCrossingResults->setcheck1(0,theFindDesignPointAlgorithm->get_check1_conv());
//						theOutCrossingResults->setcheck2(0,theFindDesignPointAlgorithm->get_check2_conv());
//						theOutCrossingResults->setcheck1_init(0,theFindDesignPointAlgorithm->get_check1_init());
//						theOutCrossingResults->setcheck2_init(0,theFindDesignPointAlgorithm->get_check2_init());
//						theOutCrossingResults->setiresult(0,iresult);
//						theOutCrossingResults->printSinglePoint(outputFile,lsf,0);
//						theOutCrossingResults->outtoFile();
						numAna=0;
	 					if(iresult<0){ 
							opserr << "Warning!! FORMAnalysis::analyze() - failed while finding the" << endln;
						}
						(*xcompleted)=theFindDesignPointAlgorithm->get_x();
						tmpcompletedThreshold=xresp;
						ampdresp=1.0;
						idivide=0;
						drespremain=threshold-tmpcompletedThreshold;
						if(print){
						output << "\n";
						output << "try to find des point\n";
						output << "threshold "<<xresp<<"\n";
						output << "iresult "<<iresult<<"\n";
						output << "tmpcompletedThreshold "<<tmpcompletedThreshold<<"\n";
						output << "ampdresp "<<ampdresp<<"\n";
						output << "drespremain "<<drespremain<<"\n";
						output.flush();
						}
					}
				}else{
					ampdresp*=0.5;
					idivide++;
					if(print){
					output << "\n";
					output << "divide\n";
					output << "ampdresp "<<ampdresp<<"\n";
					output.flush();
					}
				}
			}	
		}
		theGFunEvaluator->setThreshold(threshold);
		return (*xinitial);
	}
}
void ThresholdIncInitialPointBuilder::setPrevResults
(int prevSteps, double prevTime, double prevThreshold, const Vector& designPoint)
{
	if(prevDesign !=0){delete prevDesign;prevDesign=0;}
	prevDesign=new Vector(*xmean);
	(*prevDesign)=designPoint;
	completedSteps=prevSteps;
	completedTime=prevTime;
	completedThreshold=prevThreshold;
}
void ThresholdIncInitialPointBuilder::setOutCrossingResults
(ofstream& passedoutputFile, OutCrossingResults* passedOutCrossingResults)
{
//	outputFile=passedoutputFile;
	theOutCrossingResults=passedOutCrossingResults;
}	
void ThresholdIncInitialPointBuilder::setMirrorImageExcitation(int i,Vector v)
{
opserr << "ThresholdIncInitialPointBuilder::setMirrorImageExcitation" << endln;
opserr << "This function should not be called" << endln;
}
