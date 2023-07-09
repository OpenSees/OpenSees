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
** ****************************************************************** */
                                                                        
// $Revision: 1.19 $
// $Date: 2010-06-10 18:53:31 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/analysis/SystemAnalysis.cpp,v $


//
// Written by Kevin Mackie (kmackie@mail.ucf.edu)
//

#include <SystemAnalysis.h>
#include <ReliabilityDomain.h>
#include <ReliabilityAnalysis.h>
#include <FunctionEvaluator.h>
#include <LimitStateFunction.h>
#include <LimitStateFunctionIter.h>
#include <Cutset.h>
#include <CutsetIter.h>
#include <RandomNumberGenerator.h>
#include <CStdLibRandGenerator.h>
#include <CorrelatedStandardNormal.h>
#include <NormalRV.h>

#include <math.h>
#include <float.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include <fstream>
#include <iomanip>
#include <iostream>
using std::ifstream;
using std::ios;

SystemAnalysis::SystemAnalysis(ReliabilityDomain *passedReliabilityDomain, 
                               FunctionEvaluator *passedEvaluator, 
                               TCL_Char *passedBeta, TCL_Char *passedRho)
	:ReliabilityAnalysis()
{
	theReliabilityDomain = passedReliabilityDomain;
    theFunctionEvaluator = passedEvaluator;
	strcpy(betaFile,passedBeta);
	strcpy(rhoFile,passedRho);
	
	minLowerBound = 1.0;
	maxUpperBound = 0.0;
	
	int result = initialize();
	if (result < 0) {
		opserr << "SystemAnalysis::SystemAnalysis() ERROR - SystemAnalysis failed to initialize" << endln;
		exit(-1);
	}
	
	sets = new Matrix(2,2);
	permutedBetas = new Vector(2);
	permutedRhos = new Matrix(2,2);
	
}

SystemAnalysis::~SystemAnalysis()
{
	if (allBetas != 0)
		delete allBetas;
	if (rhos != 0 )
		delete rhos;
	if (allPf1s != 0)
		delete allPf1s;
	if (Pmn != 0 )
		delete Pmn;
		
	if (sets != 0 )
		delete sets;
	if (permutedBetas != 0)
		delete permutedBetas;
	if (permutedRhos != 0)
		delete permutedRhos;
}


int
SystemAnalysis::initialize()
{
	// Initial declarations
	double beta;
    NormalRV uRV(1, 0.0, 1.0);
	
	if (strlen(betaFile) > 1 && strlen(rhoFile) > 1) {
		// read data for beta and rho from file, do not take from reliability domain
		ifstream inputBeta( betaFile, ios::in );
		if (!inputBeta) {
			opserr << "SystemAnalysis::initialize ERROR - could not open beta file " << betaFile << endln;
			return -1;
		}
		ifstream inputRho( rhoFile, ios::in );
		if (!inputRho) {
			opserr << "SystemAnalysis::initialize ERROR - could not open rho file " << rhoFile << endln;
			return -1;
		}
		
		numLsf = 0;
		while (inputBeta >> beta)
			numLsf++;
		inputBeta.clear();
		inputBeta.seekg(0);
		
		// allocate arrays
		allBetas = new Vector(numLsf);
		allPf1s = new Vector(numLsf);
		rhos = new Matrix(numLsf,numLsf);
		
		// read data from file into arrays
		for (int i=0; i < numLsf; i++ ) {
			inputBeta >> (*allBetas)(i);
			(*allPf1s)(i) = 1.0 - uRV.getCDFvalue( (*allBetas)(i) );
			
			for (int j=0; j < numLsf; j++ ) {
				if (!inputRho.eof() )
					inputRho >> (*rhos)(i,j);
				else {
					opserr << "SystemAnalysis::initialize ERROR - rho file does not contain the same data size as the beta file"
						   << endln;
					return -1;
				}
			}
		}
		
		inputBeta.close();
		inputRho.close();
		
	} else {
		// proceed with reliability domain data, requires previous FORM analysis objects

		// Number of limit-state functions
		numLsf = theReliabilityDomain->getNumberOfLimitStateFunctions();

		// Number of random variables
		int numRV = theReliabilityDomain->getNumberOfRandomVariables();

		// Allocate vectors to store ALL the betas and alphas
		allBetas = new Vector(numLsf);
		allPf1s = new Vector(numLsf);
		Matrix allAlphas(numRV,numLsf);

		// Loop over number of limit-state functions and collect results
		LimitStateFunction *theLimitStateFunction;
		LimitStateFunctionIter &lsfIter = theReliabilityDomain->getLimitStateFunctions();
		//for (i=0; i<numLsf; i++ ) {
		while ((theLimitStateFunction = lsfIter()) != 0) {
            int lsfTag = theLimitStateFunction->getTag();
            //int i = theLimitStateFunction->getIndex();
            int i = theReliabilityDomain->getLimitStateFunctionIndex(lsfTag);

            // get betas
            (*allBetas)(i) = theFunctionEvaluator->getResponseVariable("betaFORM", lsfTag);
            (*allPf1s)(i) = 1.0 - uRV.getCDFvalue( (*allBetas)(i) );
            
            // get alpha vectors
            for (int j = 0; j < numRV; j++) {
                RandomVariable *theRV = theReliabilityDomain->getRandomVariablePtrFromIndex(j);
                int rvTag = theRV->getTag();
                allAlphas(j,i) = theFunctionEvaluator->getResponseVariable("alphaFORM", lsfTag, rvTag);
            }
		}
		
		// Compute vector of 'rhos', that is, dot products between the alphas
		rhos = new Matrix(numLsf,numLsf);
		
		double dotProduct;
		for (int i = 0; i<numLsf; i++ ) {
			for (int j = i+1; j<numLsf; j++ ) {
				dotProduct = 0.0;
				for (int k = 0; k < numRV; k++ ) {
					dotProduct += allAlphas(k,i)*allAlphas(k,j);
				}
				(*rhos)(i,j) = dotProduct;
				(*rhos)(j,i) = dotProduct;
			}
		}
	}

	//opserr << "B vector:" << *allBetas;
	//opserr << "R matrix:" << *rhos;
	
	// Compute the bi-variate parallel probability for all the pairs
	// Note this is upper-diagonal only
	Pmn = new Matrix(numLsf,numLsf);
	CorrelatedStandardNormal phi2(0.0);
	for (int m = 0; m<numLsf; m++ ) {
		for (int n = m+1; n<numLsf; n++) {
			phi2.setCorrelation((*rhos)(m,n));
			(*Pmn)(m,n) = phi2.getCDF(-(*allBetas)(m),-(*allBetas)(n));
		}
	}
	
	return 0;
}


int
SystemAnalysis::computeBounds(int aType)
{	
	double estimate, lowerBound, upperBound;
	int i, j, result;
			
	if (aType == 0) {
		// parallel system
		
		// uni-component lower bound from Boole parallel system
		estimate = 0;
		for (i=0; i<numLsf; i++)
			estimate += (*allPf1s)(i);
		estimate -= numLsf-1;
		minLowerBound = estimate;

		// uni-component upper bound from Boole parallel system
		upperBound = (*allPf1s)(0);
		for (i=1; i<numLsf; i++) {
			if ((*allPf1s)(i) < upperBound)
				upperBound = (*allPf1s)(i);
		}
		maxUpperBound = upperBound;

		if (minLowerBound < 0.0)
			minLowerBound = 0;
		if (maxUpperBound > 1.0)
			maxUpperBound = 1;
		
	} else if (aType == 1) {
		// series system

		// use simulation approach to limit factorial possibilities for the ordering
		RandomNumberGenerator *theRandomNumberGenerator = new CStdLibRandGenerator();
		result = theRandomNumberGenerator->generate_nIndependentUniformNumbers(numLsf,0,numLsf-1,time(NULL));
		if (result < 0)
			opserr << "SystemAnalysis::computeBounds() - could not generate random numbers for simulation." << endln;

		long int numPerms = factorial(numLsf+1);
		int multiplicationFactor = 500;
		if (numPerms > multiplicationFactor*numLsf)
			numPerms = multiplicationFactor*numLsf;
		ID permute(numLsf);
		
		for (long int arr = 1; arr <= numPerms; arr++ ) {
			arrange(numLsf,theRandomNumberGenerator,permute);
			
			// bi-component lower bound from KHD series system
			lowerBound = (*allPf1s)(permute(0));
			for (i=2; i <= numLsf; i++) {
				estimate = (*allPf1s)(permute(i-1));
				
				for (j=1; j <= i-1; j++)
					estimate -= (*Pmn)(permute(i-1),permute(j-1));
				
				if (estimate > 0)
					lowerBound += estimate;
			}

			// bi-component upper bound from KHD series system
			upperBound = (*allPf1s)(permute(0));
			for (i=2; i<= numLsf; i++) {
				estimate = 0;
				
				for (j=1; j < i; j++) {
					if ((*Pmn)(permute(i-1),permute(j-1)) > estimate)
						estimate = (*Pmn)(permute(i-1),permute(j-1));
				}
				
				upperBound += (*allPf1s)(permute(i-1)) - estimate;
			}

			// Update bounds
			if (lowerBound < minLowerBound && lowerBound >= 0.0)
				minLowerBound = lowerBound;
			if (upperBound > maxUpperBound && upperBound <= 1.0)
				maxUpperBound = upperBound;
		}
		
		delete theRandomNumberGenerator;
	
	} else {
	
	}
	
	return 0;
		
}


int
SystemAnalysis::setCutsets(void)
{
	Cutset *theCutset;
	CutsetIter &cutIter = theReliabilityDomain->getCutsets();
	while ((theCutset = cutIter()) != 0) {
		const Vector &allComps = theCutset->getComponents();
		int nc = theCutset->getNumberOfComponents();
		Vector cutBeta(nc);
		Matrix cutRho(nc,nc);
		
		// remember that the cutsets can refer to LSF numbers only because the use of cutset in domain 
		// checks for corresponding LSF when created. However, allow computation of upper bound using 
		// both methods here.  
		for (int i = 0; i < nc; i++) {
			int actual = 0;
			
			if (strlen(betaFile) > 1 && strlen(rhoFile) > 1) {
				// assume LSF are numbered 1:numLsf
				actual = allComps(i);
				cutBeta(i) = (*allBetas)(actual-1);
				for (int j = 0; j < nc; j++)
					cutRho(j,i) = (*rhos)(allComps(j),actual-1);
			}
			else {
				// get actual LSF index from tag
				int iComp = allComps(i);
				actual = theReliabilityDomain->getLimitStateFunctionIndex( abs(iComp) );
				cutBeta(i) = (*allBetas)(actual) * sign(iComp);
					
				for (int j = 0; j < nc; j++) {
					int jComp = allComps(j);
					int actualj = theReliabilityDomain->getLimitStateFunctionIndex( abs(jComp) );
					cutRho(j,i) = (*rhos)(actualj,actual) * sign(iComp) * sign(jComp);
				}
			}
		}
		
		theCutset->setBetaCutset(cutBeta);
		theCutset->setRhoCutset(cutRho);
		
		//opserr << cutBeta;
		//opserr << cutRho;
	}
	
	return 0;
}


double 
SystemAnalysis::getLowerBound(void)
{
	return minLowerBound;
}


double 
SystemAnalysis::getUpperBound(void)
{
	return maxUpperBound;
}

int 
SystemAnalysis::getNumberOfLimitStateFunctions(void)
{
	// note that since you can read beta and rho from a file, this should not simply
	// query the reliability domain
	return numLsf;
}

const Vector &
SystemAnalysis::getBeta(void)
{
	return *allBetas;
}

const Matrix &
SystemAnalysis::getRho(void)
{
	return *rhos;
}

long int
SystemAnalysis::factorial(int num)
{
	if (num == 0)
		return 1;
	else if (num < 0)
		return -1;
	
	int i = num-1;
	long int result = num;

	while (i > 0) {
		result = result * i;
		i--;
	}

	return result;
}

int 
SystemAnalysis::sign(int v)
{
	return v > 0 ? 1 : (v < 0 ? -1 : 0);
}

int
SystemAnalysis::arrange(int num, RandomNumberGenerator *randNum, ID &permutation)
{
	// randomly arranges the input pf and Pmn matrix entries
	int result;
	
	Vector index(num);
	
	randNum->generate_nIndependentUniformNumbers(num,0,num-1);
	const Vector &randomArray = randNum->getGeneratedNumbers();

	for (int i=0; i<num; i++) {
	  result = int(randomArray(i) + 0.5); // rounds to nearest int
		if (index(result) > 0) {
			// repeat found
			for (int j=0; j<num; j++) {
				if (index(j) == 0.0) {
					result = j;
					break;
				}
			}
		} 
		// unique assignment
		permutation(i) = result;
		index(result) = 1;
	}
	
	return 0;
}

int
SystemAnalysis::getNumPermutations(int k, int n)
{
	if (k > n) {
		opserr << "k must be less than n for n choose k permutations" << endln;
		return -1;
	}
	// binomial
	return factorial(n)/factorial(k)/factorial(n-k);
}

int 
SystemAnalysis::setPermutations(int k, int n)
{
	// early exit if user specified rho and beta files (hence no cutsets)
	if (strlen(betaFile) > 1 && strlen(rhoFile) > 1)
		return -1;
		
	// create arrays of vectors and matrices that contain beta and rho for each permutation
	int perms = getNumPermutations(k,n);
	//opserr << "set " << n << " given " << k << " permutations = " << perms << endln;
	
	// vector 1:n of cut set indices
	Vector v(n);
	for (int i = 1; i <= n; i++)
		v(i-1) = i;
	
	// put realizations of k possible permutations of n into matrix
	sets->resize(perms,k);
	sets->Zero();
	Vector tmp(k);
	int is, js;
	
	// generate actual combinations
	if (n == k) {
		for (js = 0; js < k; js++)
			(*sets)(0,js) = v(js);
	}
	else if (n == k + 1) {
		for (is = 0; is < perms; is++) {
			tmp.Zero();
			int tmp_indx = 0;
			for (int j = 1; j <= n; j++) {
				if (j != n-is) {
					tmp(tmp_indx) = v(j-1);
					tmp_indx++;
				}
			}
			for (js = 0; js < k; js++)
				(*sets)(is,js) = tmp(js);
		}
	}
	else if (k == 1) {
		for (is = 0; is < perms; is++) {
			tmp(0) = v(is);
			(*sets)(is,0) = tmp(0);
		}
	}
	else if (n < 17 && (k > 3 || n-k < 4)) {
		long int rows = (long int)pow(2.0,n);
		long int ncycles = rows;
		long int nreps, im, jm, cnter;
		int count;

		Matrix x(rows,n);
		for (count = 1; count <= n; count++) {
			ncycles = ncycles/2;
			nreps = rows/(2*ncycles);
			Matrix settings(2*nreps,ncycles);
			for (im = 0; im < nreps; im++) {
				for (jm = 0; jm < ncycles; jm++)
					settings(im,jm) = 1;
			}
			
			// put settings matrix column by column into x matrix
			cnter = 0;
			for (jm = 0; jm < ncycles; jm++) {
				for (im = 0; im < 2*nreps; im++) {
					x(cnter,n-count) = settings(im,jm);
					cnter++;
				}
			}
		}
		
		Vector unitvec(n);
		Vector rowsum(rows);
		rowsum.addMatrixVector(0,x,unitvec+1,1);
		int nrows = 0;
		for (im = 0; im < rows; im++) {
			if (rowsum(im) == k)
				nrows++;
		}
		
		// note that nrows better equal perms
		Matrix idx(nrows,n);
		cnter = 0;
		for (im = 0; im < rows; im++) {
			if (rowsum(im) == k) {
				for (jm = 0; jm < n; jm++)
					idx(cnter,jm) = x(im,jm);
				cnter++;
			}
		}
		
		// generate actual combinations in tmp vector and store in sets
		for (is = 0; is < perms; is++) {
			tmp.Zero();
			int tmp_indx = 0;
			for (jm = 0; jm < n; jm++) {
				if ( idx(is,jm) == 1 ) {
					tmp(tmp_indx) = v(jm);
					tmp_indx++;
				}
			}
			for (js = 0; js < k; js++)
				(*sets)(is,js) = tmp(js);
		}
	}
	else {
		int ks;
		is = 0;
		
		for (int idx = 1; idx <= n-k+1; idx++) {
			Vector vpass(n-idx);
			Matrix Q( getNumPermutations(k-1, n-idx), k-1 );
			for (js = 0; js < n-idx; js++)
				vpass(js) = v(idx+js);
				
			combinations( vpass, k-1, Q );
			for (ks = 0; ks < getNumPermutations(k-1, n-idx); ks++) {
				(*sets)(is,0) = v(idx-1);
				for (js = 0; js < k-1; js++)
					(*sets)(is,js+1) = Q(ks,js);
					
				is++;
			}
		}
	}
	//opserr << *sets;
	
	return 0;
}

int 
SystemAnalysis::combinations(Vector &v, int m, Matrix &Q)
{
	int n = v.Size();
	int is = 0, js, ks;
	Q.Zero();
	
	if (n == m) {
		for (js = 0; js < n; js++)
			Q(0,js) = v(js);
			
	} else if (m == 1) {
		for (is = 0; is < n; is++)
			Q(is,0) = v(is);
			
	} else {
		if (m < n && m > 1) {
			for (int k = 1; k <= n-m+1; k++) {
				Vector newv(n-k);
				Matrix newQ( getNumPermutations(m-1, n-k), m-1 );
				for (int nv = 0; nv < n-k; nv++)
					newv(nv) = v(nv + k);
				
				combinations( newv,m-1,newQ );
				for (ks = 0; ks < getNumPermutations(m-1, n-k); ks++) {
					Q(is,0) = v(k-1);
					for (js = 0; js < m-1; js++)
						Q(is,js+1) = newQ(ks,js);
						
					is++;
				}
			}
		}
	}
	
	return 0;
}

int 
SystemAnalysis::setPermutedComponents(int k, int i)
{
	Cutset *theCutset;
	Vector single_set(k);
	int js;
	for (js = 0; js < k; js++)
		single_set(js) = (*sets)(i,js);

	// need to determine size of new beta and rho, taking unique components ONLY
	int cutLen = 0;
	for (js = 0; js < k; js++) {
		theCutset = theReliabilityDomain->getCutsetPtrFromIndex( single_set(js)-1 );
		if (theCutset == 0) {
			opserr << "error getting theCutset" << endln;
			return -1;
		}
		else
			cutLen += theCutset->getNumberOfComponents();
	}
	
	// fill out component numbers to determine which ones are unique
	Vector tmp_comp(cutLen);
	int curi = 0;
	for (js = 0; js < k; js++) {
		theCutset = theReliabilityDomain->getCutsetPtrFromIndex( single_set(js)-1 );
		const Vector &allComps = theCutset->getComponents();
		for (int jk = 0; jk < theCutset->getNumberOfComponents(); jk++) {
			int isunique = 1;
			for (int kl = curi; kl >= 0; kl--) {
				if ( tmp_comp(kl) == allComps(jk) )
					isunique = 0;
			}
			
			if ( isunique == 1 ) {
				tmp_comp(curi) = allComps(jk);
				curi++;
			}
		}
	}
	cutLen = curi;
	
	Vector permutedComps(cutLen);
	for (js = 0; js < cutLen; js++)
		permutedComps(js) = tmp_comp(js);
		
	//opserr << permutedComps;
	
	permutedBetas->resize(cutLen);
	permutedBetas->Zero();
	permutedRhos->resize(cutLen,cutLen);
	permutedRhos->Zero();
	
	// now flesh out vectors with data from each component in cutset
	for (int j = 0; j < cutLen; j++) {
		int actual = theReliabilityDomain->getLimitStateFunctionIndex( fabs(permutedComps(j)) );
		(*permutedBetas)(j) = (*allBetas)(actual) * sign(permutedComps(j));
		for (int jk = 0; jk < cutLen; jk++) {
			int actualj = theReliabilityDomain->getLimitStateFunctionIndex( fabs(permutedComps(jk)) );
			(*permutedRhos)(jk,j) = (*rhos)(actualj,actual) * sign(permutedComps(j)) * sign(permutedComps(jk));
		}
	}
		
	return 0;
}

const Vector &
SystemAnalysis::getBetaPermutation()
{
	//opserr << *permutedBetas;
	return *permutedBetas;
}

const Matrix &
SystemAnalysis::getRhoPermutation()
{
	//opserr << *permutedRhos;
	return *permutedRhos;
}
