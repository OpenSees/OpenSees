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


// $Revision: 1.1 $
// $Date: 2008-02-29 19:43:53 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/telm/InitialStaticAnalysis.cpp,v $
                                                                     
#include <InitialStaticAnalysis.h>
#include <Node.h>
#include <NodeIter.h>
#include <fstream>
#include <iomanip>
#include <iostream>
using std::ifstream;
using std::ios;
using std::setw;
using std::setprecision;

InitialStaticAnalysis::InitialStaticAnalysis
					  (ReliabilityDomain* passedReliabilityDomain,
					   Domain* passedDomain,
					   bool passedprint)
{
	theReliabilityDomain = passedReliabilityDomain;
	theDomain=passedDomain;
	print=passedprint;
	if(print){
		output.open("InitialStaticAnalysis.txt", ios::out);
	}
}
InitialStaticAnalysis::~InitialStaticAnalysis()
{
}
void InitialStaticAnalysis::printResult()
{
	output<<"\n";
	output<<" ---- result of initial static analysis -----\n";
	output<<"\n";
	output.setf(ios::right);
	output.setf(ios::scientific, ios::floatfield);
	NodeIter& theNodes = theDomain->getNodes();
	Node* theNode;
	while((theNode = theNodes()) != 0){
		int tag=theNode->getTag();
		int ndof=theNode->getNumberDOF();
		Vector Disp=theNode->getDisp();
		output << setw(10) << tag;
		for(int i=0; i<ndof; i++)
			output << setw(15) << setprecision(5) << Disp(i);
		output<<"\n";
	}
	output.flush();
}

