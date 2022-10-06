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

// $Revision: 1.0 $
// $Date: 2020-11-27 18:07:11 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/analysis/ResponseSpectrumAnalysis.h,v $

// Written: Massimo Petracca (ASDEA Software) 
// Created: Fri Nov 27 18:07:11: 2020
// Revision: A
//
// Description: This file contains the class definition for ResponseSpectrumAnalysis.
//
// What: "@(#) ResponseSpectrumAnalysis.h, revA"

#ifndef ResponseSpectrumAnalysis_h
#define ResponseSpectrumAnalysis_h

#include <vector>
class AnalysisModel;
class TimeSeries;

class ResponseSpectrumAnalysis
{
public:
	ResponseSpectrumAnalysis(
		AnalysisModel* theModel,
		TimeSeries* theFunction,
		const std::vector<double>& Tn,
		const std::vector<double>& Sa,
		int theDirection,
		double scale
	);
	~ResponseSpectrumAnalysis();

public:
	int analyze();
	int analyze(int mode_id);

private:
	int check();
	int beginMode();
	int endMode();
	int solveMode();
	double getSa(double T) const;

private:
	// the model
	AnalysisModel* m_model;
	// the response spectrum function
	TimeSeries* m_function;
	std::vector<double> m_Tn;
	std::vector<double> m_Sa;
	// the direction 1 to 3 (for 2D models) or 1 to 6 (for 3D models)
	int m_direction;
	// the scale factor for the computed displacement field
	double m_scale;
	// current mode
	int m_current_mode;
};

#endif
