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
// $Date: 2017-09-15 10:45:28 $
// $Source: /usr/local/cvs/OpenSees/SRC/recorder/MPCORecorder.h,v $


// Written: Massimo Petracca 
// Created: 01/01
// Revision: A
//
// Description: This file contains the class definition for MPCORecorder.

/*
MPCO recorder - Written by ASDEA Software Technology: M.Petracca, G.Camata\n
ASDEA Software Technology: https://asdeasoft.net
STKO (Scientific ToolKit for OpenSees): https://asdeasoft.net/stko/
If you use this tool, please cite us:
Petracca, M., Candeloro, F., & Camata, G. (2017). \"STKO user manual\". ASDEA Software Technology, Pescara Italy.

The MPCORecorder writes data to an HDF5 database.

MPCORecorder uses the HDF5 library version 1.10 from the HDF Group:
The HDF Group, "Hierarchical Data Format, version 5," 1997-2017. [Online]. Available: http://www.hdfgroup.org/HDF5/.

The HDF5 shared libary is loaded at run-time and should be available on your machine if you use this
recorder.
For more information about HDF5 library, licensing terms and download:
http://www.hdfgroup.org
http://www.hdfgroup.org/HDF5/
https://www.hdfgroup.org/downloads/hdf5/
*/

#ifndef MPCORecorder_h
#define MPCORecorder_h

#include <Recorder.h>

class MPCORecorder : public Recorder
{
	friend void* OPS_MPCORecorder();
private:
	class private_data;

public:
	MPCORecorder();
	~MPCORecorder();
	int record(int commitTag, double timeStamp);
	virtual int restart(void);
	virtual int domainChanged(void);
	virtual int setDomain(Domain &theDomain);
	virtual int sendSelf(int commitTag, Channel &theChannel);
	virtual int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

private:
	int initialize();

	int writeModel();
	int writeModelNodes();
	int writeModelElements();
	int writeModelLocalAxes();
	int writeSections();
	int writeSets();

	int initNodeRecorders();
	int clearNodeRecorders();
	int initElementRecorders();
	int clearElementRecorders();

	int recordResultsOnNodes();
	int recordResultsOnElements();

protected:

private:
	private_data *m_data;
};


#endif // MPCORecorder_h