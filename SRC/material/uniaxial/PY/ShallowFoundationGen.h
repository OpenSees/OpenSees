//////////////////////////////////////////////////////////////////////
//																	//
//																	//
//  Written by: Prishati Raychowdhury (prishati@gmail.com)			//					
//              Graduate Student, UC San Diego						//
//              November, 2007										//
//////////////////////////////////////////////////////////////////////
//$Revision: 1.1 $
//$Date: 2009-01-07 22:43:24 $
//$Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/PY/ShallowFoundationGen.h,v $

#include <fstream>
#include <cmath>
#include <iostream>
#include <ID.h>

class ShallowFoundationGen
{	
	
public:

	// Public member functions accessed from TclModelBuilder.cpp
	void GetShallowFoundation(const char *FoundationID, const char *ConnectingNode, const char *InputDataFile, const char *FoundationMatType);

	ShallowFoundationGen();
	~ShallowFoundationGen();
};
