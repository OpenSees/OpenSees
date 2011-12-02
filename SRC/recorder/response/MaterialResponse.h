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
// $Date: 2000-12-18 11:35:47 $
// $Source: /usr/local/cvs/OpenSees/SRC/recorder/response/MaterialResponse.h,v $
                                                                        
// Written: MHS 
// Created: Oct 2000
//
// Description: This file contains the MaterialResponse class interface

#ifndef MaterialResponse_h
#define MaterialResponse_h

#include <Response.h>
#include <Information.h>

class Material;

class ID;
class Vector;
class Matrix;
class Tensor;

class MaterialResponse : public Response
{
public:
	MaterialResponse(Material *mat, int id);
	MaterialResponse(Material *mat, int id, int val);
	MaterialResponse(Material *mat, int id, double val);
	MaterialResponse(Material *mat, int id, const ID &val);
	MaterialResponse(Material *mat, int id, const Vector &val);
	MaterialResponse(Material *mat, int id, const Matrix &val);
	MaterialResponse(Material *mat, int ID, const Tensor &val);
	~MaterialResponse();

	int getResponse(void);
	void Print(ostream &s, int flag = 0);

private:
	Material *theMaterial;
	int responseID;
	Information matInfo;
};

#endif
