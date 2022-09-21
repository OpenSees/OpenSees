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

// Jose Abell (UANDES, github.com/jaabell)
// Massimo Petracca - ASDEA Software, Italy (2022)
//
// Test interface between eigen and opensees matrices. 
//
// testing: copyToNewMatrix
//


#include <iostream>
#include "../EigenAPI.h"
#include <StandardStream.h>

StandardStream sserr;
OPS_Stream *opserrPtr = &sserr;


using Eigen::MatrixXd;
using Eigen::Matrix2d;
using Eigen::Matrix3d;
using Eigen::Matrix4d;

template <typename Derived>
void test_copyToNewMatrix(const Eigen::DenseBase<Derived> &m)
{
	std::cout << std::endl;

	std::cout << "Eigen matrix :\n";
	std::cout << m << std::endl;

	std::cout << std::endl;

	Matrix ops_mat = copyToNewMatrix(m);

	*opserrPtr << "OPS matrix:" ;
	*opserrPtr << ops_mat ;

	return ;
}


int main()
{
	test_copyToNewMatrix(Matrix2d::Random());
	test_copyToNewMatrix(Matrix3d::Random());
	test_copyToNewMatrix(Matrix4d::Random());
	test_copyToNewMatrix(MatrixXd::Random(3,8));  
	
	typedef Eigen::Matrix<double, 3, 8> myCustomMatrix;
	test_copyToNewMatrix(myCustomMatrix::Random());
}