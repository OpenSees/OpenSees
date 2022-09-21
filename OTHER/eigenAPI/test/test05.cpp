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
// testing: copyToVectorReference
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
void test_copyToMatrixReference(const Eigen::DenseBase<Derived> &m, Matrix &ops_mat)
{
	std::cout << std::endl;

	*opserrPtr << "OPS matrix before copy:" ;
	*opserrPtr << ops_mat ;

	std::cout << "Eigen matrix :\n";
	std::cout << m << std::endl;

	std::cout << std::endl;

	copyToMatrixReference(m, ops_mat);

	*opserrPtr << "OPS matrix after copy:" ;
	*opserrPtr << ops_mat ;

	return ;
}


int main()
{
	Matrix ops_m2x2(2,2);
	test_copyToMatrixReference(Matrix2d::Random().eval(), ops_m2x2);
	test_copyToMatrixReference(Matrix2d::Random().eval(), ops_m2x2);

	Matrix ops_m3x3(3,3);
	test_copyToMatrixReference(Matrix3d::Random().eval(), ops_m3x3);
	test_copyToMatrixReference(Matrix3d::Random().eval(), ops_m3x3);

	Matrix ops_m4x4(4,4);
	test_copyToMatrixReference(Matrix4d::Random().eval(), ops_m4x4);

	Matrix ops_m2x5(2,6);
	test_copyToMatrixReference(MatrixXd::Random(2,5).eval(), ops_m2x5);

}