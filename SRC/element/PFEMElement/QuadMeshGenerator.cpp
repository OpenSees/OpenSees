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
                                                                        
//
// Written: Minjie Zhu
//
// Description: A class for quad mesh generator
//

#include "QuadMeshGenerator.h"
#include <Matrix.h>


QuadMeshGenerator::QuadMeshGenerator()
    :point(), line(), pointout(), quadout()
{
}

QuadMeshGenerator::~QuadMeshGenerator()
{
}

// mesh
int
QuadMeshGenerator::mesh(double size)
{
    if (size <= 0) return -1;
    if (line.size() != 4) {
	opserr << "WARNING: must have four lines to mesh quad\n";
	return -1;
    }

    // copy point to pointout
    clearOutput();
    pointout = point;

    // the lines are assumed to be continues in a loop
    int m = line[0].Size()-1;
    int n = line[1].Size()-1;
    if (m != line[2].Size()-1) {
	opserr << "WARNING: opposite lines must have same number of points\n";
	return -1;
    }

    if (n != line[3].Size()-1) {
	opserr << "WARNING: opposite lines must have same number of points\n";
	return -1;
    }

    // order of the lines
    int order[4] = {0,0,0,0};
    int next[4] = {0,0,0,0};
    if (line[0][0] == line[1][0]) {
	order[0] = 1;
	order[1] = 0;
	next[0] = 0;
	next[1] = n;
    } else if (line[0][0] == line[1][n]) {
	order[0] = 1;
	order[1] = 1;
	next[0] = 0;
	next[1] = 0;
    } else if (line[0][m] == line[1][0]) {
	order[0] = 0;
	order[1] = 0;
	next[0] = m;
	next[1] = n;
    } else if (line[0][m] == line[1][n]) {
	order[0] = 0;
	order[1] = 1;
	next[0] = m;
	next[1] = 0;
    } else {
	opserr << "WARNING: line 0 and 1 are not connected\n";
	return -1;
    }

    if (line[1][next[1]] == line[2][0]) {
	order[2] = 0;
	next[2] = m;
    } else if (line[1][next[1]] == line[2][m]) {
	order[2] = 1;
	next[2] = 0;
    } else {
	opserr << "WARNING: line 1 and 2 are not connected\n";
	return -1;
    }

    if (line[2][next[2]] == line[3][0]) {
	order[3] = 0;
	next[3] = n;
    } else if (line[2][next[2]] == line[3][n]) {
	order[3] = 1;
	next[3] = 0;
    } else {
	opserr << "WARNING: line 2 and 3 are not connected\n";
	return -1;
    }

    if (order[0] == 0) {
	if (line[3][next[3]] != line[0][0]) {
	    opserr << "WARNING: line 0 and 3 are not connected -- func\n";
	    return -1;
	}
    } else {
	if (line[3][next[3]] != line[0][m]) {
	    opserr << "WARNING: line 0 and 3 are not connected -- func\n";
	    return -1;
	}
    }

    // create points
    if (m<=1 || n<=1) return 0;
    Matrix grid(m+1, n+1);

    for (int i=0; i<m+1; ++i) {
	int pt1 = line[2][m-i];
	int pt2 = line[0][i];
	if (order[2] == 1) {
	    pt1 = line[2][i];
	}
	if (order[0] == 1) {
	    pt2 = line[0][m-i];
	}
	
	for (int j=0; j<n+1; ++j) {
	    if (i == 0) {
		grid(i,j) = line[3][n-j];
		if (order[3] == 1) {
		    grid(i,j) = line[3][j];
		} 
	    } else if (i == m) {
		grid(i,j) = line[1][j];
		if (order[1] == 1) {
		    grid(i,j) = line[1][n-j];
		}
	    } else if (j == 0) {
		grid(i,j) = line[0][i];
		if (order[0] == 1) {
		    grid(i,j) = line[0][m-i];
		}
	    } else if (j == n) {
		grid(i,j) = line[2][m-i];
		if (order[2] == 1) {
		    grid(i,j) = line[2][i];
		}
	    } else {
		double N1 = (double)j/n;
		double N2 = 1.0 - N1;
		Vector crds = point[pt1];
		crds.addVector(N1, point[pt2], N2);

		grid(i,j) = (int)pointout.size();
		pointout.push_back(crds);
	    }
	    
	}
    }

    // create quads
    for (int i=0; i<m; ++i) {
	for (int j=0; j<n; ++j) {
	    ID pts(4);
	    pts(0) = grid(i,j);
	    pts(1) = grid(i+1,j);
	    pts(2) = grid(i+1,j+1);
	    pts(3) = grid(i,j+1);
	    quadout.push_back(pts);
	}
    }
    
    return 0;
}

// inputs
int
QuadMeshGenerator::addPoint(const Vector& crds)
{
    point.push_back(crds);
    return 0;
}

int
QuadMeshGenerator::addLine(const ID& pts)
{
    line.push_back(pts);
    return 0;
}

// outputs
int
QuadMeshGenerator::getNumPoints() const
{    
    return (int)pointout.size();
}

void
QuadMeshGenerator::getPoint(int i, Vector& crds)
{
    if (i>=0 && i<getNumPoints()) {
	crds = pointout[i];
    }
}

int
QuadMeshGenerator::getNumQuads() const
{
    return (int)quadout.size();
}

void
QuadMeshGenerator::getQuad(int i, ID& pts)
{
    if (i>=0 && i<getNumQuads()) {
	pts = quadout[i];
    }
}

void
QuadMeshGenerator::clear()
{
    point.clear();
    line.clear();
    clearOutput();
}

void
QuadMeshGenerator::clearOutput()
{
    pointout.clear();
    quadout.clear();
}
