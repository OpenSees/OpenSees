// Calculate Shape functions and parametric derivatives for a
//   3-D finite element with varied nen between 8 and 27.

// By Jinchi Lu, 04/15/2004
// (Based on a Fortran code from R. Brockman 30 May 1986)

#include <stdio.h>
#include <math.h>

void shap3dv(double *R, int *NP, double Q[27][4]){

// Formal parameters:
//
//  R (Input) = Parametric coordinates, in the interval [-1,1],
// at which shape functions are to be evaluated
//  NP (Input) = List of flags for each local node of the finite
// element. (Usually this is just the connectivity list)
// = 0 if node is absent
// > 0 if node is present
// * Q (Output) = Shape functions and derivatives, as follows:
// * Q[i][3] is the shape function for local node 'i'
// * Q[i][0] is the r-derivative
// * Q[i][1] is the s-derivative
// * Q[i][2] is the t-derivative
//
// Local Node Pattern for Three-Dimensional Elements:
//
// * Nodes 1 - 4 Lower surface, counterclockwise
// * Nodes 5 - 8 Upper surface, counterclockwise
// * Nodes 9 - 12 Midsides of edges 1-2, 2-3, 3-4, 4-1
// * Nodes 13 - 16 Midsides of edges 5-6, 6-7, 7-8, 8-5
// * Nodes 17 - 20 Midsides of edges 1-5, 2-6, 3-7, 4-8
// * Nodes 21 - 26 Mid-face nodes on +r, +s, +t, -r, -s, -t
// * Node 27 Centroid node
//

	double G[3][3], D[3][3], C;
	int L[27] = {3,1,1,3,3,1,1,3,2,1,2,3,2,1,2,3,3,1,1,3,1,2,2,3,2,2,2},
	    M[27] = {3,3,1,1,3,3,1,1,3,2,1,2,3,2,1,2,3,3,1,1,2,1,2,2,3,2,2},
	    N[27] = {3,3,3,3,1,1,1,1,3,3,3,3,1,1,1,1,2,2,2,2,2,2,1,2,2,3,2};

	int i, j, LR, LS, LT;
	for( i = 0; i < 3; i++ ) {
         G[0][i] = 0.5 + 0.5*R[i];
         G[1][i] = 1.0 - R[i]*R[i];
         G[2][i] = 0.5 - 0.5*R[i];
         D[0][i] =  0.5	;
         D[1][i] = -2.0*R[i];
         D[2][i] = -0.5;
	}
	
	// Construct basic three-dimensional quadratic shape functions 

	for( i = 0; i < 27; i++ ) {
         LR = L[i]-1; 
         LS = M[i]-1;
         LT = N[i]-1;
         Q[i][0] = D[LR][0] * G[LS][1] * G[LT][2]; 
         Q[i][1] = G[LR][0] * D[LS][1] * G[LT][2];
         Q[i][2] = G[LR][0] * G[LS][1] * D[LT][2]; 
         Q[i][3] = G[LR][0] * G[LS][1] * G[LT][2];
//		 printf("%d %15.6e %15.6e %15.6e %15.6e\n", i, Q[i][0], Q[i][1], Q[i][2], Q[i][3]);
	}

	//   Modify basic shape functions to account for omitted nodes 

	for(j = 0; j < 4; j++ ) {
		if( NP[26] == 0 ) Q[26][j] = 0.;	
        C = -0.5*Q[26][j]; 
		for( i = 20; i < 26; i++ ) {
            Q[i][j] +=  C ;
			if( NP[i] == 0 ) Q[i][j] = 0.;
		}

		C = 0.5*C; 
		Q[ 8][j] +=  - 0.5* (Q[25][j] + Q[24][j]) + C; 
		Q[ 9][j] +=  - 0.5* (Q[25][j] + Q[20][j]) + C; 
		Q[10][j] +=  - 0.5* (Q[25][j] + Q[21][j]) + C; 
		Q[11][j] +=  - 0.5* (Q[25][j] + Q[23][j]) + C; 
		Q[12][j] +=  - 0.5* (Q[24][j] + Q[22][j]) + C; 
		Q[13][j] +=  - 0.5* (Q[20][j] + Q[22][j]) + C; 
		Q[14][j] +=  - 0.5* (Q[21][j] + Q[22][j]) + C; 
		Q[15][j] +=  - 0.5* (Q[23][j] + Q[22][j]) + C; 
		Q[16][j] +=  - 0.5* (Q[23][j] + Q[24][j]) + C; 
		Q[17][j] +=  - 0.5* (Q[24][j] + Q[20][j]) + C; 
		Q[18][j] +=  - 0.5* (Q[20][j] + Q[21][j]) + C; 
		Q[19][j] +=  - 0.5* (Q[21][j] + Q[23][j]) + C; 

		for( i = 8; i < 20; i++ ) {
			if( NP[i] == 0 ) Q[i][j] = 0.;
		}
 
        C = 0.5*C;
		for( i = 0; i < 8; i++ ) {
			Q[i][j] += C;
		}

         Q[0][j] +=  - 0.50* (Q[16][j] + Q[11][j] + Q[ 8][j]) - 0.25* (Q[23][j] + Q[24][j] + Q[25][j]); 
         Q[1][j] +=  - 0.50* (Q[17][j] + Q[ 9][j] + Q[ 8][j]) - 0.25* (Q[20][j] + Q[24][j] + Q[25][j]); 
         Q[2][j] +=  - 0.50* (Q[10][j] + Q[ 9][j] + Q[18][j]) - 0.25* (Q[25][j] + Q[20][j] + Q[21][j]); 
         Q[3][j] +=  - 0.50* (Q[19][j] + Q[11][j] + Q[10][j]) - 0.25* (Q[25][j] + Q[21][j] + Q[23][j]); 
         Q[4][j] +=  - 0.50* (Q[16][j] + Q[15][j] + Q[12][j]) - 0.25* (Q[22][j] + Q[23][j] + Q[24][j]); 
         Q[5][j] +=  - 0.50* (Q[17][j] + Q[12][j] + Q[13][j]) - 0.25* (Q[24][j] + Q[20][j] + Q[22][j]); 
         Q[6][j] +=  - 0.50* (Q[18][j] + Q[13][j] + Q[14][j]) - 0.25* (Q[20][j] + Q[21][j] + Q[22][j]); 
         Q[7][j] +=  - 0.50* (Q[19][j] + Q[15][j] + Q[14][j]) - 0.25* (Q[21][j] + Q[22][j] + Q[23][j]); 
	}

}

int brcshl(double shl[4][20][27], double w[27], int nint, int nen) {
/*
     PROGRAM TO CALCULATE INTEGRATION-RULE WEIGHTS, SHAPE FUNCTIONS
        AND LOCAL DERIVATIVES FOR A EIGHT-NODE BRICK ELEMENT

             R,S,T = LOCAL ELEMENT COORD ("XI", "ETA", "ZETA" RESP.)
        SHL(1,I,L) = LOCAL ("XI") DERIVATIVE OF SHAPE FUNCTION
        SHL(2,I,L) = LOCAL ("ETA") DERIVATIVE OF SHAPE FUNCTION
        SHL(3,I,L) = LOCAL ("ZETA") DERIVATIVE OF SHAPE FUNCTION
        SHL(4,I,L) = LOCAL  SHAPE FUNCTION
              W(L) = INTEGRATION-RULE WEIGHT
                 I = LOCAL NODE NUMBER
                 L = INTEGRATION POINT NUMBER
              NINT = NUMBER OF INTEGRATION POINTS, EQ. 1, 8 OR 27
*/

	double RA[27] = {-0.50, 0.50, 0.50,-0.50,-0.50, 0.50, 0.50,-0.50,
     		 0.00, 0.50, 0.00,-0.50, 0.00, 0.50, 0.00,-0.50,
     		-0.50, 0.50, 0.50,-0.50, 0.50, 0.00, 0.00,-0.50,
     		 0.00, 0.00, 0.00};

	double SA[27] = {-0.50,-0.50, 0.50, 0.50,-0.50,-0.50, 0.50, 0.50,
     		-0.50, 0.00, 0.50, 0.00,-0.50, 0.00, 0.50, 0.00,
     		-0.50,-0.50, 0.50, 0.50, 0.00, 0.50, 0.00, 0.00,
     		-0.50, 0.00, 0.00};     

	double TA[27] = {-0.50,-0.50,-0.50,-0.50, 0.50, 0.50, 0.50, 0.50,
     		-0.50,-0.50,-0.50,-0.50, 0.50, 0.50, 0.50, 0.50,
     		 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.50, 0.00,
     		 0.00,-0.50, 0.00};		
	double five9 = 0.5555555555555556, eight9 = 0.8888888888888889;
	double G = 0., R[3], Q[27][4];
	int i, j, L, NP[27];
	w[0] = 8;
	if( nint == 8) {
		G = 2/sqrt(3.0);
		for(i = 0; i < nint; i++ )
			w[i] = 1.;
	}
	else if ( nint == 27 ) {
		G = 2.0*sqrt(3./5.);
		w[0] = five9 * five9 * five9;
		for( i = 1; i < 8; i++ ) 
			w[i] = w[0];
		w[8] = five9 * five9 * eight9;
		for( i = 9; i < 20; i++ )
			w[i] = w[8];
		w[20] = five9 * eight9 * eight9;
		for( i = 21; i < 26; i++ ) 
			w[i] = w[20];
		w[26] = eight9 * eight9 * eight9;
	}
	else {
		//printf("invalid nint %d\n", nint);
		//exit(-1);
		return -1;
	}

	for( i = 0; i < 27; i ++ )
		NP[i] = 1;
	if( nen < 27) {
		for( i = nen; i < 27; i++ ) 
			NP[i] = 0;
	}
	else if( nen < 8 ) {
		//printf("invalid nen %d\n", nen);
		//exit(-1);
		return -1;
	}

	for( L = 0; L < nint; L++ ) {
		R[0] = G * RA[L];
		R[1] = G * SA[L];
		R[2] = G * TA[L];

		shap3dv(R, NP, Q);

		for( j = 0; j < nen; j++) {
			for( i = 0; i < 4; i++ ) {
				shl[i][j][L] = Q[j][i];
			}
		//printf("%5d %5d %15.6e %15.6e %15.6e %15.6e\n", L+1, j+1,
		//	shl[0][j][L],shl[1][j][L],shl[2][j][L],shl[3][j][L]);
		}

	}
	return 0;
}





















