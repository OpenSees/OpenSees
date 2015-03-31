#ifndef _OPS_TriDiagonalMatrixF
#define _OPS_TriDiagonalMatrixF
class TriDiagonalMatrixF
{		
	public:
		int Lenght;
		/// <summary>
		/// The values for the sub-diagonal. A[0] is never used.
		/// </summary>
		double *A;

		/// <summary>
		/// The values for the main diagonal.
		/// </summary>
		double *B;

		/// <summary>
		/// The values for the super-diagonal. C[C.Length-1] is never used.
		/// </summary>
		double *C;

		int N();

		~TriDiagonalMatrixF(void);
		
		/// <summary>
		/// Construct an NxN matrix.
		/// </summary>
		TriDiagonalMatrixF(int n);
		void SetMat(int row, int col, double value);
		double GetMat(int row, int col);
		double* Solve(double* d, int dLength);
};
#endif

