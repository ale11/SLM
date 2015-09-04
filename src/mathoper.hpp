#ifndef MATHOPER_H
#define MATHOPER_H

#include <math.h>
#include <armadillo>

inline double vecDotVec3d(const double *v1, const double *v2) {

  double dotproduct = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];

  return(dotproduct);
}

inline void normVec3d(double *v1)
{
  double mag = vecDotVec3d(v1, v1);

  mag = sqrt(mag);

  v1[0] /= mag;
  v1[1] /= mag;
  v1[2] /= mag;
}

inline void vecCrossVec3d(double *v1, double *v2, double *v3)
{
  v1[0] = v2[1]*v3[2] - v2[2]*v3[1];
  v1[1] = v2[2]*v3[0] - v2[0]*v3[2];
  v1[2] = v2[0]*v3[1] - v2[1]*v3[0];
}

inline double vecDotMatDotVec3d(double *v1, double *v2, double (*M)[3])
{
  double result = 0.0;

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      result += v1[i]*M[i][j]*v2[j];

  return result;
}

inline void matDotVec3d(double *v1, double *v2, double (*M)[3])
{
  for (int i = 0; i < 3; i++)
  {
    v1[i] = 0.0;
    for (int j = 0; j < 3; j++)
      v1[i] += M[i][j]*v2[j];
  }
}

inline void matTransDotVec3d(double *v1, double *v2, double (*M)[3])
{
  for (int i = 0; i < 3; i++)
    {
      v1[i] = 0.0;
      for (int j = 0; j < 3; j++)
	v1[i] += M[j][i]*v2[j];
    }
}

inline void matPlusMat3d(double (*M1)[3], double (*M2)[3], double (*M3)[3])
{
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      M1[i][j] = M2[i][j] + M3[i][j];
}

inline void matTimesMat3d(double (*M1)[3], double (*M2)[3], double (*M3)[3])
{
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
    {
      M1[i][j] = 0.0;
      for (int k = 0; k < 3; k++)
	M1[i][j] += M2[i][k]*M3[k][j];
    }
}

inline double traceMatTimesMat3d(double (*M1)[3], double (*M2)[3])
{
  double result = 0.0;

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      result += M1[i][j]*M2[j][i];
  
  return result;
}

inline void GaussSolver( const int n, double *M, double *v)
{
	// form lower triangular system
	for (int col = (n-1); col > 0; col--)
	{
		// check for singular matrix
		if ( M[col*n+col] != 0.0 )
		{
			for ( int row = 0; row < col; row++)
			{
				double ratio = M[row*n+col]/M[col*n+col];
				// update vector
				v[row] -= ratio*v[col];
				// update matrix
				for (int i = 0; i < col; i++)
					M[row*n+i] -= ratio*M[col*n+i];
			}
		}
		else std::cout << "Singular matrix" << std::endl;
	}

	// compute the back solution
	if ( M[0] != 0 )
	{
		v[0] /= M[0];
		for ( int row = 1; row < n; row++)
		{
			double result = v[row];
			for ( int col = 0; col < row; col++)
				result -= M[row*n+col]*v[col];
			v[row] = result/M[row*n+row];
		}
	}
	else std::cout << "Singular matrix" << std::endl;
}

inline void VecOrthoVec3d( double *v1, const double * v2)
{
  double v1Dotv2 = vecDotVec3d(v1, v2);
  for (int i = 0; i < 3; i++)
  	v1[i] = v1[i] - v1Dotv2*v2[i];
}

inline void SymMatOrthoVec3d( double M[6], const double *v)
{
	double A[3][3], b[3];

	A[0][0] = 1.0 + v[0]*v[0];   A[0][1] =     + v[0]*v[1];   A[0][2] =     + v[0]*v[2];
	A[1][0] =     + v[1]*v[0];   A[1][1] = 1.0 + v[1]*v[1];   A[1][2] =     + v[1]*v[2];
	A[2][0] =     + v[2]*v[0];   A[2][1] =     + v[2]*v[1];   A[2][2] = 1.0 + v[2]*v[2];

	b[0] = M[0]*v[0] + M[3]*v[1] + M[4]*v[2];
	b[1] = M[3]*v[0] + M[1]*v[1] + M[5]*v[2];
	b[2] = M[4]*v[0] + M[5]*v[1] + M[2]*v[2];

	GaussSolver(3, A[0], b);

	M[0] -= b[0]*v[0] + b[0]*v[0];
	M[1] -= b[1]*v[1] + b[1]*v[1];
	M[2] -= b[2]*v[2] + b[2]*v[2];
	M[3] -= b[0]*v[1] + b[1]*v[0];
	M[4] -= b[0]*v[2] + b[2]*v[0];
	M[5] -= b[1]*v[2] + b[2]*v[1];
}

inline void CholeskyDecomp(const int n, double *A, double *L)
{
	for (int i = 0; i < n; i++)
		for (int j = 0; j < (i+1); j++)
		{
			double s = 0;
			for (int k = 0; k < j; k++)
				s += L[i * n + k] * L[j * n + k];

			L[i * n + j] = (i == j) ?
					sqrt(A[i * n + i] - s) : (1.0 / L[j * n + j] * (A[i * n + j] - s));
		}
}

inline void CholeskySolver(const int n, double *M, double *v)
{
	double *L = new double [n*n];
	CholeskyDecomp(n,M,L);

	// compute back solution
	if ( L[0] != 0 )
	{
		v[0] /= L[0];
		for ( int row = 1; row < n; row++)
		{
			double result = v[row];
			for ( int col = 0; col < row; col++)
				result -= L[row*n+col]*v[col];
			v[row] = result/L[row*n+row];
		}
	}
	else std::cout << "Singular matrix" << std::endl;

	// compute forward solution
	if ( L[0] != 0 )
	{
		v[n-1] /= L[n*n-1];
		for ( int row = n-2; row >= 0; row--)
		{
			double result = v[row];
			for ( int col = (n-1); col > row; col--)
				result -= L[col*n+row]*v[col];
			v[row] = result/L[row*n+row];
		}
	}
	else std::cout << "Singular matrix" << std::endl;

}

/*// TO TEST THE GAUSS SOLVER
double dc[6], Ac[6][6];
Ac[0][0] = 4.0;  Ac[0][1] = 2.0;  Ac[0][2] = 1.0;
Ac[0][3] = 4.0;  Ac[0][4] = 5.0;  Ac[0][5] = 6.0;

Ac[1][0] = 8.0;  Ac[1][1] = 3.0;  Ac[1][2] = 10.0;
Ac[1][3] = 11.0; Ac[1][4] = 12.0; Ac[1][5] = 13.0;

Ac[2][0] = 14.0; Ac[2][1] = 7.0; Ac[2][2] = 16.0;
Ac[2][3] = 17.0; Ac[2][4] = 18.0; Ac[2][5] = 19.0;

Ac[3][0] = 13.0; Ac[3][1] = 11.0; Ac[3][2] = 2.0;
Ac[3][3] = 6.0; Ac[3][4] = 4.0; Ac[3][5] = 5.0;

Ac[4][0] = 2.0; Ac[4][1] = 7.0; Ac[4][2] = 4.0;
Ac[4][3] = 9.0; Ac[4][4] = 0.0; Ac[4][5] = 1.0;

Ac[5][0] = 3.0; Ac[5][1] = 4.0; Ac[5][2] = 5.0;
Ac[5][3] = 5.0; Ac[5][4] = 2.0; Ac[5][5] = 9.0;

dc[0] = 5.0; dc[1] = 8.0; dc[2] = 4.0; dc[3] = 3.0; dc[4] = 2.0; dc[5] = 1.0;

GaussSolver(6, Ac[0], dc);

cout << dc[0] << ", " << dc[1] << ", " << dc[2] << ", " << dc[3] << ", "
     << dc[4] << ", " << dc[5] << endl;*/

#endif
