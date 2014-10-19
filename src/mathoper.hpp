#ifndef MATHOPER_H
#define MATHOPER_H

#include <math.h>

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

#endif
