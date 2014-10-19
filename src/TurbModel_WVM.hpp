#include <iostream>
#include <cstdlib>
#include <string.h>
#include <ctime>
#include "mathoper.hpp"
using namespace std;

class prm
{
  // Constructors & Destructor
public: 
  prm();
  ~prm();

  // Data memebers
protected:
  double e[3];   // unit wave vector
  double a[3];   // eddy-axis vector
  double u[3];   // velocity vector

  // Member functions
public:
  void get_eau(double *eref, double *aref, double *uref);
  void set_eau(double *eref, double *aref, double *uref);
  void calcRhs(double *erhs, double *arhs, double *urhs, double (*G)[3], 
	       double *eref, double *aref, double *uref);
  void update(double dt, double (*Gn)[3], double (*Gnph)[3], double (*Gnp1)[3]);
};

class TurbModel_WVM
{
  // Constructors & Destructor
public:
  TurbModel_WVM();
  ~TurbModel_WVM();

  // Data memebers
public:
  int nshells;    // number of shells
  int nmodes;  // number of modes per shell

  prm (*par)[1];     // particles

  // Memeber functions
public:
  void initialHookScalarRansTurbModel(double *rey);
  void calcReStress(double *rey, double (*Gn)[3], double (*Gnph)[3], 
                    double (*Gnp1)[3], double dt);
};
