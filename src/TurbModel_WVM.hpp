#include <iostream>
#include <cstdlib>
#include <string.h>
#include <ctime>
#include <random>
#include "mathoper.hpp"
using namespace std;

class TurbModel_WVM
{
  // Constructors & Destructor
public:
  TurbModel_WVM(const char *name);
  ~TurbModel_WVM();

  // Data memebers
public:
  int nshells; // number of shells
  int nmodes;  // number of modes per shell
  int npar;    // number of particles
  
  double (*e)[3];
  double (*a)[3];
  double (*u)[3];

  default_random_engine generator;
  normal_distribution<double> distribution;

  // Memeber functions
public:
  void initialHookScalarRansTurbModel(double *rey);

  void bkeuler(int ipar, double dt, double (*Gn)[3]);
  
  void rk4(int ipar, double dt, double (*Gn)[3], double (*Gnph)[3], 
           double (*Gnp1)[3]);
  
  void calcRhs(double *erhs, double *arhs, double *urhs, double (*G)[3],
	       double *eref, double *aref, double *uref);

  void calcReStress(double *rey, double (*Gn)[3], double (*Gnph)[3], 
                    double (*Gnp1)[3], double dt);
};
