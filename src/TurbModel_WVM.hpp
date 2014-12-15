#include "UgpWithCvCompFlow.hpp"
#include <random>

class TurbModel_WVM: public UgpWithCvCompFlow
{
  // Constructors & Destructor
public:
  TurbModel_WVM(const char *name);
  virtual ~TurbModel_WVM();

  // Data memebers
public:
  int nshells; // number of shells
  int nmodes;  // number of modes per shell
  int npar;    // number of particles
  
  double (*e)[3];
  double (*a)[3];
  double (*u)[3];

  double eps, tau;

  default_random_engine generator;
  normal_distribution<double> distribution;

  double Cn, Cv, Ceps1, Ceps2, Ceps3;

  // Memeber functions
public:
  void initialHookScalarRansTurbModel(double *struc, double &eps_init);

  void bkeuler(double dt, double (*Gn)[3]);
  
  void rk4(double dt, double (*Gn)[3], double (*Gnph)[3], double (*Gnp1)[3]);
  
  void calcRhs(double *erhs, double *arhs, double *urhs, double (*G)[3],
	       double *eref, double *aref, double *uref, double dt);

  double updateDissipation(double (*Gn)[3], double dt);

  void calcTurbTimeScale();

  void calcReStress(double *struc, double (*Gn)[3], double (*Gnph)[3], 
                    double (*Gnp1)[3], double dt);
};
