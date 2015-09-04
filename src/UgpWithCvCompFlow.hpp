#ifndef UGPWITHCVCOMPFLOW_H
#define UGPWITHCVCOMPFLOW_H

#include <iostream>
#include <cstdlib>
#include <string.h>
#include <ctime>
#include <iomanip>
#include "mathoper.hpp"

using namespace std;
using namespace arma;

class UgpWithCvCompFlow
{
  // Constructors & Destructor
public:
  UgpWithCvCompFlow();
  virtual ~UgpWithCvCompFlow();
  
  // Data members
public:
  double rey[6], dim[6], cir[6];
  double tke;
};

#endif
