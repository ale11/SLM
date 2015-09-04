main: UgpWithCvCompFlow.o TurbModel_WVM.o main.o
	$(CPP) -o $(BINDIR)/main UgpWithCvCompFlow.o TurbModel_WVM.o main.o -L/usr/local/lib -larmadillo -framework Accelerate

UgpWithCvCompFlow.o: UgpWithCvCompFlow.cpp UgpWithCvCompFlow.hpp mathoper.hpp
	$(CPP) -c UgpWithCvCompFlow.cpp -I. -I/usr/local/include -O2

TurbModel_WVM.o: TurbModel_WVM.cpp TurbModel_WVM.hpp UgpWithCvCompFlow.hpp mathoper.hpp
	$(CPP) -c TurbModel_WVM.cpp -I. -I/usr/local/include -O2

main.o: main.cpp TurbModel_WVM.hpp UgpWithCvCompFlow.hpp mathoper.hpp
	$(CPP) -c main.cpp -I. -I/usr/local/include -O2

prm_vars: prm_vars.o
	$(CPP) -o $(BINDIR)/prm_vars prm_vars.o

prm_vars.o: prm_vars.cpp
	$(CPP) -c prm_vars.cpp