main: main.o TurbModel_WVM.o UgpWithCvCompFlow.o
	$(CPP) -o $(BINDIR)/main main.o TurbModel_WVM.o UgpWithCvCompFlow.o

main.o: main.cpp TurbModel_WVM.hpp UgpWithCvCompFlow.hpp
	$(CPP) -c main.cpp -I.

TurbModel_WVM.o: TurbModel_WVM.cpp TurbModel_WVM.hpp UgpWithCvCompFlow.hpp
	$(CPP) -c TurbModel_WVM.cpp -I.

UgpWithCvCompFlow.o: UgpWithCvCompFlow.cpp UgpWithCvCompFlow.hpp mathoper.hpp
	$(CPP) -c UgpWithCvCompFlow.cpp -I.

prm_vars: prm_vars.o
	$(CPP) -o $(BINDIR)/prm_vars prm_vars.o

prm_vars.o: prm_vars.cpp
	$(CPP) -c prm_vars.cpp