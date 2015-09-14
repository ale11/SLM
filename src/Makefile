main: UgpWithCvCompFlow.o TurbModel_WVM.o main.o
	$(CPP) -o $(BINDIR)/main UgpWithCvCompFlow.o TurbModel_WVM.o main.o $(LIBS)

UgpWithCvCompFlow.o: UgpWithCvCompFlow.cpp UgpWithCvCompFlow.hpp mathoper.hpp
	$(CPP) -c UgpWithCvCompFlow.cpp -I. $(INCS) $(FLAGS)

TurbModel_WVM.o: TurbModel_WVM.cpp TurbModel_WVM.hpp UgpWithCvCompFlow.hpp mathoper.hpp
	$(CPP) -c TurbModel_WVM.cpp -I. $(INCS) $(FLAGS)

main.o: main.cpp TurbModel_WVM.hpp UgpWithCvCompFlow.hpp mathoper.hpp
	$(CPP) -c main.cpp -I. $(INCS) $(FLAGS)

prm_vars: prm_vars.o
	$(CPP) -o $(BINDIR)/prm_vars prm_vars.o

prm_vars.o: prm_vars.cpp
	$(CPP) -c prm_vars.cpp