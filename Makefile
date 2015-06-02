all: EffAct testT

testT.o: testT.c intT.h Makefile
	/usr/local/cuda/bin/nvcc -c -g testT.c -L/home/mazur/software/lib -lgsl -lgslcblas -lm
Tint.o:  Tint.c intT.h
	/usr/local/cuda/bin/nvcc -c -g Tint.c -L/home/mazur/software/lib -lgsl -lgslcblas -lm

Integrand.o: Integrand.c intT.h Makefile
	/usr/local/cuda/bin/nvcc -c -g Integrand.c -L/home/mazur/software/lib -lgsl -lgslcblas -lm

getwl.o: getwl.c Makefile
	/usr/local/cuda/bin/nvcc -c getwl.c

intEV.o: intEV.cu intT.h Makefile
	/usr/local/cuda/bin/nvcc -c -g -arch sm_13 -prec-div=false --ptxas-options="-v" -prec-sqrt=false intEV.cu -lm
intEVff.o: intEVff.cu intT.h Makefile
	/usr/local/cuda/bin/nvcc -c -g -arch sm_13 -prec-div=false --ptxas-options="-v" -prec-sqrt=false intEVff.cu -lm

EffAct.o: EffAct.c intT.h Makefile
	/usr/local/cuda/bin/nvcc -c -g EffAct.c

periodic.o: periodic.c intT.h Makefile
	/usr/local/cuda/bin/nvcc -c -g periodic.c -L/home/mazur/software/lib -lgsl -lgslcblas -lm



EffAct: EffAct.o intEV.o getwl.o Tint.o Integrand.o periodic.o Makefile
	/usr/local/cuda/bin/nvcc -arch sm_13 -prec-div=false -prec-sqrt=false --ptxas-options="-v" -g -o EffAct EffAct.o intEV.o getwl.o Tint.o Integrand.o periodic.o -L/home/mazur/software/lib -lcuda -lcudart -lgsl -lgslcblas -lm

testT:  testT.o getwl.o Tint.o Integrand.o intEV.o
	/usr/local/cuda/bin/nvcc -arch sm_13 -g -o testT testT.o intEV.o getwl.o Tint.o Integrand.o -L/home/mazur/software/lib -lcuda -lcudart -lgsl -lgslcblas -lm
clean:
	rm *.o EffAct
