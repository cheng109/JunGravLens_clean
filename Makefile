WORK_DIR := $(shell pwd)
CXX := g++
FC 	:= gfortran
INC := -I$(WORK_DIR)/lib/cfitsio -I$(WORK_DIR)/lib/eigen_include
CXXFLAGS := -Wall -O3 -g -std=c++11 $(INC) -fopenmp
LDFLAGS := -L. -L$(WORK_DIR)/lib/cfitsio -L/usr/local/gfortran/lib -fopenmp
LDLIBS := -lcfitsio -lfortranstuff -lgfortran 
#-lgsl -lgslcblas -lm -lfortranstuff

SRCS := main.cpp Image.cpp commons.cpp Model.cpp parafit.cpp Kmeans.cpp parafit2.cpp mc.cpp
OBJS := $(patsubst %.cpp, %.o, $(SRCS))

## ARG1:   working directory. 
## ARG2:   configuration file. 
## ARG3:   output files. 

ARG1 :=sie_test/
ARG2 :=conf.txt
ARG3 :=output_new.txt

appname := junGL 

all: $(appname)
	./$(appname) $(ARG1) $(ARG2) $(ARG3)
	@#valgrind --tool=memcheck --leak-check=full --verbose --log-file=memcheck.log --show-leak-kinds=all --track-origins=yes ./$(appname) $(ARG1) $(ARG2) $(ARG3)

$(appname): libfortranstuff.a $(OBJS) 
	$(CXX) $(LDFLAGS) -o $(appname) $(OBJS) $(LDLIBS)

#depend: .depend

.depend: $(SRCS)
	rm -f ./.depend
	$(CXX) $(CXXFLAGS) -MM $^>>./.depend;

clean:
	rm -f $(OBJS) $(appname) output.txt libfortranstuff.a *.o

dist-clean: clean
	rm -f *~ .depend

libfortranstuff.a:
	@$(FC) -O -c slatec/src/*.f fastell.f 
	@ar -r libfortranstuff.a *.o
	rm *.o

plot:
	@./plotScript

kmeans: 
	@$(CXX) $(CXXFLAGS) Kmeans.cpp && ./a.out
	@rm a.out

include .depend



