CC                   = gcc # icc
CXX 		      	 = g++ # icpc
OPTIMIZE_FLAGS   	 = -O3 -finline-functions # -fargument-noalias-global -restrict -qopt-report=1 -qopt-report-phase=vec -fno-alias
CXXFLAGS             = -DNO_DEFAULT_MAPPING -std=c++0x -Wall
INCS                 = -I $(FF_HOME)
LDFLAGS 			 = -pthread

SRC_HOME			 = src
TARGET_FAST		     = jacobi_seq_mmic jacobi_par_mmic gauss_seq_mmic gauss_par_mmic gaussrb_seq_mmic jacobi_seq jacobi_par gauss_seq gaussrb_seq gauss_par	

SOURCES				 = $(SRC_HOME)/PoissonMain.cpp $(SRC_HOME)/_2DDirichlet.cpp $(SRC_HOME)/_2DGrid.cpp $(SRC_HOME)/_2DPoissonEquation.cpp

all: $(TARGET)

jacobi_seq_mmic: $(SOURCES) $(SRC_HOME)/_2DSequentialJacobi.cpp
	$(CXX) $(CXXFLAGS) $(INCS) -mmic -DSEQ_JAC -o $@ $^ $(LDFLAGS)

jacobi_par_mmic: $(SOURCES) $(SRC_HOME)/_2DParallelJacobi.cpp
	$(CXX) $(INCS) $(CXXFLAGS) $(OPTIMIZE_FLAGS) -mmic -DPAR_JAC -o $@ $^ $(LDFLAGS)

gauss_seq_mmic: $(SOURCES) $(SRC_HOME)/_2DSequentialGaussSeidel.cpp
	$(CXX) $(INCS) -std=c++11 -mmic -DSEQ_GAUSS -o $@ $^ $(LDFLAGS)

gaussrb_seq_mmic: $(SOURCES) $(SRC_HOME)/_2DSequentialGaussSeidelRedBlack.cpp
	$(CXX) $(INCS) -std=c++11 -mmic -DSEQ_RBGAUSS -o $@ $^ $(LDFLAGS)

gauss_par_mmic: $(SOURCES) $(SRC_HOME)/_2DParallelGaussSeidelRedBlack.cpp
	$(CXX) $(INCS) -std=c++11 -mmic -DPAR_GAUSS $(CXXFLAGS) $(OPTIMIZE_FLAGS) -o $@ $^ $(LDFLAGS)

jacobi_seq: $(SOURCES) $(SRC_HOME)/_2DSequentialJacobi.cpp
	$(CXX) $(CXXFLAGS) $(INCS) -DSEQ_JAC -o $@ $^ $(LDFLAGS)

gaussrb_seq: $(SOURCES) $(SRC_HOME)/_2DSequentialGaussSeidelRedBlack.cpp
	$(CXX) $(INCS) -std=c++11 -DSEQ_RBGAUSS -o $@ $^ $(LDFLAGS)

gauss_seq: $(SOURCES) $(SRC_HOME)/_2DSequentialGaussSeidel.cpp
	$(CXX) $(INCS) -std=c++11 -DSEQ_GAUSS -o $@ $^ $(LDFLAGS)	

jacobi_par: $(SOURCES) $(SRC_HOME)/_2DParallelJacobi.cpp
	$(CXX) $(INCS) $(CXXFLAGS) $(OPTIMIZE_FLAGS) -DPAR_JAC -o $@ $^ $(LDFLAGS)

gauss_par: $(SOURCES) $(SRC_HOME)/_2DParallelGaussSeidelRedBlack.cpp
	$(CXX) $(INCS) -std=c++11 -DPAR_GAUSS $(CXXFLAGS) $(OPTIMIZE_FLAGS) -o $@ $^ $(LDFLAGS)	



clean: 
	-rm -fr *.o *~
cleanall: clean
	-rm -fr $(TARGET) *.d 


