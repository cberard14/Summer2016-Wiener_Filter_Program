include ./include/ElVars

all: WF1D

WF1D: WF1D.o WienerFilter1D.o MatrixSetUpFunctions.o VectorFunctions.o
	$(CXX) $(EL_LINK_FLAGS) -o WF1D WF1D.o WienerFilter1D.o MatrixSetUpFunctions.o VectorFunctions.o $(EL_LIBS)

WF1D.o: WF1D.cpp ./include/VectorFunctions.hpp ./include/MatrixSetUpFunctions.hpp ./include/WienerFilter1D.hpp
	$(CXX) $(EL_COMPILE_FLAGS) -c WF1D.cpp

WienerFilter1D.o: ./src/WienerFilter1D.cpp ./include/VectorFunctions.hpp ./include/MatrixSetUpFunctions.hpp
	$(CXX) $(EL_COMPILE_FLAGS) -c ./src/WienerFilter1D.cpp

VectorFunctions.o: ./src/VectorFunctions.cpp ./include/VectorFunctions.hpp
	$(CXX) $(EL_COMPILE_FLAGS) -c ./src/VectorFunctions.cpp

MatrixSetUpFunctions.o: ./src/MatrixSetUpFunctions.cpp ./include/VectorFunctions.hpp ./include/MatrixSetUpFunctions.hpp
	$(CXX) $(EL_COMPILE_FLAGS) -c ./src/MatrixSetUpFunctions.cpp

clean:
	rm *o WF1D