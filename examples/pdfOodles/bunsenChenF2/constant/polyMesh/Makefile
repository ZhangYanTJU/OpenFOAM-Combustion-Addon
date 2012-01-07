
blockMeshDict: a.out
	./a.out

a.out: generate.cpp $(FOAM_UTILITIES)/LTT/blockMeshToolkit/infra.H
	g++ generate.cpp -I$(FOAM_UTILITIES)/LTT/blockMeshToolkit -I$(OCTAVE_ARCH_PATH)/include/octave -L$(OCTAVE_ARCH_PATH)/lib -loctave -lcruft

