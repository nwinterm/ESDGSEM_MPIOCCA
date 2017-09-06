ifndef OCCA_DIR
ERROR:
	@echo "Error, environment variable [OCCA_DIR] is not set"
endif

include ${OCCA_DIR}/scripts/makefile
sPath = src
oPath = obj
#compilerFlags += -Ddfloat=float -Ddfloat4=float4 -g -O0 -DdfloatString=\"float\" -Ddfloat4String=\"float4\"
compilerFlags += -Ddfloat=double -Ddfloat4=double4 -g -O0 -DdfloatString=\"double\" -Ddfloat4String=\"double4\"
links += -llapack -lblas
#links += ./liblapack.so libblas.so


#---[ COMPILATION ]-------------------------------
headers = $(wildcard $(iPath)/*.hpp) $(wildcard $(iPath)/*.tpp)
sources = $(wildcard $(sPath)/*.cpp)

objects = $(subst $(sPath)/,$(oPath)/,$(sources:.cpp=.o))



main: $(objects) $(headers) main.cpp
	$(compiler) $(compilerFlags) -o main $(flags) $(objects) main.cpp $(paths) $(links)



$(oPath)/%.o:$(sPath)/%.cpp $(wildcard $(subst $(sPath)/,$(iPath)/,$(<:.cpp=.hpp))) $(wildcard $(subst $(sPath)/,$(iPath)/,$(<:.cpp=.tpp)))
	$(compiler) $(compilerFlags) -o $@ $(flags) -c $(paths) $<

clean:
	rm -f $(oPath)/*;
	rm -f main
#=================================================
