
# The libraries for pythia, fastjet, and root ( "LIB_TRI" )
LIB_PYTH=-I${PYTHIA8}/include -L${PYTHIA8}/lib -lpythia8
LIB_FASTJET=`${FASTJET3}/fastjet-config --cxxflags --libs`
LIB_ROOT=`root-config --cflags --glibs`
LIB_TRI=${LIB_PYTH} ${LIB_FASTJET} ${LIB_ROOT}

# compilation option
CC=g++
CFLAGS=-O3 -Wno-deprecated

main_: main.cc
	${CC} ${CFLAGS} -o main main.cc ${LIB_TRI}

# or, more tersly
	# ${CC} ${CFLAGS} -o $@ $^ ${LIB_TRI}
