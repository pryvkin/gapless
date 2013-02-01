PROG = gapless

CXX = g++

# libraries

# compiler flags
CFLAGS = $(LIB_CFLAGS) -O3 -Wall
LDFLAGS = $(LIB)

PROG_DEPENDS = $(OBJECTS)

SRCS = read_fasta.cpp gapless.cpp

OBJECTS = $(SRCS:.cpp=.o)

all: $(PROG)

clean:
	rm -f *.o $(PROG)

$(PROG): $(PROG_DEPENDS)
	${CXX} ${LDFLAGS} $(CFLAGS) $(OBJECTS) -o $(PROG)

.cpp.o:
	${CXX} ${CFLAGS} -c $<


