CXX=mpicxx -std=c++20
PACS_ROOT?=../pacs-examples/Examples
CXXFLAG= -fopenmp -DMPI -I$(PACS_ROOT)/include

# Usually could be disable since I'm trouble with installation of muParser
MUPARSER_LIBDIR=-L$(PACS_ROOT)/lib
LIBRARY_NAME=muparser
DYNAMIC_LIBFILE=lib$(LIBRARY_NAME).so
LDLIBS+=-L. -l$(LIBRARY_NAME) $(MUPARSER_LIBDIR) -lmuparser -L$(PACS_ROOT)/lib -lpacs
CXXFLAGS+=$(MUPARSER_INCLUDE)

EXEC=main
SRC=main.cpp problem.cpp
OBJECTS=$(SRC:.cpp=.o)
HEADERS=problem.hpp

$(EXEC) : $(OBJECTS)
	$(CXX) $(CXXFLAG) $^ -o $@ $(LDLIBS)

%.o: %.cpp $(HEADERS)
	$(CXX) $(CXXFLAG) -c $< -o $@

$()

clean:
	rm -f *.o $(EXEC)

optimize: CXXFLAGS += -O3
optimize: $(EXEC)
