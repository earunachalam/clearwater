include ~/makefile.$(OS)

CXXFLAGS=-Wall -pedantic -O0 -g -ggdb -std=c++11 -march=native -fopenmp $(PLATFORMSPECIFICCXXFLAGS)
INCLUDES=$(PLATFORM_SPECIFIC_INCLUDES)
LIBS=$(PLATFORM_SPECIFIC_LIBS)
LINKFLAGS=
EXEC=../bin/read-lammpstrj

ALLOBJS=$(patsubst %.cpp,%.o,$(wildcard *.cpp))

$(EXEC): $(ALLOBJS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(ALLOBJS) $(LIBS) -o $(EXEC) $(LINKFLAGS) 

%.o:%.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $<

clean:
	rm -f *.o *.gch $(EXEC)

.phony: all clean
