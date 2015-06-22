objects= Part.o rand.o Cont.o

molecularDynamic: $(objects)
	g++-4.6 -fopenmp -Wall $(objects) molecularDynamic.cpp -o $@
	
$(objects):%.o:%.cpp %.hpp
	g++-4.6 -fopenmp -Wall -c $< -o $@

clean:
	rm -f *.o
	
