ROOT=brownTown3d
default: gnu

debug:
	g++ -ggdb -Wall $(ROOT).C -o $(ROOT) -lm -lgsl -lgslcblas

gnu:
	g++ -O3 -Wall -Wno-strict-overflow -Wno-unused-result -march=native -mfpmath=sse $(ROOT).C -o $(ROOT) -lm -lgsl -lgslcblas

clean:
	rm -f $(ROOT) $(ROOT)
