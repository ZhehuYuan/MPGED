default:
	g++ -c src/ged.cpp -O3 -Wall -std=c++11 -I~/local/tbb/include -L~/local/tbb/lib64 -ltbb -g -o src/ged.o
	g++ -o ged src/ged.o -pthread -L~/local/tbb/lib64 -ltbb -ltbbmalloc
#Tune
	g++ -c src/gedTune.cpp -O3 -Wall -std=c++11 -I~/local/tbb/include -L~/local/tbb/lib64 -ltbb -g -o src/ged.o
	g++ -o gedTune src/ged.o -pthread -L~/local/tbb/lib64 -ltbb -ltbbmalloc
#TimeGo
	g++ -c src/gedTimeGo.cpp -O3 -Wall -std=c++11 -I~/local/tbb/include -L~/local/tbb/lib64 -ltbb -g -o src/ged.o
	g++ -o gedTimeGo src/ged.o -pthread -L~/local/tbb/lib64 -ltbb -ltbbmalloc
#Count
	g++ -c src/gedCount.cpp -O3 -Wall -std=c++11 -I~/local/tbb/include -L~/local/tbb/lib64 -ltbb -g -o src/ged.o
	g++ -o gedCount src/ged.o -pthread -L~/local/tbb/lib64 -ltbb -ltbbmalloc
#Ub
	g++ -c src/gedUb.cpp -O3 -Wall -std=c++11 -I~/local/tbb/include -L~/local/tbb/lib64 -ltbb -g -o src/ged.o
	g++ -o gedUb src/ged.o -pthread -L~/local/tbb/lib64 -ltbb -ltbbmalloc

clean:
	rm src/*.o
	rm ged*
	rm Results/*
	rm ResultsTune/*

