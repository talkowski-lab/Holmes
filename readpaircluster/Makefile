#CFLAGS=-I /scr/talkowski/tools/bin/samtools-current -L /scr/talkowski/tools/bin/samtools-current -lbam -L ~/local/lib -lz
#CFLAGS=-I /data/talkowski/tools/bin/samtools-0.1.18 -L /data/talkowski/tools/bin/samtools-0.1.18 -lbam -lz
#CFLAGS=-I /apps/source/samtools-0.1.18 -L /apps/source/samtools-0.1.18 -I /apps/source/zlib/1.2.8/include -lbam -L /apps/source/zlib/1.2.8/lib -lz 

readpairCluster: cluster.o readPair.o
	g++ -O3 -o readPairCluster cluster.o readPair.o $(CFLAGS)

cluster.o: cluster.cpp simpleGraph.hpp
	g++ -O3 -c cluster.cpp $(CFLAGS)

readPair.o: readPair.cpp readPair.hpp 
	g++ -O3 -c readPair.cpp $(CFLAGS)
