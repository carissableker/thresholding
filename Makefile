CC=g++
FLAGS= -std=c++11 -g -I /home/cbleker/my_lib/igraphc/include/igraph/ 
LIBS=-L /home/cbleker/my_lib/igraphc/lib/ -ligraph 
EXECUTABLE=threshold

threshold: thresholding.cpp
	$(CC) $(FLAGS) $(LIBS) thresholding.cpp -o $(EXECUTABLE)

clean:
	rm -rf core* *.o $(EXECUTALBE)


