CC	= g++
FLAGS   = -std=c++11 -g  
LIBS	= -ligraph 
LIB_DIR = -L/usr/local/lib
INC     = -I/usr/local/include/igraph

EXECUTABLE = threshold

threshold: ./src/thresholding.cpp
	$(CC) $(FLAGS) ./src/thresholding.cpp $(INC) $(LIB_DIR) $(LIBS) -o $(EXECUTABLE)

clean:
	rm -rf core* *.o $(EXECUTALBE)


