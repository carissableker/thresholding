CC=g++
FLAGS= -std=c++11 -g -I /storage/langston/biodata/carissa/thresholding_project/thresholding/my_lib/igraphc/include/igraph/ 
LIBS=-L /storage/langston/biodata/carissa/thresholding_project/thresholding/my_lib/igraphc/lib/ -ligraph 
EXECUTABLE=threshold

threshold: thresholding.cpp
	$(CC) $(FLAGS) $(LIBS) thresholding.cpp -o $(EXECUTABLE)

clean:
	rm -rf core* *.o $(EXECUTALBE)


