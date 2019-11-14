CC	    = g++
LIBS	= -ligraph 
LIBDIR  = -L/usr/local/lib -L$(PWD)/include/igraph/lib
INC     = -I./include -I$(PWD)/include/alglib -I$(PWD)/include/igraph/include/igraph
FLAGS   = -std=c++11 -g -Wl,-rpath=$(PWD)/include/igraph/lib,$(LIBDIR)

TARGET = bin/threshold

SRCDIR = src
INCLUDEDIR = include
SRCEXT = cpp

SOURCES =  $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))
SOURCES += $(shell find $(INCLUDEDIR)/alglib -type f -name *.$(SRCEXT))

BUILDDIR = build
OBJECTS = $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))

$(info $(SOURCES))
$(info $(OBJECTS))


$(TARGET): $(OBJECTS)
	$(CC) -o $@ $^ $(FLAGS) $(LIBS)

$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT)
	@mkdir -p $(BUILDDIR)
	$(CC) -c -o $@ $< $(FLAGS) $(INC) $(LIBDIR) $(LIBS)

clean:
	$(RM) -r $(BUILDDIR) $(TARGET)
