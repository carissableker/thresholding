CC	    = g++
FLAGS   = -std=c++11 -g  
LIBS	= -ligraph 
LIBDIR  = -L/usr/local/lib 
INC     = -I/usr/local/include/igraph -I./include

TARGET = bin/threshold

SRCDIR = src
SRCEXT = cpp
SOURCES = $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))

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
