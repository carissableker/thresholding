CC	   	 = g++
LIBS	 = -ligraph 
LIBDIRS  = -L/usr/local/lib -L$(PWD)/lib
INC      = -I./include -I$(PWD)/include/alglib -I$(PWD)/include/igraph
FLAGS    = -std=c++11 -g -Wl,-rpath=$(PWD)/lib,$(LIBDIRS)

TARGET   = bin/threshold

SRCDIR   = src
SRCEXT   = cpp

INCLUDEDIR = $(PWD)/include
LIBDIR     = $(PWD)/lib

IGRAPHDIR  = $(PWD)/external/igraph-0.7.1

SOURCES  =  $(shell find $(SRCDIR) -type f -name "*".$(SRCEXT))
SOURCES += $(shell find $(INCLUDEDIR)/alglib -type f -name "*".$(SRCEXT))

BUILDDIR = $(PWD)/build
OBJECTS = $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))

$(info $(SOURCES))
$(info $(OBJECTS))

$(TARGET): $(OBJECTS)
	$(CC) -o $@ $^ $(FLAGS) $(LIBS)

$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT) igraph
	@mkdir -p $(BUILDDIR)
	$(CC) -c -o $@ $< $(FLAGS) $(INC) $(LIBDIRS) $(LIBS)

.PHONY: igraph 
igraph:
	cd $(IGRAPHDIR) && \
	./configure --libdir=$(LIBDIR) --prefix=$(BUILDDIR) && \
	$(MAKE) && \
	$(MAKE) install && \
	cp $(BUILDDIR)/include/igraph/* $(INCLUDEDIR)/igraph/

.PHONY: clean
clean:
	$(RM) -r $(BUILDDIR) $(TARGET)
	#$(RM) $(LIBDIR)/lib*
	$(MAKE) --directory=$(IGRAPHDIR) clean


# http://nuclear.mutantstargoat.com/articles/make/#building-sub-projects 
# Doesn't work, see https://github.com/igraph/igraph/issues/1263#issue-523750694
# ./configure --includedir=$(INCLUDEDIR)/igraph --libdir=$(LIBDIR) --prefix=$(BUILDDIR) 


