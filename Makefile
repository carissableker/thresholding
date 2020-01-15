CC		= g++
LIBS	= -ligraph 
LIBDIRS	= -L/usr/local/lib -L$(PWD)/lib
INC		= -I$(PWD)/include -I$(PWD)/include/alglib -I$(PWD)/include/igraph
FLAGS	= -std=c++11 -g -Wl,-rpath=$(PWD)/lib

TARGET	= $(PWD)/bin/threshold

SRCDIR	= $(PWD)/src
SRCEXT	= cpp

INCLUDEDIR	= $(PWD)/include
LIBDIR 		= $(PWD)/lib
BUILDDIR 	= $(PWD)/build
EXTERNALDIR = $(PWD)/external/
IGRAPHDIR   = igraph-0.7.1

SOURCES  =  $(shell find $(SRCDIR) -xtype f -name "*".$(SRCEXT))
OBJECTS = $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))

$(info $(SOURCES))
$(info $(OBJECTS))


$(TARGET): $(OBJECTS) 
	$(CC) -o $@ $^ $(FLAGS) $(LIBDIRS) $(LIBS) 

$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT)  #igraph
	@mkdir -p $(BUILDDIR)
	$(CC) -c -o $@ $< $(FLAGS) $(INC) $(LIBDIRS) $(LIBS) 

.PHONY: igraph 
igraph:
	cd $(EXTERNALDIR) && \
	tar zxvf $(IGRAPHDIR).tar.gz && \
	cd $(IGRAPHDIR) && \
	./configure --libdir=$(LIBDIR) --prefix=$(BUILDDIR) && \
	$(MAKE) && \
	$(MAKE) install #&& \
	#cp $(BUILDDIR)/include/igraph/* $(INCLUDEDIR)/igraph/

.PHONY: clean
clean:
	$(RM) -r $(BUILDDIR) $(TARGET)
	#$(RM) $(LIBDIR)/lib*
	#$(MAKE) --directory=$(IGRAPHDIR) clean


# http://nuclear.mutantstargoat.com/articles/make/#building-sub-projects 

# Doesn't work, see https://github.com/igraph/igraph/issues/1263#issue-523750694
# ./configure --includedir=$(INCLUDEDIR)/igraph --libdir=$(LIBDIR) --prefix=$(BUILDDIR) 


