# location of the Python header files
 
PYTHON_VERSION = 2.7
PYTHON_INCLUDE = /usr/include/python$(PYTHON_VERSION)
 
# location of the Boost Python include files and library
 
BOOST_INC = /usr/include
BOOST_LIB = /usr/lib
 
# target files
TARGET = core/test core/aligner modules/module data/container data/nucSeq

#flags
CC=gcc
CCFLAGS= -Wall -fPIC -std=c++11
LDFLAGS= -shared -Wl,--export-dynamic -std=c++11 -L$(BOOST_LIB) -lboost_python -L/usr/lib/python$(PYTHON_VERSION)/config -lpython$(PYTHON_VERSION)
LDFLAGS2= -L$(BOOST_LIB) -lboost_python -std=c++11 -L/usr/lib/python$(PYTHON_VERSION)/config -lpython$(PYTHON_VERSION)

aligner/%.so: obj/%.o
	$(CC) $(LDFLAGS) $< -o $@
 
obj/%.o: src/%.cpp inc/%.h
	$(CC) $(CCFLAGS) -I$(PYTHON_INCLUDE) -I$(BOOST_INC) -Iinc -c $< -o $@

all: $(addprefix aligner/,$(addsuffix .so,$(TARGET))) $(addprefix obj/,$(addsuffix .o,$(TARGET)))

clean:
	rm -f -r $(addprefix aligner/,$(addsuffix .so,$(TARGET))) $(addprefix obj/,$(addsuffix .o,$(TARGET)))

.Phony: all clean