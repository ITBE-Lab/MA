# location of the Python header files
 
PYTHON_VERSION = 2.7
PYTHON_INCLUDE = /usr/include/python$(PYTHON_VERSION)
 
# location of the Boost Python include files and library
 
BOOST_INC = /usr/include
BOOST_LIB = /usr/lib
 
# target files
TARGET = $(subst .cpp,,$(subst src/,,$(wildcard src/*.cpp)))
TARGET_OBJ = $(addprefix obj/,$(addsuffix .o,$(TARGET)))

#flags
CC=gcc
CCFLAGS= -Wall -fPIC -std=c++11
LDFLAGS= -shared -Wl,--export-dynamic -std=c++11 -L$(BOOST_LIB) -lboost_python -L/usr/lib/python$(PYTHON_VERSION)/config -lpython$(PYTHON_VERSION)
LDFLAGS2= -L$(BOOST_LIB) -lboost_python -std=c++11 -L/usr/lib/python$(PYTHON_VERSION)/config -lpython$(PYTHON_VERSION)

LAuS_Aligner.so: $(TARGET_OBJ)
	$(CC) $(LDFLAGS) $(TARGET_OBJ) -o $@
 
obj/%.o: src/%.cpp inc/%.h
	$(CC) $(CCFLAGS) -I$(PYTHON_INCLUDE) -I$(BOOST_INC) -Iinc -c $< -o $@

all: LAuS_Aligner.so

install:
	pip install . --upgrade

distrib:
	python setup.py sdist
	python setup.py bdist_wheel

clean:
	rm -f -r $(wildcard obj/*.o) *.so
	rm -r -f dist *.egg-info build

.Phony: all clean install distrib