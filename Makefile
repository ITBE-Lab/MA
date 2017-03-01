# location of the Python header files
 
PYTHON_VERSION = 2.7.9
PYTHON_INCLUDE = /usr/include/python$(PYTHON_VERSION)
 
# location of the Boost Python include files and library
 
BOOST_INC = /usr/include
BOOST_LIB = /usr/lib
 
# target files
TARGET = core/test core/aligner modules/module data/container data/nucSeq

#flags
CC=gcc
CCFLAGS= -Wall -fPIC
LDFLAGS= -shared -Wl --export-dynamic -L$(BOOST_LIB) -lboost_python-$(PYTHON_VERSION) -L/usr/lib/python$(PYTHON_VERSION)/config -lpython$(PYTHON_VERSION)
INCLUDES= $(wildcard ./inc/*) $(PYTHON_INCLUDE) $(BOOST_INC)
OBJECTS= $(SOURCES:.cpp=.o) 

 
$(TARGET).so: $(TARGET).o
	$(CC) $(LDFLAGS) obj/$(TARGET).o -o aligner/_$(TARGET).so
 
$(TARGET).o: $(TARGET).cpp
	$(CC) -I$(INCLUDES) -c src/$(TARGET).cpp -o obj/$(TARGET).o

all: $(TARGET)

clean:
	rm -f -r *.so *.o 

.phony all clean