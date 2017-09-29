# location of the Python header files
 
PYTHON_VERSION = 3.5
PYTHON_INCLUDE = C:\Python35-64\include
 
# location of the Boost Python include files and library
 
BOOST_INC = C:/boost/boost_1_65_1  C:/Python35-32/include
BOOST_LIB_PATH = C:/boost/boost_1_65_1/stage/lib
BOOST_LIB_PRE = boost_system boost_python3 boost_thread boost_log boost_log_setup boost_filesystem boost_program_options boost_regex boost_iostreams
BOOST_LIB = $(addsuffix -mgw48-mt-1_65_1,$(BOOST_LIB_PRE))
LIB = python35
LIB_PATH = C:\Python35-64 C:\cpp-tools\lib64
 
# target files
TARGET = $(subst .cpp,,$(subst src/,,$(wildcard src/*.cpp)))
CTARGET = $(subst .c,,$(subst src/,,$(wildcard src/*.c)))
TARGET_OBJ = $(addprefix obj/,$(addsuffix .o,$(TARGET)))
CTARGET_OBJ = $(addprefix obj/,$(addsuffix .co,$(CTARGET)))

#flags
CC=g++
CCFLAGS= -Wall -std=c++11 -DBOOST_ALL_DYN_LINK -Werror -fno-threadsafe-statics
CFLAGS= -Wall -DBOOST_ALL_DYN_LINK -Werror
LDFLAGS= -shared -std=c++11 $(addprefix -L,$(LIB_PATH)) $(addprefix -l,$(LIB)) $(addprefix -L,$(BOOST_LIB_PATH)) $(addprefix -l,$(BOOST_LIB)) -pthread

all: LAuS.so

LAuS.so: $(TARGET_OBJ) $(CTARGET_OBJ)
	$(CC) $(LDFLAGS) $(TARGET_OBJ) $(CTARGET_OBJ) -o $@
 
obj/%.o: src/%.cpp inc/%.h
	$(CC) $(CCFLAGS) -I$(PYTHON_INCLUDE) $(addprefix -isystem,$(BOOST_INC)) -Iinc -c $< -o $@

obj/%.co: src/%.c inc/%.h
	$(CC) $(CFLAGS) -I$(PYTHON_INCLUDE) $(addprefix -isystem,$(BOOST_INC)) -Iinc -c $< -o $@

LAuS.html: $(wildcard src/*.cpp) $(wildcard inc/*.h)
	python -m pydoc -w LAuS

html/index.html: $(wildcard src/*.cpp) $(wildcard inc/*.h)
	doxygen doxygen.config

install: LAuS.so
	pip install . --upgrade

distrib:
	python setup.py sdist bdist_egg bdist_wheel

clean:
	rm -f -r $(wildcard obj/*.o) $(wildcard obj/*.co) *.so
	rm -r -f dist *.egg-info build
	rm -r -f html
	rm -r -f LAuS.html

docs: html/index.html LAuS.html

.Phony: all clean install distrib docs