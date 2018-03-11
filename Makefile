# location of the Boost Python include files and library

BOOST_LIB_PATH = $(BOOST_ROOT)/stage/lib/
BOOST_LIB = boost_python3 boost_iostreams boost_log boost_filesystem boost_system boost_program_options
 
# target files
TARGET = $(subst .cpp,,$(subst src/,,$(wildcard src/*/*.cpp)))
CTARGET = $(subst .c,,$(subst src/,,$(wildcard src/*/*.c)))
TARGET_OBJ = $(addprefix obj/,$(addsuffix .o,$(TARGET)))
CTARGET_OBJ = $(addprefix obj/,$(addsuffix .co,$(CTARGET)))

#flags
CC=gcc
CCFLAGS=-Wall -DBOOST_ALL_DYN_LINK -Werror -g -fPIC -std=c++11 -mavx2
CFLAGS=-Wall -DBOOST_ALL_DYN_LINK -Werror -g -fPIC -mavx2
LDSFLAGS=-shared -Wl,--export-dynamic
LDFLAGS=-g -std=c++11
LDLIBS=$(PYTHON_LIB) -L$(BOOST_LIB_PATH) $(addprefix -l,$(addsuffix $(BOOST_SUFFIX),$(BOOST_LIB))) -lm -lpthread -lstdc++ 

all: ma

ma: libMA.so src/cmdMa.cpp
	$(CC) $(CCFLAGS) src/cmdMa.cpp -isystem$(PYTHON_INCLUDE)/ -isystem$(BOOST_ROOT)/ -Iinc $(LDLIBS) libMA.so -o ma

libMA.so: $(TARGET_OBJ) $(CTARGET_OBJ)
	$(CC) $(LDFLAGS) $(LDSFLAGS) $(TARGET_OBJ) $(CTARGET_OBJ) $(LDLIBS) -o $@

obj/%.o: src/%.cpp inc/%.h
	$(CC) $(CCFLAGS) -isystem$(PYTHON_INCLUDE)/ -isystem$(BOOST_ROOT)/ -Iinc -c $< -o $@

obj/%.co: src/%.c inc/%.h
	$(CC) $(CFLAGS) -isystem$(PYTHON_INCLUDE)/ -isystem$(BOOST_ROOT)/ -Iinc -c $< -o $@

html/index.html: $(wildcard inc/*) $(wildcard inc/*/*) $(wildcard src/*) $(wildcard src/*/*) $(wildcard MA/*.py) doxygen.config
	doxygen doxygen.config

install: all
	pip3 install . --upgrade --no-cache-dir
	cp libMA.so /usr/lib
	pip3 show MA

##@todo remove me
vid:
	gource -f --seconds-per-day 0.1

distrib:
	python setup.py sdist bdist_egg bdist_wheel

clean:
	rm -f -r $(wildcard obj/*.o) $(wildcard obj/*/*.o) $(wildcard obj/*.co) $(wildcard obj/*/*.co) libMA.so
	rm -r -f dist *.egg-info build
	rm -r -f html

docs: html/index.html

.Phony: all clean install distrib docs vid