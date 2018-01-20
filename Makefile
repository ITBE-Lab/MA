# location of the Boost Python include files and library

BOOST_LIB_PATH = $(BOOST_ROOT)/stage/lib/
BOOST_LIB = boost_python3 boost_iostreams boost_log boost_filesystem boost_system
 
# target files
TARGET = $(subst .cpp,,$(subst src/,,$(wildcard src/*.cpp))) $(subst .cpp,,$(subst src/,,$(wildcard src/*/*.cpp)))
CTARGET = $(subst .c,,$(subst src/,,$(wildcard src/*.c))) $(subst .c,,$(subst src/,,$(wildcard src/*/*.c)))
TARGET_OBJ = $(addprefix obj/,$(addsuffix .o,$(TARGET)))
CTARGET_OBJ = $(addprefix obj/,$(addsuffix .co,$(CTARGET)))

#flags
CC=gcc
CCFLAGS=-Wall -DBOOST_ALL_DYN_LINK -Werror -g -std=c++11 -pthread
CFLAGS=-Wall -DBOOST_ALL_DYN_LINK -Werror -g -pthread
LDSFLAGS=-shared
LDFLAGS=-g -std=c++11
LDLIBS=$(PYTHON_LIB) -L$(BOOST_LIB_PATH) $(addprefix -l,$(addsuffix $(BOOST_SUFFIX),$(BOOST_LIB))) -lm -lpthread -lstdc++ 

ALL=ml.exe
# add flags required for wondows or linux specifically
ifeq ($(OS),Windows_NT)
	LDSFLAGS+= -Wl,--export-all-symbols
	ALL+=
else
	CCFLAGS+= -fPIC -mavx2 
	CFLAGS+= -fPIC
	LDSFLAGS+= -Wl,--export-dynamic
	ALL+= libMABS.so
endif

all: $(ALL)

ml.exe: $(TARGET_OBJ) $(CTARGET_OBJ)
	$(CC) $(LDFLAGS) $(TARGET_OBJ) $(CTARGET_OBJ) -o $@ $(LDLIBS)

libMABS.so: $(TARGET_OBJ) $(CTARGET_OBJ)
	$(CC) $(LDFLAGS) $(LDSFLAGS) $(TARGET_OBJ) $(CTARGET_OBJ) $(LDLIBS) -o $@

obj/%.o: src/%.cpp inc/%.h
	$(CC) $(CCFLAGS) -isystem$(PYTHON_INCLUDE)/ -isystem$(BOOST_ROOT)/ -Iinc -c $< -o $@

obj/%.co: src/%.c inc/%.h
	$(CC) $(CFLAGS) -isystem$(PYTHON_INCLUDE)/ -isystem$(BOOST_ROOT)/ -Iinc -c $< -o $@

html/index.html: $(wildcard inc/*.h) $(wildcard inc/*/*.h) $(wildcard MABS/*.py) doxygen.config
	doxygen doxygen.config

install: all
	pip3 install . --upgrade --no-cache-dir
	cp libMABS.so /usr/lib
	pip3 show MABS

#@todo remove me
vid:
	gource -f --seconds-per-day 0.1

distrib:
	python setup.py sdist bdist_egg bdist_wheel

clean:
	rm -f -r $(wildcard obj/*.o) $(wildcard obj/*/*.o) $(wildcard obj/*.co) $(wildcard obj/*/*.co) libMABS.so
	rm -r -f dist *.egg-info build
	rm -r -f html

docs: html/index.html

.Phony: all clean install distrib docs vid