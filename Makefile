# location of the Python header files
 
PYTHON_VERSION = 3.5
PYTHON_INCLUDE = /usr/include/python$(PYTHON_VERSION)
 
# location of the Boost Python include files and library
 
BOOST_INC = $(BOOST_ROOT)
BOOST_LIB_PATH = $(BOOST_ROOT)/stage/lib
BOOST_LIB = boost_python3-mt dl rt z boost_system-mt boost_thread-mt boost_log-mt boost_log_setup-mt boost_filesystem-mt boost_program_options-mt boost_regex-mt boost_iostreams-mt
 
# target files
TARGET = $(subst .cpp,,$(subst src/,,$(wildcard src/*.cpp))) $(subst .cpp,,$(subst src/,,$(wildcard src/*/*.cpp)))
CTARGET = $(subst .c,,$(subst src/,,$(wildcard src/*.c))) $(subst .c,,$(subst src/,,$(wildcard src/*/*.c)))
TARGET_OBJ = $(addprefix obj/,$(addsuffix .o,$(TARGET)))
CTARGET_OBJ = $(addprefix obj/,$(addsuffix .co,$(CTARGET)))

#flags
CC=gcc
CCFLAGS= -Wall -fPIC -std=c++11 -DBOOST_ALL_DYN_LINK -Werror -g -mavx2
CFLAGS= -Wall -fPIC -DBOOST_ALL_DYN_LINK -Werror -g
LDFLAGS= -shared -Wl,--export-dynamic -std=c++11 -L$(BOOST_LIB_PATH) $(addprefix -l,$(BOOST_LIB)) -L/usr/lib/python$(PYTHON_VERSION)/config-x86_64-linux-gnu -lpython$(PYTHON_VERSION) -pthread -g

all: libMABS.so

libMABS.so: $(TARGET_OBJ) $(CTARGET_OBJ)
	$(CC) $(LDFLAGS) $(TARGET_OBJ) $(CTARGET_OBJ) -o $@

obj/%.o: src/%.cpp inc/%.h
	$(CC) $(CCFLAGS) -I$(PYTHON_INCLUDE) $(addprefix -isystem,$(BOOST_INC)) -Iinc -c $< -o $@

obj/%.co: src/%.c inc/%.h
	$(CC) $(CFLAGS) -I$(PYTHON_INCLUDE) $(addprefix -isystem,$(BOOST_INC)) -Iinc -c $< -o $@

html/index.html: $(wildcard inc/*.h) $(wildcard inc/*/*.h) $(wildcard MABS/*.py) doxygen.config
	doxygen doxygen.config

install: all
	pip3 install . --upgrade --no-cache-dir
	cp libMABS.so /usr/lib
	pip3 show MABS

#@todo  remove me
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