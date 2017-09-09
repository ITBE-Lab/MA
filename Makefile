# location of the Python header files
 
PYTHON_VERSION = 2.7
PYTHON_INCLUDE = /usr/include/python$(PYTHON_VERSION)
 
# location of the Boost Python include files and library
 
BOOST_INC = /opt/dev/boost_1_60_0__
BOOST_LIB_PATH = /opt/dev/boost_1_60_0__/stage/lib /opt/dev/boost_1_60_0__/lib
BOOST_LIB = boost_python dl rt z boost_system-mt boost_thread-mt boost_log-mt boost_log_setup-mt boost_filesystem-mt boost_program_options-mt boost_regex-mt boost_iostreams
 
# target files
TARGET = $(subst .cpp,,$(subst src/,,$(wildcard src/*.cpp)))
CTARGET = $(subst .c,,$(subst src/,,$(wildcard src/*.c)))
TARGET_OBJ = $(addprefix obj/,$(addsuffix .o,$(TARGET)))
CTARGET_OBJ = $(addprefix obj/,$(addsuffix .co,$(CTARGET)))

#flags
CC=gcc
CCFLAGS= -Wall -fPIC -std=c++11 -DBOOST_ALL_DYN_LINK -Werror -fprofile-arcs -ftest-coverage
CFLAGS= -Wall -fPIC -DBOOST_ALL_DYN_LINK -Werror -fprofile-arcs -ftest-coverage
LDFLAGS= -shared -Wl,--export-dynamic -std=c++11 $(addprefix -L,$(BOOST_LIB_PATH)) $(addprefix -l,$(BOOST_LIB)) -L/usr/lib/python$(PYTHON_VERSION)/config-x86_64-linux-gnu -lpython$(PYTHON_VERSION) -pthread -fprofile-arcs

all: LAuS.so LAuS.html

LAuS.so: $(TARGET_OBJ) $(CTARGET_OBJ)
	$(CC) $(LDFLAGS) $(TARGET_OBJ) $(CTARGET_OBJ) -o $@
 
obj/%.o: src/%.cpp inc/%.h
	$(CC) $(CCFLAGS) -I$(PYTHON_INCLUDE) $(addprefix -isystem,$(BOOST_INC)) -Iinc -c $< -o $@

obj/%.co: src/%.c inc/%.h
	$(CC) $(CFLAGS) -I$(PYTHON_INCLUDE) $(addprefix -isystem,$(BOOST_INC)) -Iinc -c $< -o $@

LAuS.html: $(wildcard src/*.cpp) $(wildcard inc/*.h)
	pydoc -w LAuS


install: LAuS.so
	pip install . --upgrade

distrib:
	python setup.py sdist bdist_egg bdist_wheel

clean:
	rm -f -r $(wildcard obj/*.o) $(wildcard obj/*.co) *.so
	rm -r -f dist *.egg-info build

docs: LAuS.html

.Phony: all clean install distrib docs