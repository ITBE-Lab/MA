# location of the Boost Python include files and library
BOOST_LIB_PATH = $(BOOST_ROOT)/stage/lib/
BOOST_LIB = boost_python3 boost_iostreams boost_filesystem boost_system boost_program_options
CUDA_PATH = /usr/local/cuda-9.1/lib64/
 
# target files
TARGET = $(subst .cpp,,$(subst src/,,$(wildcard src/*/*.cpp)))

TARGET_OBJ=$(addprefix obj/,$(addsuffix .o,$(TARGET))) \
	obj/ksw/ksw2_dispatch.co obj/ksw/ksw2_extz2_sse2.co obj/ksw/ksw2_extz2_sse41.co \
	obj/container/qSufSort.co

DEBUG_OBJ=$(addprefix dbg/,$(addsuffix .o,$(TARGET))) \
	obj/ksw/ksw2_dispatch.co obj/ksw/ksw2_extz2_sse2.co obj/ksw/ksw2_extz2_sse41.co \
	obj/container/qSufSort.co

#flags
CC=gcc
CCFLAGS=-Wall -DBOOST_ALL_DYN_LINK -Werror -fPIC -std=c++11 -mavx2 -O3 \
	-DWITH_PYTHON # this flag enables the python library
CFLAGS=-Wall -Werror -fPIC -O3
LDSFLAGS=-shared -Wl,--export-dynamic
LDFLAGS= -std=c++11
LDLIBS= \
	$(PYTHON_LIB) \
	-L$(BOOST_LIB_PATH) \
	-L$(PARSAIL_HOME)/build \
	-L$(LIBGABA_HOME) \
	$(addprefix -l,$(addsuffix $(BOOST_SUFFIX),$(BOOST_LIB))) \
	-L$(CUDA_PATH) \
	-lm \
	-lpthread \
	-lstdc++ \
	-lparasail \
	-lcudart \
	-lgaba

all: ma

debug: TARGET_OBJ=$(addprefix dbg/,$(addsuffix .o,$(TARGET))) \
	obj/ksw/ksw2_dispatch.co obj/ksw/ksw2_extz2_sse2.co obj/ksw/ksw2_extz2_sse41.co \
	obj/container/qSufSort.co

debug: CFLAGS = -Wall -Werror -fPIC -g -DDEBUG_LEVEL=1
debug: $(DEBUG_OBJ)
	$(CC) $(LDFLAGS) $(LDSFLAGS) -g $(DEBUG_OBJ) $(LDLIBS) sw_gpu.o -o libMA.so

ma: libMA src/cmdMa.cpp
	$(CC) $(CCFLAGS) src/cmdMa.cpp -isystem$(PYTHON_INCLUDE)/ -isystem$(BOOST_ROOT)/ -isystem$(PARSAIL_HOME)/ -isystem$(LIBGABA_HOME)/ -Iinc $(LDLIBS) libMA.so -o $@

libMA: $(TARGET_OBJ)
	$(CC) $(LDFLAGS) $(LDSFLAGS) $(TARGET_OBJ) $(LDLIBS) sw_gpu.o -o libMA.so

#special targets for the ksw2 library
obj/ksw/ksw2_dispatch.co:src/ksw/ksw2_dispatch.c inc/ksw/ksw2.h
		$(CC) -c $(CFLAGS) -Iinc -DKSW_CPU_DISPATCH $< -o $@

obj/ksw/ksw2_extz2_sse2.co:src/ksw/ksw2_extz2_sse.c inc/ksw/ksw2.h
		$(CC) -c $(CFLAGS) -Iinc -msse2 -mno-sse4.1 -DKSW_CPU_DISPATCH -DKSW_SSE2_ONLY $< -o $@

obj/ksw/ksw2_extz2_sse41.co:src/ksw/ksw2_extz2_sse.c inc/ksw/ksw2.h
		$(CC) -c $(CFLAGS) -Iinc -msse4.1 -DKSW_CPU_DISPATCH $< -o $@

obj/container/qSufSort.co:src/container/qSufSort.c inc/container/qSufSort.h
		$(CC) -c $(CFLAGS) -Iinc $< -o $@

dbg/%.o: src/%.cpp inc/%.h
	$(CC) -Wall -DBOOST_ALL_DYN_LINK -Werror -fPIC -std=c++11 -mavx2 -g -DDEBUG_LEVEL=1 -isystem$(PYTHON_INCLUDE)/ -isystem$(BOOST_ROOT)/ -isystem$(PARSAIL_HOME)/ -isystem$(LIBGABA_HOME)/ -Iinc -c $< -o $@

obj/%.o: src/%.cpp inc/%.h
	$(CC) $(CCFLAGS) -isystem$(PYTHON_INCLUDE)/ -isystem$(BOOST_ROOT)/ -isystem$(PARSAIL_HOME)/ -isystem$(LIBGABA_HOME)/ -Iinc -c $< -o $@

html/index.html: $(wildcard inc/*) $(wildcard inc/*/*) $(wildcard src/*) $(wildcard src/*/*) $(wildcard MA/*.py) doxygen.config
	doxygen doxygen.config

install: all
	#pip3 install . --upgrade --no-cache-dir #no pip installation at the moment
	cp libMA.so /usr/lib
	#pip3 show MA #no pip installation at the moment

##@todo remove me
vid:
	gource -f --seconds-per-day 0.1

distrib:
	python setup.py sdist bdist_egg bdist_wheel

clean:
	rm -f -r obj/*.o dbg/*.o obj/*/*.o dbg/*/*.o obj/*.co dbg/*.co obj/*/*.co dbg/*/*.co libMA.so
	rm -r -f dist *.egg-info build
	rm -r -f html

docs: html/index.html

.Phony: all clean install distrib docs vid libMA debug