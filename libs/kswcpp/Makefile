
# target files
STATIC=lib/libkswcpp.a

TARGET_OBJ= \
	obj/cpu_info.o \
	obj/kswcpp_avx2.o \
	obj/kswcpp_sse_xx.o

# flags
CC=gcc
CCFLAGS= -Wall -O3 -std=c++17 -fPIC
LDFLAGS= -std=c++17

# primary target
all: dirs lib/libkswcpp

# create build directories if not present
dirs:
	mkdir -p obj lib

obj/cpu_info.o:
	$(CC) -c $(CCFLAGS) -Iinc src/cpu_info.cpp -o $@

obj/kswcpp_avx2.o:
	$(CC) -c $(CCFLAGS) -Iinc src/kswcpp_avx2.cpp -Wno-sign-compare -Wno-strict-aliasing -mavx2 -o $@

obj/kswcpp_sse_xx.o:
	$(CC) -c $(CCFLAGS) -Iinc src/kswcpp_sse_xx.cpp -Wno-sign-compare -Wno-strict-aliasing -msse4.1 -o $@

# library target
lib/libkswcpp: $(TARGET_OBJ)
	@echo "[Link (Static)]" $(TARGET_OBJ)
	ar -rcs -o lib/libkswcpp.a $(TARGET_OBJ)

clean:
	rm -f -r obj ma libkswcpp.a

.Phony: all dirs clean