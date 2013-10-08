CC=g++
CPPFLAG	:= -g
pdc: main.cpp add.h add.cpp
	$(CC) $(CPPFLAG)  $^ -o pdc
clean:
	rm -rf pdc
