PROG=isChIP
COPT=-c #-D_NO_ZLIB# uncomment last macro if no ZLIB on your system
LOPT=-lpthread -lz# comment last option if no ZLIB on your system
SRC=$(wildcard *.cpp)
HDR=$(wildcard *.h)
OBJ=$(SRC:.cpp=.o)
EXEC=$(PROG)
CC=g++
#CC=icpc
.PHONY: all clean

all: $(HDR) $(SRC) $(EXEC)

$(EXEC): $(OBJ)
	$(CC) $(LOPT) $(OBJ) -o $@
	@echo "$(PROG) compilation complete."
#	cp $@ ..

.cpp.o:
	$(CC) $(COPT) $< -o $@

clean:
	rm *o
