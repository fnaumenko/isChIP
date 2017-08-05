#.PHONY: all clean

PROG=isChIP
CC=g++
#CC=icpc
COPT=-c #-D_NO_HASHTABLE #-Wall
LOPT=-lpthread -lz
SRC=$(PROG).cpp ImitFile.cpp Imitator.cpp common.cpp Data.cpp TxtFile.cpp
HDR=$(PROG).h def.h ImitFile.h Imitator.h common.h Data.h TxtFile.h
OBJ=$(SRC:.cpp=.o)
EXEC=$(PROG)

all: $(HDR) $(SRC) $(EXEC)

$(EXEC): $(OBJ)
	$(CC) $(LOPT) $(OBJ) -o $@
	cp $@ ..

.cpp.o:
	$(CC) $(COPT) $< -o $@

clean:
	rm *o
