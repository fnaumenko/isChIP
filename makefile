PROG=isChIP
COPT=-c -O3
LOPT=-lpthread
ifeq ($(shell whereis zlib | wc -l),0)
	COPT+= -D_NO_ZLIB
	WARNING="WARNING: zlib is not installed! isChIP will not read/write .gz files!\n---------------------------------------------------------------------\n"
else
	LOPT+= -lz
endif
SRC_DIR=src
SRC=$(wildcard $(SRC_DIR)/*.cpp)
HDR=$(wildcard $(SRC_DIR)/*.h)
OBJ=$(SRC:.cpp=.o)
EXEC=$(PROG)
CC=g++
#CC=icpc
.PHONY: all clean

all: print_warning $(HDR) $(SRC) $(EXEC)

$(EXEC): $(OBJ)
	$(CC) $(LOPT) $(OBJ) -o $@
	@echo "$(PROG) compilation complete."

.cpp.o:
	$(CC) $(COPT) $< -o $@

print_warning:
	@echo -n -e $(WARNING)
.PHONY: print_warning

clean:
	rm $(SRC_DIR)/*o
