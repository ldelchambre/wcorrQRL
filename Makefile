CC=gcc
CFLAGS=-W -Wall
INCLUDES=-I./includes/
LDFLAGS=
LDLIBS=-lm -lfftw3
EXEC=demo
OBJ=
BASE_DIR=src

.PHONY: clean dist-clean

all: $(EXEC)

base:
	@$(MAKE) -C $(BASE_DIR)

%.o : %.c
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $^
	
clean:
	@make -C $(BASE_DIR) clean
	@rm -rf $(OBJ) *.o
	
dist-clean: clean
	@make -C $(BASE_DIR) dist-clean
	@rm -rf $(EXEC)

demo : $(OBJ) | base
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ demo.c $(BASE_DIR)/*.o $^ $(LDFLAGS) $(LDLIBS)
	
run: all
	./demo
	
doc:
	@doxygen
