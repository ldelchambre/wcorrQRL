CC=gcc
CFLAGS=-W -Wall
INCLUDES=-I../includes
SRC= fft.c ldl.c wchi2.c wcorr.c
OBJ= $(SRC:.c=.o)

.PHONY: clean dist-clean

all: $(OBJ)

%.o : %.c
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $^
	
clean:
	@rm -rf $(OBJ)
	
dist-clean: clean 
