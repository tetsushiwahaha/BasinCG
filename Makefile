INCLUDES = -I/usr/local/include -I/opt/local/include
CFLAGS = -O2 -I. $(INCLUDES)
LIBS = -L/opt/local/lib
CC = gcc

all: Basin BasinPutColor

.c.o:
	$(CC) $(CFLAGS) -c $<
	

Basin: Basin.o BasinZrw.o BasinUtils.o BasinFunc.o
	$(CC) $(CFLAGS) $(LIBS) -o Basin Basin.o BasinZrw.o BasinFunc.o BasinUtils.o -lz -lm

BasinPutColor:	BasinPutColor.o BasinZrw.o
	$(CC) $(CFLAGS) $(LIBS) -o BasinPutColor BasinPutColor.o BasinZrw.o BasinFunc.o -lz -lpng  -lm

clean:
	\rm *.o *.bak *core
