CC      = gcc
CFLAGS  = -Wall -Werror 
CFLAGS += -g
## For large files support :
CFLAGS += -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE 

LDFLAGS = -lm
BINARIES =  count_kmer
OBJECTS = count_lib.o count_kmer.o count_read.o count_print.o

all : $(BINARIES) $(OBJECTS)

###------------------------------
### Entries
###------------------------------------------------------------
count_kmer: $(OBJECTS) 
	$(CC)  -o $@ $^ $(CFLAGS) $(LDFLAGS)

###------------------------------
### Dependencies
###------------------------------------------------------------
# thanks to 'make depend'
count_kmer.o: count_kmer.c count_lib.h count_read.h count_print.h
count_lib.o: count_lib.c count_lib.h
count_read.o: count_read.c count_read.h count_lib.h
count_print.o: count_print.c count_print.h count_lib.h

###------------------------------
### Misc.
###------------------------------------------------------------
.PHONY: clean realclean depend
clean::
	$(RM) $(BINARIES) $(OBJECTS)
realclean:: clean 
	$(RM) $(OBJECTS) 
depend: 
	$(CC) $(CFLAGS) -MM *.c 

