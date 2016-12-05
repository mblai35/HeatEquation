####################################################################
#
# Makefile for CompMethII homework. For now it only works for 
# theta1d.c
#
# Group Member: Mallory, Geeta, Xiukun
#
# Modification:
# 	None
#
# How to use this:
# 	
# 	make 			Compile theta1D.c.
# 	make run [ARG=<Name>]	Compile and run with specific data.
# 	make clean 		Delete all object file *.o
# 	make distclean		Delete all object file *.o, output 
# 				files, debug files and executables.
# 	
####################################################################

PROJECT = CM2hw2

TFILE 	= $(PROJECT).tgz
PACKFILES = $(TARGET1D).c $(TARGET2D).c HeatEquation.c HeatEquation.h Makefile

CC 	= gcc
LFLAGS 	= -g -O3 -lm

TARGET1D	= theta1D
TARGET2D	= theta2D

DATA 	= XIUKUN
LENGTH  = 3.0
HEIGHT  = .5
TIME 	= 70.0
ALPHA 	= 1.0

MACRO	= -D L=$(LENGTH) -D H=$(HEIGHT) -D T=$(TIME) -D ALPHA=$(ALPHA) -D $(DATA) 

all:		$(TARGET1D) $(TARGET2D)

1D:		$(TARGET1D)
	./$(TARGET1D)

2D:		$(TARGET2D)
	./$(TARGET2D)

run:		1D 2D

$(TARGET1D): 	$(TARGET1D).o HeatEquation.o Makefile
	$(CC) $(LFLAGS) -o $(TARGET1D) $(TARGET1D).o HeatEquation.o

$(TARGET2D):	$(TARGET2D).o HeatEquation.o Makefile
	$(CC) $(LFLAGS) -o $(TARGET2D) $(TARGET2D).o HeatEquation.o

HeatEquation.o:	HeatEquation.c Makefile
	$(CC) $(LFLAGS) $(MACRO) -c -o HeatEquation.o HeatEquation.c
	
%.o:	%.c Makefile
	$(CC) $(LFLAGS) $(MACRO) -c -o $@ $<

clean:
	- /bin/rm -f *.o $(TARGET1D) $(TARGET2D) $(TFILE) *.dSYM

distclean: clean
	- /bin/rm -f theta*.txt theta*.bin

pack:
	tar -zcvf $(TFILE) $(PACKFILES)
