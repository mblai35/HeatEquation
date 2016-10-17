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
# 	
####################################################################

CC 	= gcc
LFLAGS 	= -g -O3

TARGET	= theta1D

Geeta 	=c 81.09 -.09036 92.93 -.002168 b 80.35 -.1156 93.69 -.002442
Mallory =c 90.49 -.06361 80.13 -.001023 b 73.97 -.08249 81.56 -.001303
Xiukun  =

ARG 	=Xiukun

all:		$(TARGET)

run:		all
	./$(TARGET) $($(ARG))

$(TARGET): 	$(TARGET).o Makefile
	$(CC) $(LFLAGS) -o $(TARGET) $(TARGET).o

%.o:	%.c Makefile
	$(CC) $(LFLAGS) -c -o $@ $<

clean:
	- /bin/rm -f *.o

distclean:
	- /bin/rm -rf *.o Theta1D.txt *.dSYM $(TARGET)
