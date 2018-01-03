CC      = g++
CFLAGS  = -Wall -O2
LIBS    = -L/usr/local/lib -lnort -lm
.PHONY : all test
all :
	eli 'front%dtarget >icosfit.exe'
test :
	eli 'front%dtarget >/usr/local/bin/icosfit_test.exe'

include front/Makefile.tail
