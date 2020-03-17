# COMPILATEUR
cc=gcc
# ARCHITECTURE
ARCH=native
# OPTIONS DE COMPILATION
ifeq ($(ARCH), )
CFLAGS= -O3
else
CFLAGS= -O3 -march=$(ARCH) -mtune=$(ARCH)
endif

libMieGruneisen.so : miegruneisen.o
	$(cc) -shared $(CFLAGS) miegruneisen.o -o libMieGruneisen.so

miegruneisen.o : miegruneisen.c miegruneisen.h
	$(cc) -c -fPIC $(CFLAGS) miegruneisen.c

clean : 
	rm libMieGruneisen.so miegruneisen.o


