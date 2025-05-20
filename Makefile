# define the C compiler to use
CC=gcc

# define any compile-time flags
CFLAGS=-Wall -Ofast -funroll-loops 


# define any libraries to link into executable:
LIBS=-lm


md: EDMC_poly.c  EDMC_poly.h
	$(CC) $(CFLAGS) EDMC_poly.c -o md  $(LIBS)
