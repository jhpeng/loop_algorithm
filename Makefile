# 'make'	build executable file
# 'make clean'	removes all *.o and executalbe file

# define the C compiler
APP = hb_afm_2d_normal

CC	= gcc

# define any compile-time flags
CFLAGS = -Wall -g -fPIC -O0 -std=c99

# define openmp flags
OPENMP  = -fopenmp
CUOPENMP  = -Xcompiler -fopenmp

# define the direction containing header file
INCLUDES= -I/usr/local/include -I./include

# define the library path
LFLAGS	= -L/usr/local/lib -L./lib

# define any libraries to link to executable
LIBS	= -lm -lgsl -lgslcblas -lloop


#define the directory for object
OBJSDIR = object

# define the executable file
MAIN	= $(APP)

all: $(MAIN)

$(MAIN): $(APP).o
	$(CC) $(CFLAGS) -o $(MAIN) $(APP).o $(LIBS) $(LFLAGS) $(INCLUDES) 
	rm $(APP).o 

$(APP).o: $(APP).c
	$(CC) $(CFLAGS) $(INCLUDES) -c $(APP).c


# clean the executable file and object files
clean:
	$(RM) $(APP).o $(MAIN)
