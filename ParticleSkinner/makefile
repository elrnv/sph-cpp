#!smake 

#-----------------------------------------
#Basic Stuff -----------------------------
CC          = g++ 
cc          = gcc

#-----------------------------------------
#Misc Flags ------------------------------
MISC_FLAGS = -Wall

#-----------------------------------------
#Optimization ----------------------------
OPT   = -O3  #-fast

#-----------------------------------------
# Debug mode -----------------------------
DEBUG = -g -pg

LIB_MOD = 

TARGETS = particleskinner

#-----------------------------------------
#-----------------------------------------

OBJECTS = main.o smoothingGrid.o marchingTet.o kdTree.o slMatrix.o slVector.o slUtil.o

#-----------------------------------------

LIBS =  
INCS = 

CCOPTS = $(OPT) $(INCS)
LDOPTS = $(OPT) $(LIBS)

#-----------------------------------------
#-----------------------------------------

default: $(TARGETS)



clean: 
	/bin/rm -rf *.o particleskinner

#-----------------------------------------

particleskinner: $(OBJECTS)
	$(CC) $(LDOPTS) $(OBJECTS) $(LIBS) -o particleskinner

#-----------------------------------------
.C.o: 
	$(CC) $(CCOPTS) -c $< 

.cpp.o: 
	$(CC) $(CCOPTS) -c $< 

.c.o: 
	$(cc) $(CCOPTS) -c $< 

.o: $(OBJECTS)
	$(CC)  $(LDOPTS) $(OBJS) $(OBJECTS) $< $(LIBS) -o $@ 

.C: $(OBJECTS)
	$(CC) $(LDOPTS)   $(OBJS) $(OBJECTS) $< $(LIBS) $(FOR_LIB) -o $@

#-----------------------------------------
#-----------------------------------------














