GCC = gcc
GPP = g++
LINK = g++
CPP  = cpp
LIBTM_DIR = ./include
CC_FLAGS = -ggdb -pg #-msse2 -mfpmath=sse -O3 -g

#the program I am currently working on
all: water


water:
	$(GPP) -c $(CC_FLAGS) data.cc
	$(GPP) -c $(CC_FLAGS) bndry.cc
	$(GPP) -c $(CC_FLAGS) cnstnt.cc
	$(GPP) -c $(CC_FLAGS) cshift.cc
	$(GPP) -c $(CC_FLAGS) initia.cc
	$(GPP) -c $(CC_FLAGS) interf.cc
	$(GPP) -c $(CC_FLAGS) intraf.cc
	$(GPP) -c $(CC_FLAGS) kineti.cc
	$(GPP) -c $(CC_FLAGS) poteng.cc
	$(GPP) -c $(CC_FLAGS) predcor.cc
	$(GPP) -c $(CC_FLAGS) mdmain.cc
	$(GPP) -c $(CC_FLAGS) syscons.cc
	$(GPP) -c $(CC_FLAGS) water.cc
	$(GPP) *.o -o water -pg -lpthread -lm -l_tm -L$(LIBTM_DIR)
	rm *.o

#
clean:
	rm *.o water core*



