INCDIR = -I$(ROOTSYS)/include -Iinc 

TARGET = main	
LIBS = -lstdc++ -L/usr/local/lib -lm $(shell root-config --evelibs) 

CC = gcc
CFLAGS = -g -ansi -pedantic -W -Wall -Wextra -std=c++11 -m64
CFLAGS = -g -ansi -pedantic -W -Wall -Wextra -m64

SRC = $(wildcard src/*.cc)
OBJECTS = $(subst src,obj,$(patsubst %.cc,%.o,$(SRC)))
HEADERS = $(wildcard $(INCDIR)/*.h)

.PHONY: all clean

all: $(TARGET)
	@echo "DONE!"

obj/%.o: src/%.cc $(HEADERS) | obj
	$(CC) $(CFLAGS) $(INCDIR) -c $< -o $@

$(TARGET): $(OBJECTS)
	$(CC) $(CFLAGS) $(LIBS) $(OBJECTS) -o $@

obj:
	mkdir -p $@

clean:
	rm -r ./obj
	rm -r $(TARGET)
