SHELL = /bin/sh
CC    = gcc
FLAGS        = # -std=gnu99 -Iinclude
CFLAGS       = -fPIC -g # -pedantic -Wall -Wextra -march=native -ggdb3
LDFLAGS      = -shared
DEBUGFLAGS   = -O0 -D _DEBUG
RELEASEFLAGS = -O2 -D NDEBUG -combine -fwhole-program
INCLUDEFLAGS = -lm

TARGET  = libemath.so
SOURCES = $(shell echo src/*.c)
HEADERS = $(shell echo include/*.h)
OBJECTS = $(SOURCES:.c=.o)

PREFIX = $(DESTDIR)/usr/local
BINDIR = $(PREFIX)/bin

all: $(TARGET)

$(TARGET): $(OBJECTS) 
	$(CC) $(FLAGS) $(CFLAGS) $(DEBUGFLAGS) -o $(TARGET) $(OBJECTS) $(INCLUDEFLAGS)
