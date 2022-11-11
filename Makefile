CC = g++
CFLAGS = -Wall -g -std=c++14
PACKAGE = `pkg-config --cflags`
LIBS = `pkg-config`
INC=-I./include

default:	main

main:	src/Main.cpp
	$(CC) $(PACKAGE) $(CFLAGS) $(INC) -o build/kontig src/*.cpp

clean:
	rm -r build/*
