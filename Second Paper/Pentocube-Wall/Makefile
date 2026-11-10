.PHONY: all venv

CC ?= cc
CFLAGS ?= -O2 -Wall -Wextra

all: ssxcc

venv:
	python3 -m venv venv
	venv/bin/python -m pip install polycubetools numpy

ssxcc: src/ssxcc.c src/gb_flip.c src/gb_flip.h
	$(CC) $(CFLAGS) -o $@ src/ssxcc.c src/gb_flip.c
