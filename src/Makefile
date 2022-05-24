# Makefile for HiGHS package for Octave
#
# Copyright (C) 2022 Markus Mützel

MKOCTFILE ?= mkoctfile
PKG_CONFIG ?= pkg-config

CXX ?= $(shell $(MKOCTFILE) -p CXX)
CXX += -g3 -std=gnu++11

all:
	CXX="$(CXX)" $(MKOCTFILE) $(shell $(PKG_CONFIG) --cflags highs) $(shell $(PKG_CONFIG) --libs highs) __highs__.cc
