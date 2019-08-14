#=================================================
# Makefile settings
#=================================================
.PHONY: create_dir # specify phony target, do not generate file as target name

.DEFAULT_GOAL=default # specify default target if there is no target specified

VPATH=. # specify directory where to check target file update to date or not

SHELL=/bin/bash

#=================================================
# Compiler settings
#=================================================
CC         = gcc
LINKER     = $(CC)
PRECISION  = 80
DEFINES    = -DPRECISION_$(PRECISION)_BIT
BUILD_TYPE = -O2
CFLAGS     = $(BUILD_TYPE) $(DEFINES) -Wfatal-errors -Wall -std=gnu99 -fPIC 
MAT_SOLVER_PATH = 
LIBS       = lm 
INCLUDE    = -I. 
MOVE       = mv -f

#=================================================
# Build target
#=================================================
LIBRARY_DIR	  := library
BIN_DIR       	  := bin
SO_NAME           := numerical_util
TARGET            := $(SO_NAME) test

#=================================================
# Compile implicit rules
#=================================================
numerical_util : numerical_util.c numerical_util.i 
	$(CC) -shared $< -o $@ $(CFLAGS) $(INCLUDE) 
	@$(MOVE) $@ $(LIBRARY_DIR)/$(SO_NAME)_$(PRECISION).so
test: test.o
	@$(MOVE) $^ $(LIBRARY_DIR)
	$(LINKER) $(LIBRARY_DIR)/$(notdir $@.o) -o $(BIN_DIR)/$@ $(LIBRARY_DIR)/$(SO_NAME)_$(PRECISION).so

%.o:%.c
	$(CC) -c $< -o $@ $(CFLAGS) $(INCLUDE) 

%.i:%.c
	$(CC) -E $< -o $@ $(CFLAGS) $(INCLUDE)
	@$(MOVE) $@ $(LIBRARY_DIR)

#=================================================
# Target rules
#=================================================
default: create_dir $(TARGET)

debug: BUILD_TYPE:=-g3 
debug: default

create_dir:
	@mkdir -p $(LIBRARY_DIR) $(BIN_DIR)

clean:
	rm -f $(LIBRARY_DIR)/*


