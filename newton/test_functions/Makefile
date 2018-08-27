#=================================================
# Makefile settings
#=================================================
.PHONY: clean # specify phony target, do not generate file as target name

.DEFAULT_GOAL=default # specify default target if there is no target specified

VPATH=../library # specify directory where to check target file update to date or not

SHELL=/bin/bash

#=================================================
# Compiler settings
#=================================================
CC         = gcc
LINKER     = $(CC)
DEFINES    = 
BUILD_TYPE = -O3
CFLAGS     = $(BUILD_TYPE) -Wfatal-errors -Wall -std=gnu99 -fPIC 
LIBS       = -lm
INCLUDE    = -I.
MOVE       = mv -f

#=================================================
# Build target
#=================================================
TEST_FUNCTION_OBJ := f1.so f2.so f3.so

#=================================================
# Compile implicit rules
#=================================================
%.so:%.c
	$(CC) -shared $< -o $@ $(CFLAGS) $(INCLUDE) 


#=================================================
# Target rules
#=================================================
default: $(TEST_FUNCTION_OBJ)

clean:
	rm -rf $(TEST_FUNCTION_OBJ)


