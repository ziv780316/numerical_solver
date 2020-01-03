#=================================================
# Makefile settings
#=================================================
.PHONY: clean # specify phony target, do not generate file as target name

.DEFAULT_GOAL=default # specify default target if there is no target specified

VPATH=library # specify directory where to check target file update to date or not

SHELL=/bin/bash

#=================================================
# Compiler settings
#=================================================
CC         = gcc
LINKER     = $(CC)
DEFINES    = 
BUILD_TYPE = -O0
CFLAGS     = $(BUILD_TYPE) -Wfatal-errors -Wall -Wextra -std=gnu99 -fPIC -Wno-unused-parameter
LIBS       = -lm
INCLUDE    = -I.
MOVE       = mv -f

#=================================================
# Build target
#=================================================
TEST_FUNCTION_OBJ := r_net1.so r_net1_partition.so \
		     r_net2.so r_net2_partition.so \
		     ir_net1.so ir_net1_partition.so \
		     ir_net2.so ir_net2_partition.so \
		     ids_net1.so ids_net1_partition.so \

#=================================================
# Compile implicit rules
#=================================================
%.so:%.c
	$(CC) -shared $< -o $@ $(CFLAGS) $(INCLUDE) ${LIBS}


#=================================================
# Target rules
#=================================================
default: $(TEST_FUNCTION_OBJ)

debug: BUILD_TYPE:=-g3 
debug: default

clean:
	rm -rf $(TEST_FUNCTION_OBJ)


