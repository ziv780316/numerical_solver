#=================================================
# Makefile settings
#=================================================
.PHONY: default create_dir clean # specify phony target, do not generate file as target name

.DEFAULT_GOAL=default # specify default target if there is no target specified

VPATH=library # specify directory where to check target file update to date or not

SHELL=/bin/bash

#=================================================
# Compiler settings
#=================================================
CC         = gcc
CPP        = g++
LINKER     = $(CPP)
CPPLINKER  = $(CPP)
DEFINES    = 
BUILD_TYPE = -O0 -g
CFLAGS     = $(BUILD_TYPE) -Wfatal-errors -std=gnu99 -fPIC 
CPPFLAGS   = $(BUILD_TYPE) -Wfatal-errors -std=c++98 -fPIC
LIBS       = -lm -ldl
INCLUDE    = -I.
MOVE       = mv -f

#=================================================
# Build target
#=================================================
BIN_DIR       	  := bin
LIB_DIR       	  := library
OBJS              := main.o rcr.o util.o
ASMS              := main.s rcr.s util.s
HEADER            := main.i rcr.i util.i
TARGET            := rcr

#=================================================
# Compile implicit rules
#=================================================
%.o:%.c
	$(CC) $(INCLUDE) $(LIBS) $(CFLAGS) -c $< -o $@
	$(MOVE) $@ $(LIB_DIR)

%.o:%.cpp
	$(CPP) $(INCLUDE) $(LIBS) $(CPPFLAGS) -c $< -o $@
	$(MOVE) $@ $(LIB_DIR)

%.i:%.c
	$(CC) -E $< -o $@ $(CFLAGS) $(INCLUDE)
	@$(MOVE) $@ $(LIB_DIR) # @ means do not echo command

%.i:%.cpp
	$(CPP) -E $< -o $@ $(CPPFLAGS) $(INCLUDE)
	@$(MOVE) $@ $(LIB_DIR) # @ means do not echo command

%.s:%.c
	$(CC) -S $< -o $@ $(CFLAGS) $(INCLUDE)
	@$(MOVE) $@ $(LIB_DIR) # @ means do not echo command

%.s:%.cpp
	$(CPP) -S $< -o $@ $(CPPFLAGS) $(INCLUDE)
	@$(MOVE) $@ $(LIB_DIR) # @ means do not echo command

#=================================================
# Target rules
#=================================================
$(TARGET): $(OBJS) $(HEADER) $(ASMS)
	cd $(LIB_DIR);\
	$(LINKER) -o $@ $(OBJS) $(LIBS);\
	$(MOVE) $@ ../$(BIN_DIR);\
	cd ..;\
	objdump -dC $(BIN_DIR)/$(TARGET) > $(TARGET).asm;\
	objdump -tTC $(BIN_DIR)/$(TARGET) > $(TARGET).sym

#=================================================
# Phony rules
#=================================================
default: create_dir $(TARGET)
	

create_dir:
	@echo "create dir ..."
	@mkdir -p $(LIB_DIR)
	@mkdir -p $(BIN_DIR)

clean:
	rm -rf $(LIB_DIR)/* $(BIN_DIR)/*

cscope:
	find . -name "*.c" -o -name "*.h" > cscope.files 
	cscope -Rbq -i cscope.files

test: 
	@echo "---------------------------------------"
	@echo "TEST"
	@echo "---------------------------------------"
	valgrind --error-limit=no ./$(BIN_DIR)/$(TARGET) >& test.log;\
	@echo "---------------------------------------"
	@echo "DONE"
	@echo "---------------------------------------"

	
	



