#=================================================
# Makefile settings
#=================================================
.PHONY: tags

.DEFAULT_GOAL=default # specify default target if there is no target specified

VPATH=library # specify directory where to check target file update to date or not

SHELL=/bin/bash

#=================================================
# Compiler settings
#=================================================
FORTRAN    = gfortran
LINKER     = $(FORTRAN)
DEFINES    = 
BUILD_TYPE = -O3
FFLAGS     = $(BUILD_TYPE) -Wfatal-errors -Wall -Wextra -std=f95 -fPIC 
LIBS       = 
INCLUDE    = 
MOVE       = mv -f

#=================================================
# Build target
#=================================================
LIBRARY_DIR       := library
BIN_DIR       	  := bin
OBJS              := main.o utils.o
TARGET            := test

#=================================================
# Compile implicit rules
#=================================================
%.o:%.f
	$(FORTRAN) -c $< -o $@ $(FFLAGS) $(INCLUDE)
	@$(MOVE) $@ $(LIBRARY_DIR) # @ means do not echo command

#=================================================
# Target rules
#=================================================
default: $(TARGET)

$(TARGET): create_dir tags $(OBJS) 
	cd $(LIBRARY_DIR);\
	$(LINKER) -o $@ $(OBJS) $(LIBS);\
	$(MOVE) $@ ../$(BIN_DIR)

debug: BUILD_TYPE:=-g3 
debug: default

create_dir:
	@mkdir -p $(LIBRARY_DIR)
	@mkdir -p $(BIN_DIR)

clean:
	rm -rf $(LIBRARY_DIR)/*

tags: 
	rm -f tags
	find . -name "*.f" | xargs ctags --language-force=fortran -f tags
