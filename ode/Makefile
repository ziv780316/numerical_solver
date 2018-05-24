#=================================================
# Makefile settings
#=================================================
-include dependency.mk # - do not abort if file does not exist

.PHONY: create_dir verbose_info depend default all clean cscope qa# specify phony target, do not generate file as target name

.DEFAULT_GOAL=default # specify default target if there is no target specified

VPATH=include library bin unit_test # specify directory where to check target file update to date or not

SHELL=/bin/bash

#=================================================
# Compiler settings
#=================================================
CC         = gcc
CXX        = g++
FLEX       = flex 
FLEX_FLAG  = 
BISON      = bison 
LINKER     = $(CXX)
DEFINES    = 
IVP        = IVP1
BUILD_TYPE = -O3
CFLAGS     = $(BUILD_TYPE) -Wfatal-errors -Wall -std=gnu99 -D$(IVP) -fPIC 
CXXFLAGS   = $(BUILD_TYPE) -Wfatal-errors -Wall -std=c++11 -D$(IVP) -fPIC
LIBS       = -lm -L/usr/lib/x86_64-linux-gnu -lfl -lreadline
INCLUDE    = -I./include
MOVE       = mv -f

#=================================================
# Build target
#=================================================
BIN_DIR           := bin
UNIT_TEST_DIR     := unit_test
LIBRARY_DIR       := library
ODE_SOLVER_OBJ    := adams_bashforth.o adams_moulton.o bdf.o opts.o diff.o main.o
LEX_YACC_TARGET   := 
UNIT_TEST_OBJ     := 
SO_TARGET         :=
ODE_SOLVER_TARGET := ode_solver
UNIT_TEST_TARGET  := 

#=================================================
# Compile implicit rules
#=================================================
%.o:%.c
	$(CC) -c $< -o $@ $(CFLAGS) $(INCLUDE)
	@$(MOVE) $@ $(LIBRARY_DIR) # @ means do not echo command

%.e:%.c
	$(CC) -E $< -o $@ $(CFLAGS) $(INCLUDE)
	@$(MOVE) $@ $(LIBRARY_DIR)

%.o:%.cxx
	$(CXX) -c $< -o $@ $(CXXFLAGS) $(INCLUDE)
	@$(MOVE) $@ $(LIBRARY_DIR)

%.e:%.cxx
	$(CXX) -E $< -o $@ $(CXXFLAGS) $(INCLUDE)
	@$(MOVE) $@ $(LIBRARY_DIR)

%.l.c:%.l
	$(FLEX) $(FLEX_FLAG) -o $@ $< 

%.y.c:%.y
	$(BISON) -v -y -Wall -o $@ $<

%.lxx.c:%.lxx
	$(FLEX) $(FLEX_FLAG) --c++ -o $@ $< 

%.yxx.c:%.yxx
	$(BISON) -v -y -Wall -o $@ $<

#=================================================
# Target rules
#=================================================
default: $(ODE_SOLVER_TARGET) $(SO_TARGET) 

all: $(ODE_SOLVER_TARGET) $(SO_TARGET) $(UNIT_TEST_TARGET)

debug: BUILD_TYPE:=-g3 
debug: default

create_dir:
	@mkdir -p $(LIBRARY_DIR) $(BIN_DIR) $(UNIT_TEST_DIR)

verbose_info:
	@echo IVP=$(IVP)

$(ODE_SOLVER_TARGET): verbose_info create_dir $(LEX_YACC_TARGET) $(ODE_SOLVER_OBJ)
	cd $(LIBRARY_DIR);\
	$(LINKER) -o $@ $(ODE_SOLVER_OBJ) $(LIBS);\
	$(MOVE) $@ ../$(BIN_DIR)

$(UNIT_TEST_TARGET): create_dir $(UNIT_TEST_OBJ)
	cd $(LIBRARY_DIR);\
	cd $(LIBRARY_DIR);\
	$(LINKER) -o $@ $(UNIT_TEST_OBJ) $(LIBS);\
	$(MOVE) $@ ../$(UNIT_TEST_DIR)

clean:
	rm -f dependency.mk
	rm -rf $(LIBRARY_DIR)/*

depend:
	rm -f dependency.mk
	for file in `find . -name "*.c" -o -name "*.cxx" -o -name "*.h" -o -name "*.hxx"`; do $(CC) -M $(INCLUDE) $$file >> dependency.mk; done

cscope:
	find . -name "*.c" -o -name "*.cxx" -o -name "*.h" -o -name "*.hxx" > cscope.files
	cscope -Rbq -i cscope.files

qa:
	./bin/$(ODE_SOLVER_TARGET) > output


