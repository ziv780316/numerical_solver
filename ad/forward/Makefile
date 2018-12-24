#=================================================
# Makefile settings
#=================================================
-include dependency.mk # - do not abort if file does not exist

.PHONY: create_dir depend default clean cscope qa# specify phony target, do not generate file as target name

.DEFAULT_GOAL=default # specify default target if there is no target specified

VPATH=include library bin unit_test # specify directory where to check target file update to date or not

#=================================================
# Specify package config path
#=================================================


#=================================================
# Compiler settings
#=================================================
CC         = gcc
CPP        = g++
FLEX       = flex 
BUILD_TYPE = -O2
FLEX_BUILD_TYPE = 
FLEX_FLAG  = $(FLEX_BUILD_TYPE)
DEFINES    =
YACC       = bison -Wall
YACC_FLAG  = 
LINKER     = $(CC)
CFLAGS     = $(DEFINES) $(BUILD_TYPE) -Wfatal-errors -Wall -Wextra -std=gnu99 -fPIC
CPPFLAGS   = $(DEFINES) $(BUILD_TYPE) -Wfatal-errors -Wall -Wextra -std=c++11 -fPIC
LIBS       = -lm 
INCLUDE    = -I./include 
MOVE       = mv -f

#=================================================
# Build target
#=================================================
BIN_DIR          = bin
UNIT_TEST_DIR    = unit_test
LIBRARY_DIR      = library
LEX_YACC_TARGET  = 
LEX_YACC_OBJ     = $(LEX_YACC_TARGET:.c=.o)
COMPLEX_AD_OBJS   = main.o 
COMPLEX_AD_TARGET = test

#=================================================
# Compile implicit rules
#=================================================
%.o:%.c
	$(CC) -c $< -o $@ $(CFLAGS) $(INCLUDE)
	@$(MOVE) $@ $(LIBRARY_DIR) # @ means do not echo command

%.e:%.c
	$(CC) -E $< -o $@ $(CFLAGS) $(INCLUDE)
	@$(MOVE) $@ $(LIBRARY_DIR)

%.o:%.cpp
	$(CPP) -c $< -o $@ $(CPPFLAGS) $(INCLUDE)
	@$(MOVE) $@ $(LIBRARY_DIR) # @ means do not echo command

%.e:%.cpp
	$(CPP) -E $< -o $@ $(CPPFLAGS) $(INCLUDE)
	@$(MOVE) $@ $(LIBRARY_DIR) # @ means do not echo command

%.l.o:%.l.c
	$(CC) -c $< -o $@ $(CFLAGS) $(INCLUDE)
	@$(MOVE) $@ $(LIBRARY_DIR) # @ means do not echo command

%.l.c:%.l
	$(FLEX) $(FLEX_FLAG) -o $@ $< 

%.y.c:%.y
	$(YACC) $(YACC_FLAG) -v -y -Wall -o $@ $<

#=================================================
# Target rules
#=================================================
default: $(COMPLEX_AD_TARGET) 

debug: BUILD_TYPE:=-g3 
debug: FLEX_BUILD_TYPE:=-d
debug: default

create_dir:
	@mkdir -p $(LIBRARY_DIR) $(BIN_DIR) 

$(COMPLEX_AD_TARGET): create_dir $(COMPLEX_AD_OBJS)
	cd $(LIBRARY_DIR);\
	$(LINKER) -o $@ $(COMPLEX_AD_OBJS) $(LIBS);\
	$(MOVE) $@ ../$(BIN_DIR)

clean:
	rm -f dependency.mk
	rm -rf $(LIBRARY_DIR)/*

depend:
	rm -f dependency.mk
	for file in $(shell find . -name "*.c" -o -name "*.cxx" -o -name "*.h" -o -name "*.hxx"); do $(CC) -M $(INCLUDE) $$file >> dependency.mk; done

cscope:
	find . -name "*.c" -o -name "*.cxx" -o -name "*.h" -o -name "*.hxx" > cscope.files
	cscope -Rbq -i cscope.files

qa:

