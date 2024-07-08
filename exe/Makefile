ESCAPE_HOME := ../

TARGETS := sample generateInput generateL0seed createCGraph generateL0COO queryOnCore PLL generateSparsifiedCOO

OBJECTS := $(TARGETS:%=%.o)

include ../common.mk

INCLUDES += -I..
LDFLAGS += -L..

count_% : count_%.o ../libescape.a
	$(COMPILE_AND_LINK)
