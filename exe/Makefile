ESCAPE_HOME := ../

TARGETS := sample runQbS generateInput generateL0seed createCGraph generateL0COO queryOnCore

OBJECTS := $(TARGETS:%=%.o)

include ../common.mk

INCLUDES += -I..
LDFLAGS += -L..

count_% : count_%.o ../libescape.a
	$(COMPILE_AND_LINK)
