ESCAPE_HOME := ../

TARGETS := sample generateL0seed generateInput createCGraph sanitize sampleL0 printL2 createL0CGraph PLL

OBJECTS := $(TARGETS:%=%.o)

include ../common.mk

INCLUDES += -I..
LDFLAGS += -L..

count_% : count_%.o ../libescape.a
	$(COMPILE_AND_LINK)
