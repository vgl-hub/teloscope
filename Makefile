CXX ?= g++
INCLUDE_DIR ?= -I./include -Igfalibs/include
WARNINGS = -Wall -Wextra

CXXFLAGS = -g -std=gnu++14 -O3 $(INCLUDE_DIR) $(WARNINGS) $(CFLAGS)

TARGET = teloscope
BUILD = build/bin
SOURCE = src
INCLUDE = include
BINDIR := $(BUILD)/.o
LIBS = -lz
LDFLAGS = -pthread

GFALIBS_DIR := $(CURDIR)/gfalibs

OBJS := main teloscope input
BINS := $(addprefix $(BINDIR)/, $(OBJS))

head: $(BINS) gfalibs | $(BUILD)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $(BUILD)/$(TARGET) $(wildcard $(BINDIR)/*) $(GFALIBS_DIR)/*.o $(LIBS)
	
debug: CXXFLAGS += -DDEBUG
debug: CCFLAGS += -DDEBUG
debug: head

all: head 

$(OBJS): %: $(BINDIR)/%
	@
$(BINDIR)%: $(SOURCE)/%.cpp $(INCLUDE)/%.h | $(BINDIR)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -c $(SOURCE)/$(notdir $@).cpp -o $@
	
.PHONY: gfalibs
gfalibs:
	$(MAKE) -j -C $(GFALIBS_DIR) CXXFLAGS="$(CXXFLAGS)"
	
$(BUILD):
	-mkdir -p $@
	
$(BINDIR):
	-mkdir -p $@
          
debug: CXXFLAGS += -DDEBUG -O0
debug: head

clean:
	$(RM) -r build
	$(MAKE) -C $(GFALIBS_DIR) clean
