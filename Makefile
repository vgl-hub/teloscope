CXX ?= g++
INCLUDE_DIR ?= -I./include -Igfalibs/include
WARNINGS = -Wall -Wextra

CXXFLAGS = -g -std=gnu++17 -O3 $(INCLUDE_DIR) $(WARNINGS) $(CFLAGS)

TARGET = teloscope
TEST_TARGET = validate
GENERATE_TARGET = generate-tests
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

all: head validate regenerate

$(OBJS): %: $(BINDIR)/%
	@
$(BINDIR)%: $(SOURCE)/%.cpp $(INCLUDE)/%.h | $(BINDIR)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -c $(SOURCE)/$(notdir $@).cpp -o $@
	
.PHONY: gfalibs
gfalibs:
	$(MAKE) -j -C $(GFALIBS_DIR) CXXFLAGS="$(CXXFLAGS)"

validate: | $(BUILD)
	$(CXX) $(CXXFLAGS) -o $(BUILD)/$(TARGET)-$(TEST_TARGET) $(SOURCE)/$(TEST_TARGET).cpp $(LIBS)
	
regenerate: | $(BUILD)
	$(CXX) $(CXXFLAGS) -o $(BUILD)/$(TARGET)-$(GENERATE_TARGET) $(SOURCE)/$(GENERATE_TARGET).cpp $(LIBS)

	
$(BUILD):
	-mkdir -p $@
	
$(BINDIR):
	-mkdir -p $@
          
debug: CXXFLAGS += -DDEBUG -O0
debug: head

clean:
	$(RM) -r build
	$(MAKE) -C $(GFALIBS_DIR) clean
