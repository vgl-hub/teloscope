CXX ?= g++
INCLUDE_DIR ?= -I./include -Igfalibs/include
WARNINGS = -Wall -Wextra

TELOSCOPE_COMMIT := $(shell git rev-parse --short HEAD 2>/dev/null || printf '%s' unknown)

CXXFLAGS = -g -std=gnu++17 -O3 $(INCLUDE_DIR) $(WARNINGS) $(CFLAGS) -DTELOSCOPE_COMMIT=\"$(TELOSCOPE_COMMIT)\"
DEPFLAGS = -MMD -MP

TARGET = teloscope
TEST_TARGET = validate
GENERATE_TARGET = generate-tests
BUILD = build/bin
SOURCE = src
INCLUDE = include
BINDIR := $(BUILD)/.o
LIBS = -lz
LDFLAGS = -pthread

# Static linking on Windows to avoid DLL dependencies
ifeq ($(OS),Windows_NT)
    LDFLAGS += -static
endif

GFALIBS_DIR := $(CURDIR)/gfalibs

OBJS := main teloscope input tools read-filter bgzf bam
BINS := $(addprefix $(BINDIR)/, $(OBJS))
DEPFILES := $(addsuffix .d, $(BINS))
COMMIT_STATE := $(BINDIR)/.teloscope-commit

head: $(BINS) gfalibs | $(BUILD)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $(BUILD)/$(TARGET) $(BINS) $(GFALIBS_DIR)/*.o $(LIBS)
	
debug: CXXFLAGS += -DDEBUG -O0
debug: head

all: head validate regenerate simulate

$(OBJS): %: $(BINDIR)/%
	@
$(BINDIR)/teloscope: $(COMMIT_STATE)
$(BINDIR)%: $(SOURCE)/%.cpp $(INCLUDE)/%.h | $(BINDIR)
	$(CXX) $(CXXFLAGS) $(DEPFLAGS) $(LDFLAGS) -c $(SOURCE)/$(notdir $@).cpp -o $@

.PHONY: FORCE
FORCE:

# CXXFLAGS are not make prerequisites. Keep a content-sensitive stamp so the
# object that embeds the commit is rebuilt after HEAD changes, but not on every
# incremental build.
$(COMMIT_STATE): FORCE | $(BINDIR)
	@current="$$(test -f "$@" && sed -n '1p' "$@")"; \
	if test "$$current" != "$(TELOSCOPE_COMMIT)"; then \
		printf '%s\n' "$(TELOSCOPE_COMMIT)" > "$@"; \
	fi
	
.PHONY: gfalibs
gfalibs:
	$(MAKE) -j -C $(GFALIBS_DIR) CXXFLAGS="$(CXXFLAGS)"

validate: | $(BUILD)
	$(CXX) $(CXXFLAGS) -o $(BUILD)/$(TARGET)-$(TEST_TARGET) $(SOURCE)/$(TEST_TARGET).cpp $(LIBS)
	
regenerate: | $(BUILD)
	$(CXX) $(CXXFLAGS) -o $(BUILD)/$(TARGET)-$(GENERATE_TARGET) $(SOURCE)/$(GENERATE_TARGET).cpp $(LIBS)

simulate: | $(BUILD)
	$(CXX) $(CXXFLAGS) -o $(BUILD)/$(TARGET)-simulate $(SOURCE)/get-mock-chr.cpp

$(BUILD):
	-mkdir -p $@
	
$(BINDIR):
	-mkdir -p $@

test-gaps: head
	bash scripts/test_gaps_bed.sh

test-n50: head
	bash scripts/test_n50.sh

test-bam: head
	python3 scripts/test_bam_subset.py

.PHONY: test-filters
test-filters: head
	TELOSCOPE="$(BUILD)/$(TARGET)" python3 scripts/test_sequence_filters.py

.PHONY: test-bam-hardening test-bam-coverage test-bam-sanitize
test-bam-hardening: head
	TELOSCOPE="$(BUILD)/$(TARGET)" BAM_MUTATION_CASES="$${BAM_MUTATION_CASES:-512}" python3 scripts/test_bam_subset.py
	BUILD_DIR="$(BUILD)" CXX="$(CXX)" bash scripts/test_bgzf_faults.sh
	TELOSCOPE="$(BUILD)/$(TARGET)" python3 scripts/test_bam_samtools.py

test-bam-coverage:
	CXX="$(CXX)" bash scripts/test_bam_coverage.sh

test-bam-sanitize:
	CXX="$(CXX)" bash scripts/test_bam_sanitized.sh

gfa-oracle: head
	bash scripts/compare_to_reference.sh

clean:
	$(RM) -r build
	$(MAKE) -C $(GFALIBS_DIR) clean

-include $(DEPFILES)
