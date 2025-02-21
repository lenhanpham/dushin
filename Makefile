# Compiler and flags
FC = ifort
FFLAGS_DEBUG = -g -mcmodel=medium -m64
FFLAGS_OPT = -O2 -mcmodel=medium -m64
LDFLAGS = -g


# Directories
OBJDIR = obj
BINDIR = bin

# Common object files
COMMON_OBJS = $(OBJDIR)/subs.o \
              $(OBJDIR)/recalc-freq.o \
              $(OBJDIR)/proj0freq.o \
              $(OBJDIR)/bmatred.o \
              $(OBJDIR)/ddiag.o \
              $(OBJDIR)/dmpower.o \
              $(OBJDIR)/dmatinv.o

# Executables
EXECUTABLES = $(BINDIR)/dushin \
              $(BINDIR)/displace \
              $(BINDIR)/compare-geom

# Create directories if they don't exist
$(shell mkdir -p $(OBJDIR) $(BINDIR))

# Default target
all: $(EXECUTABLES)

# Debug objects (compiled with -g)
$(OBJDIR)/%.o: %.f
	$(FC) $(FFLAGS_DEBUG) $(FCFLAGS) -c $< -o $@

# Optimized objects (compiled with -O2)
$(OBJDIR)/recalc-freq.o: recalc-freq.f
	$(FC) $(FFLAGS_OPT) $(FCFLAGS) -c $< -o $@

$(OBJDIR)/proj0freq.o: proj0freq.f
	$(FC) $(FFLAGS_OPT) $(FCFLAGS) -c $< -o $@

$(OBJDIR)/bmatred.o: bmatred.f
	$(FC) $(FFLAGS_OPT) $(FCFLAGS) -c $< -o $@

$(OBJDIR)/ddiag.o: ddiag.f
	$(FC) $(FFLAGS_OPT) $(FCFLAGS) -c $< -o $@

$(OBJDIR)/dmpower.o: dmpower.f
	$(FC) $(FFLAGS_OPT) $(FCFLAGS) -c $< -o $@

$(OBJDIR)/dmatinv.o: dmatinv.f
	$(FC) $(FFLAGS_OPT) $(FCFLAGS) -c $< -o $@

# Keep subs-dftb.o compilation rule but don't include it in COMMON_OBJS
$(OBJDIR)/subs-dftb.o: subs-dftb.f
	$(FC) $(FFLAGS_DEBUG) $(FCFLAGS) -c $< -o $@

# Executables
$(BINDIR)/dushin: $(OBJDIR)/dushin.o $(COMMON_OBJS)
	$(FC) $(LDFLAGS) $(FCFLAGS) $^ -o $@

$(BINDIR)/displace: $(OBJDIR)/displace.o $(COMMON_OBJS)
	$(FC) $(LDFLAGS) $(FCFLAGS) $^ -o $@

$(BINDIR)/compare-geom: $(OBJDIR)/compare-geom.o $(COMMON_OBJS)
	$(FC) $(LDFLAGS) $(FCFLAGS) $^ -o $@

# Phony targets
.PHONY: all clean help

clean:
	rm -rf $(OBJDIR) $(BINDIR)

help:
	@echo "Available targets:"
	@echo "  all          - Build all executables (default)"
	@echo "  clean        - Remove all built files"
	@echo "  help         - Show this help message"
	@echo ""
	@echo "Executables will be built in $(BINDIR)/"
	@echo "Object files will be stored in $(OBJDIR)/"



