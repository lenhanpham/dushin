# Compiler and flags
FC = gfortran
FFLAGS_DEBUG = -g -fdefault-real-8 -fdefault-double-8
FFLAGS_OPT = -O2 -fdefault-real-8 -fdefault-double-8
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
$(OBJDIR)/%.o: %.for
	$(FC) $(FFLAGS_DEBUG) -c $< -o $@

# Optimized objects (compiled with -O2)
$(OBJDIR)/recalc-freq.o: recalc-freq.for
	$(FC) $(FFLAGS_OPT) -c $< -o $@

$(OBJDIR)/proj0freq.o: proj0freq.for
	$(FC) $(FFLAGS_OPT) -c $< -o $@

$(OBJDIR)/bmatred.o: bmatred.for
	$(FC) $(FFLAGS_OPT) -c $< -o $@

$(OBJDIR)/ddiag.o: ddiag.for
	$(FC) $(FFLAGS_OPT) -c $< -o $@

$(OBJDIR)/dmpower.o: dmpower.for
	$(FC) $(FFLAGS_OPT) -c $< -o $@

$(OBJDIR)/dmatinv.o: dmatinv.for
	$(FC) $(FFLAGS_OPT) -c $< -o $@

# Keep subs-dftb.o compilation rule but don't include it in COMMON_OBJS
$(OBJDIR)/subs-dftb.o: subs-dftb.for
	$(FC) $(FFLAGS_DEBUG) -c $< -o $@

# Executables
$(BINDIR)/dushin: $(OBJDIR)/dushin.o $(COMMON_OBJS)
	$(FC) $(LDFLAGS) $^ -o $@

$(BINDIR)/displace: $(OBJDIR)/displace.o $(COMMON_OBJS)
	$(FC) $(LDFLAGS) $^ -o $@

$(BINDIR)/compare-geom: $(OBJDIR)/compare-geom.o $(COMMON_OBJS)
	$(FC) $(LDFLAGS) $^ -o $@

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



