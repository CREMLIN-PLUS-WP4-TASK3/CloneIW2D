CC:=g++
CCFLAGS:=-O2 -fPIC
CCFLAGS+=-Wall -Wextra -pedantic
OUTDIR:=build
BINDIR:=bin
SRCDIR:=.
INCLUDE=-I$(SRCDIR)
INCLUDE+=-I../External_libs/ALGLIB/include
LCFLAGS:=-L../External_libs/ALGLIB/
LCFLAGS+=-lalg -lmpfr -lgmp -lm -lgsl
BINNAMES:=roundchamber wake_roundchamber flatchamber wake_flatchamber
LIBNAME:=libIW2D

define compilebin=
	$(CC) -static $(CCFLAGS) $(INCLUDE) $^ -Wl,-Bstatic $(LCFLAGS) -o $@
endef

static: $(OUTDIR) $(BINDIR) $(BINDIR)/round $(BINDIR)/wakeround $(BINDIR)/flat \
	$(BINDIR)/wakeflat $(BINDIR)/$(LIBNAME).so

$(OUTDIR)/%.o: $(SRCDIR)/%.cc
	$(CC) $(CCFLAGS) -c $(INCLUDE) $^ -o $@

$(OUTDIR):
	mkdir -p $(OUTDIR)

$(BINDIR):
	mkdir -p $(BINDIR)

$(BINDIR)/round: $(OUTDIR)/roundchamber.o
	$(compilebin)

$(BINDIR)/wakeround: $(OUTDIR)/wake_roundchamber.o
	$(compilebin)

$(BINDIR)/flat: $(OUTDIR)/flatchamber.o
	$(compilebin)

$(BINDIR)/wakeflat: $(OUTDIR)/wake_flatchamber.o
	$(compilebin)

$(BINDIR)/$(LIBNAME).so: $(OUTDIR)/wake_roundchamber.o
	$(CC) -shared $(CCFLAGS) $(INCLUDE) $^ -Wl,-Bdynamic $(LCFLAGS) -o $@

.PHONY: clean
clean:
	rm -rf $(OUTDIR) $(BINDIR)
