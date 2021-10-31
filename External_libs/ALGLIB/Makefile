CC:=g++
CCFLAGS:=-fPIC -O2 -Wall -Wextra
SRCDIR:=src
INCLUDE:=-Iinclude
LIBNAME:=libalg
OUTDIR:=builddir
SRC=$(notdir $(wildcard $(SRCDIR)/*.cpp))
TARGET=$(patsubst %.cpp,%.o,$(SRC))

$(LIBNAME): $(OUTDIR) compile

.PHONY: compile
compile: $(addprefix $(OUTDIR)/,$(TARGET))
	ar -rcs $(LIBNAME).a $^
	$(CC) $(CCFLAGS) -shared $(INCLUDE) -o $(LIBNAME).so $^

$(OUTDIR)/%.o: $(SRCDIR)/%.cpp
	$(CC) $(CCFLAGS) -c -o $@ $(INCLUDE) $^

$(OUTDIR):
	mkdir -p $(OUTDIR)

.PHONY: clean
clean:
	rm -rf $(OUTDIR) $(LIBNAME).so $(LIBNAME).a
