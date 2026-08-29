FC = gfortran
FFLAGS ?= -O2 -std=legacy -ffree-line-length-none -fno-backtrace
LDFLAGS ?=

ifeq ($(OS),Windows_NT)
  EXEEXT := .exe
  PYTHON ?= python
else
  EXEEXT :=
  PYTHON ?= python3
endif

# Pick deletion/run commands for the shell make actually uses for recipes:
# cmd.exe echoes the quotes around "x", POSIX sh strips them. A cmd /c wrapper
# must not be used here - under an sh-driven make (MSYS2/Git Bash) MSYS path
# conversion mangles /c, making the command a silent no-op with exit 0.
ifeq ($(shell echo "x"),"x")
  RM_FILES = del /q $(subst /,\,$(1))
  TEST_RUNNER := tests\test_integration$(EXEEXT)
else
  RM_FILES = rm -f $(1)
  TEST_RUNNER := ./tests/test_integration$(EXEEXT)
endif

TARGET := electron_fluxes$(EXEEXT)
OBJS := subroutines.o electron_fluxes.o
TEST_TARGET := tests/test_integration$(EXEEXT)

.PHONY: all clean debug test run-grid diagnostics

all: $(TARGET)

$(TARGET): $(OBJS)
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJS)

electron_fluxes.o: electron_fluxes.f90 subroutines.o
	$(FC) $(FFLAGS) -c $< -o $@

subroutines.o: subroutines.f90
	$(FC) $(FFLAGS) -c $< -o $@

$(TEST_TARGET): tests/test_integration.f90 subroutines.o
	$(FC) $(FFLAGS) -o $@ tests/test_integration.f90 subroutines.o

test: $(TARGET) $(TEST_TARGET)
	$(TEST_RUNNER)
	$(PYTHON) -m unittest discover -s tests -p "test_*.py" -v

# Both targets drop their objects afterwards so a later plain `make`/`make test`
# rebuilds with the default FFLAGS instead of silently reusing -O0 objects.
debug:
	$(MAKE) clean
	$(MAKE) FFLAGS="-O0 -g -std=legacy -ffree-line-length-none -fcheck=all -fbacktrace -finit-real=snan -finit-integer=-999999" $(TARGET)
	-$(call RM_FILES,$(OBJS) *.mod)

diagnostics:
	$(MAKE) clean
	$(MAKE) FFLAGS="-O0 -g -std=legacy -ffree-line-length-none -Wall -Wextra -Wimplicit-interface -Wconversion-extra" $(TARGET)
	-$(call RM_FILES,$(OBJS) *.mod)

run-grid: $(TARGET)
	$(PYTHON) run_on_grid.py

clean:
	-$(call RM_FILES,$(TARGET) $(OBJS) *.mod tests/test_integration$(EXEEXT))
