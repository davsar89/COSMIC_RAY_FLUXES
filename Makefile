FC = gfortran
FFLAGS ?= -O2 -std=legacy -ffree-line-length-none -fno-backtrace
LDFLAGS ?=

ifeq ($(OS),Windows_NT)
  EXEEXT := .exe
  PYTHON ?= python
  RM := cmd /c "del /q"
  TEST_RUNNER := cmd /c "tests\test_integration$(EXEEXT)"
else
  EXEEXT :=
  PYTHON ?= python3
  RM := rm -f
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

debug:
	$(MAKE) clean
	$(MAKE) FFLAGS="-O0 -g -std=legacy -ffree-line-length-none -fcheck=all -fbacktrace -finit-real=snan -finit-integer=-999999" $(TARGET)

diagnostics:
	$(MAKE) clean
	$(MAKE) FFLAGS="-O0 -g -std=legacy -ffree-line-length-none -Wall -Wextra -Wimplicit-interface -Wconversion-extra" $(TARGET)

run-grid: $(TARGET)
	$(PYTHON) run_on_grid.py

clean:
	-$(RM) $(TARGET) $(OBJS) *.mod
ifeq ($(OS),Windows_NT)
	-cmd /c "del /q tests\test_integration$(EXEEXT)"
else
	-$(RM) $(TEST_TARGET)
endif
