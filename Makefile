# If you want to use ifort, you need to specify it at 
# the command line. Usage: make compiler=ifort.

MAKE = make -i -r

ifeq ($(compiler),ifort)
        FC      =      ifort
        FFLAGS  =      -w -ffree-form -ffree-line-length-none -Wargument-mismatch
else 
        FC      =     gfortran 
#        FFLAGS  =      -w -ffree-form -ffree-line-length-none -Wargument-mismatch
        FFLAGS  =      -w -ffree-form -ffree-line-length-none
endif

OUT_DIR = out
RM      =      rm -rf
TARGET  =      kdm_vals.exe

kdm_kvals_code: 
	$(FC) $(FFLAGS) ck_generate.f90 -o $(TARGET)
	$(shell if [ ! -d "./$(OUT_DIR)" ];then mkdir $(OUT_DIR); fi;)

clean:
	$(RM) $(TARGET)

ultraclean:
	$(RM) $(TARGET)
	$(RM) out

