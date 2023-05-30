# RHEL7 
CERN_ROOT = /apps/cernlib/x86_64_rhel7/2005/
OTHERLIBS = -L$(CERN_ROOT)/lib -lmathlib
FFLAGS    = -C -g -w -fno-automatic -fbounds-check -ffixed-line-length-132
F77       :=gfortran
########################################

externals_all_objs = externals_all.o sigmodel_calc_hybrid.o fy.o target_info.o \
	             F1F209.o f_to_sig.o y_calc.o sig_bar_df.o nform_xem.o \
                     get_cc_info.o bdisnew4he3_ineft.o
externals_all_srcs = externals_all.f sigmodel_calc_hybrid.f fy.f target_info.f \
                     F1F209.f f_to_sig.f y_calc.f sig_bar_df.f nform_xem.f \
                     get_cc_info.f bdisnew4he3_ineft.f

#externals_all_objs = externals_all.o sigmodel_calc.o fy.o target_info.o \
	             f_to_sig.o y_calc.o sig_bar_df.o nform_xem.o \
                     get_cc_info.o bdisnew4he3.o f1f2in06.o
#externals_all_srcs = externals_all.f sigmodel_calc.f fy.f target_info.f \
                     f_to_sig.f y_calc.f sig_bar_df.f nform_xem.f \
                     get_cc_info.f bdisnew4he3.f f1f2in06.f

########################################

none: externals_all

all:  externals_all

##externals_all:
##	$(F77) $(FFLAGS) -o jtest $(externals_all_srcs) $(OTHERLIBS)

externals_all.o: externals_all.f
		 $(F77) $(FFLAGS) -c $< -o $@

fy.o: fy.f
	$(F77) $(FFLAGS) -c $< -o $@

sigmodel_calc_simple.o: sigmodel_calc_simple.f
	$(F77) $(FFLAGS) -c $< -o $@

sigmodel_calc.o: sigmodel_calc.f
	$(F77) $(FFLAGS) -c $< -o $@

sigmodel_calc_hybrid.o: sigmodel_calc_hybrid.f
	$(F77) $(FFLAGS) -c $< -o $@

bdisnew4he3_ineft.o: bdisnew4he3_ineft.f
	$(F77) $(FFLAGS) -c $< -o $@

target_info.o: target_info.f
	$(F77) $(FFLAGS) -c $< -o $@

f1f2in06.o: f1f2in06.f
	$(F77) $(FFLAGS) -c $< -o $@

F1F209.o: F1F209.f
	$(F77) $(FFLAGS) -c $< -o $@

#f1f209newmecnof1cor.o: f1f209newmecnof1cor.f
#	$(F77) $(FFLAGS) -c $< -o $@

f_to_sig.o: f_to_sig.f
	$(F77) $(FFLAGS) -c $< -o $@

y_calc.o: y_calc.f
	$(F77) $(FFLAGS) -c $< -o $@

sig_bar_df.o: sig_bar_df.f
	$(F77) $(FFLAGS) -c $< -o $@

nform_xem.o: nform_xem.f
	$(F77) $(FFLAGS) -c $< -o $@

get_cc_info.o: get_cc_info.f
	$(F77) $(FFLAGS) -c $< -o $@

#f2allm.o: f2allm.f
#	$(F77) $(FFLAGS) -c $< -o $@
#r1998.o: r1998.f
#	$(F77) $(FFLAGS) -c $< -o $@

externals_all: $(externals_all_objs) Makefile
	$(F77) -o $@  $(FFLAGS) $(externals_all_objs) $(OTHERLIBS)


clean:
	rm -f *.o externals_all
