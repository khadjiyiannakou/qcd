CC=mpicc

CFLAGS=-I./include -std=gnu99 -O3
LDFLAGS=-L./lib -L/usr/local/lib -L/home/tuj57552/install/lib/
LIBS=-lqcd -lgsl -llime -lgslcblas -lm

.PHONY: clean\
	cleanall\
	lib

#TARGETS=show_conf

TARGETS=b_minus_Dx\
	invert\
	landau\
	seq_source_dd\
	seq_source_dd_new_format\
	seq_source_uu_new_format\
	seq_source_dd_idris\
	seq_source_uu\
	seq_source_uu_idris\
	source\
	source_test\
	source_idris\
	threep_idris\
	twop\
	zfac_disc\
	sourceExp_new_format\
	zfac\
	zfac_disc_w_der\
	zfac_PDFs\
	gluonLoops
#	unit_gaugefield\
#	zfac\
#	show_conf\
#	check_prop\
#	twop_qkxTM\
#	check_seq_source\
#	seq_source_uu_new_format_001\
#	seq_source_dd_new_format_001\
#	seq_threep_fix_current_J32_correct\
#	seq_threep_fix_current_J32\
#	seq_source_fix_current\
#	seq_threep_fix_current\
#	threep_idris_wo_theta\
#	seq_threep_fix_current\
#	seq_threep_fix_current_J32\
#	twop_all_particles\
#	twop_all_particles_J32\
#	twop_all_particles_J32_dec\
#	twop_all_particles_J32_deltas\
#	twop_all_particles_J32_unproj\
#	stochastic_alloperators_allprojectors_up\
#	stochastic_alloperators_allprojectors_down\
#	stochastic_alloperators_allprojectors_up_nss\
#	stochastic_alloperators_allprojectors_down_nss\
#	stochastic_alloperators_allprojectors_up_parallel_time\
#	stochastic_alloperators_allprojectors_up_parallel_time_nss\
#	zfac_NEW_format\
#	twop_full\

#	twop_all_particles_J12 	twop_all_particles_J12_dec
all: lib ${addsuffix .exe, $(TARGETS)}

lib: lib/libqcd.a
lib/libqcd.a:
	cd lib/src &&\
	make

%.exe: %.o lib/libqcd.a
	$(CC) $< -o $@ $(LDFLAGS) $(LIBS)

%.o: %.c 
	$(CC) $(CFLAGS) -c $<

clean:		
	rm -vf ${addsuffix .o, $(TARGETS)}

cleanall: clean
	rm -vf ${addsuffix .exe, $(TARGETS)} &&\
	cd lib/src/ &&\
	make clean
