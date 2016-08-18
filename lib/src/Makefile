CC=mpicc -std=gnu99 -Wall -Wno-unused-variable -O2 -D_FILE_OFFSET_BITS=64 -D_GNU_SOURCE

AR=ar

INCLUDE=../../include/*.h
INCLDDIR= -I../../include/ -I./

TARGET=../libqcd.a

all: $(TARGET)

SOURCES=\
	qcd_rng\
	qcd_blas\
	qcd_communication\
	qcd_gamma\
	qcd_gaugeFixing\
	qcd_init\
	qcd_io\
	qcd_observables\
	qcd_smearing\
	qcd_wilson\
	qcd_stochastic\
	qcd_local_time\
	qcd_smearing_special\
	qcd_blas_special\
	qcd_test\
#	threep_seq_all_J32\
#	qcd_all_particles\
#	qcd_all_particles_J32
#	twop_contract\
#	p1_p2_p1_fun_oc\
#	p1_p2_p1_fun_dec\
#	delta_fun\
#	lambdas_fun\
#	p1_p1_p1_fun\
#	p1_p2_p3_fun_dec\
#	p1_p2_p3_fun_oc\
#	p123_p321_fun\
#	sigmas2_fun\
#	sigmas4_fun\
#	qcd_twop_all_particles_J32\
#	qcd_all_particles\
#	qcd_all_particles_J32\
#	qcd_twop_all_particles_J32\
#	qcd_twop_all_particles


$(TARGET):		${addsuffix .o, $(SOURCES)} $(INCLUDE)
			$(AR) -r $@ ${addsuffix .o, $(SOURCES)}
			ranlib $@

%.o:			%.c $(INCLUDE)
			$(CC) -c $< $(INCLDDIR)
clean:		
		rm -f *.o libqcd.a ../libqcd.a
