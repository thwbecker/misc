#CFLAGS = -g
#FFLAGS = -g
#F90FLAGS = -g

LMFIT_DIR = ../lmfit-5.1
LMFIT_LIBS = -L$(LMFIT_DIR)/lib/ -llmfit -static
LMFIT_INC = -I$(LMFIT_DIR)/include/


BDIR = $(HOME)/progs/bin
ODIR = objects
#
# GMT4 stuff
#
GMT = $(GMT4HOME)/
GMT_LIBS = -L$(GMT)/lib/ -L$(NETCDFDIR)/lib/ -lpsl -lgmt -lnetcdf

GMT_HDRS = $(GMT4HOME)/include/gmt.h
GMT_INC = -I$(GMT4HOME)/include/ -I$(NETCDFDIR)/include/ -I$(GMT4HOME)/src/meca/ -DUSE_GMT4

GGRD_INC = -I../cig/hc/ -I$(GMT4HOME)/include/ -I$(NETCDFDIR)/include/
GGRD_LIB = -L../cig/hc/objects/$(ARCH)/ -L$(GMT4HOME)/lib/  -lggrd -lhc  -lgmt \
	-L$(NETCDFDIR)/lib/  -lnetcdf 
#
# GMT6
#
GMT6_INC = -I$(GMTHOME)/include/gmt/ -I/usr/include/gdal/
GMT6_LIB = -L$(GMTHOME)lib/ -lgmt -lm

EISLIB = -L$(HOME)/progs/src/fstrack/eispack/$(ARCH)/ -lmyeis

INCLUDES = 

MT_LIBS = ../mtspec/mtspec/src/libraries/mwlib.a ../lpmtm/libjtap.a



PROGS =		$(BDIR)/grd2poisson 	$(BDIR)/grds2profile	\
	$(BDIR)/grd2poisson_my 	$(BDIR)/grd2cumuhist $(BDIR)/grd2poisson_mud  

PROGS1 =  $(BDIR)/prune_gps_data \
	$(BDIR)/erfc 	$(BDIR)/half_space_T $(BDIR)/erf $(BDIR)/make_two_layer_split_model \
	$(BDIR)/eval_two_layer_split_model $(BDIR)/read_itl_sol \
	$(BDIR)/invert_two_layer_split $(BDIR)/invert_two_layer_split_field \
	$(BDIR)/asciiint2bebin $(BDIR)/asciifloat2scaled_int \
	$(BDIR)/asciifloat2bin $(BDIR)/sumfloatbin \
	$(BDIR)/aki2cmt	\
	$(BDIR)/asciifloat2bebin $(BDIR)/ps2raster $(BDIR)/spline_sample	\
	$(BDIR)/asciiint2bin $(BDIR)/asciiuint2bin $(BDIR)/fit_tgr $(BDIR)/generate_tgr \
	$(BDIR)/dbl2ascii  $(BDIR)/assign_neighbor_code  $(BDIR)/combineoff \
	$(BDIR)/matrix2vector   $(BDIR)/lonlatpv2cv  $(BDIR)/geomean \
	$(BDIR)/misfit  $(BDIR)/random_gauss  $(BDIR)/random \
	$(BDIR)/random_lonlat  $(BDIR)/random_tensor $(BDIR)/random_power \
	$(BDIR)/random_exp  $(BDIR)/extract_dens  $(BDIR)/random_tgr \
	$(BDIR)/rtp2xyz $(BDIR)/cov_from_events $(BDIR)/random_gvel $(BDIR)/xyz2lonlatz $(BDIR)/randomize_simple_xyz \
	$(BDIR)/simple_spline $(BDIR)/extract_bin  $(BDIR)/extract_bin_file  $(BDIR)/extract_bin_column \
	 $(BDIR)/extract_bin_int $(BDIR)/spat_power_grid $(BDIR)/close2point \
	$(BDIR)/close2pointz $(BDIR)/close2points \
	$(BDIR)/determine_coeff  $(BDIR)/resample  $(BDIR)/extract_bin_int_file 
# $(BDIR)/make_pdf 	$(BDIR)/fit_exp_scan  
#	$(BDIR)/spherical_spline_invert \
#	$(BDIR)/spherical_spline_evaluate $(BDIR)/spherical_spline_construct 

PROGS2 = \
	$(BDIR)/crosscorrelate  $(BDIR)/crosscorrelate_window \
	$(BDIR)/crosscorrelate_delay  $(BDIR)/stddev  $(BDIR)/gtr1_code  $(BDIR)/gtr1_avg \
	$(BDIR)/ran_angle $(BDIR)/gutenberg $(BDIR)/eigen  $(BDIR)/eigenvalues \
	$(BDIR)/eigen3ds  $(BDIR)/eigenvalues3ds $(BDIR)/calcms  \
	$(BDIR)/invert3x3mat $(BDIR)/gpoints $(BDIR)/correl \
	$(BDIR)/rotatesmat $(BDIR)/rotatecmat \
	$(BDIR)/rotatesvec $(BDIR)/rotatecvec $(BDIR)/fft $(BDIR)/comb2col $(BDIR)/fitxy \
	$(BDIR)/fitxyee $(BDIR)/p_fisher $(BDIR)/p_student \
	$(BDIR)/pearson_corr $(BDIR)/bootstrap_corr $(BDIR)/bootstrap_spear \
	$(BDIR)/spear $(BDIR)/period $(BDIR)/student  \
	$(BDIR)/fit_plane $(BDIR)/averageirrdata_spline \
	$(BDIR)/averageirrdata_linear $(BDIR)/uudecode \
	$(BDIR)/averageirrdata_quad \
	$(BDIR)/averageirrdata_cubic $(BDIR)/banner \
	$(BDIR)/lonlat2gcpoints $(BDIR)/lonlatzcv2pv \
	$(BDIR)/gifmerge $(BDIR)/median	$(BDIR)/svdfit $(BDIR)/momt_project \
	$(BDIR)/pca $(BDIR)/polyfit $(BDIR)/plume_detect  $(BDIR)/plume_detect_d \
	$(BDIR)/paste_bin $(BDIR)/dist_on_sphere_and_azi \
	$(BDIR)/cholesky_lsq $(BDIR)/svd_lsq $(BDIR)/svd_lsq_s $(BDIR)/ftest  
#	$(BDIR)/solve_hom3d 
#	$(BDIR)/iaspei91 $(BDIR)/mpi_example

#$(BDIR)/poly_avg_splitting \
#$(BDIR)/pointingmtpoly \
#
# interact geoemtry operations
#
INTACTGEOOBJS =	\
	$(HOME)/progs/src/interact/objects/geometry.o \
	$(HOME)/progs/src/interact/objects/llgeo.o \
	$(HOME)/progs/src/interact/objects/mysincos.o \
	$(HOME)/progs/src/interact/objects/fit_plane.o

DEFS = -DCOMP_PRECISION=double

HDR_FLS = rand.h rotate.h trig.h polysort.h coordconvention.h

all: 
	make $(PROGS);\
	make $(PROGS1);\
	make $(PROGS2)


tmp:
	make	$(BDIR)/plume_detect  $(BDIR)/plume_detect_d 


really_all: all $(BDIR)/remove_net_rotation $(BDIR)/determine_net_rotation   \
	$(BDIR)/unhqx

eigen_tools: $(BDIR)/eigen  $(BDIR)/eigenvalues  \
	$(BDIR)/eigen3ds  $(BDIR)/eigenvalues3ds

clean:
	rm $(ODIR)/*.o ; rm $(PROGS1) $(PROGS2)

cproto:
	cd cproto-4.6;\
	configure --prefix $(HOME)/progs/$(ARCH)/;\
	make ; make install;

$(BDIR)/fit_plane: $(ODIR)/fit_plane.o $(INTACTGEOOBJS) \
	$(HOME)/progs/src/interact/interact.h \
	$(HOME)/progs/src/interact/auto_proto.h 
	$(CC) $(CFLAGS) $(ODIR)/fit_plane.o	\
		$(INTACTGEOOBJS) -o $(BDIR)/fit_plane \
		$(MATHLIBS) $(LDFLAGS) 

$(ODIR)/libnrtools.a:  $(ODIR)/nr_util.o $(ODIR)/svd_util.o $(ODIR)/datafit_util.o
	$(AR) rv $(ODIR)/libnrtools.a  \
	$(ODIR)/nr_util.o $(ODIR)/svd_util.o \
	$(ODIR)/datafit_util.o

$(BDIR)/gutenberg: gutenberg.c $(ODIR)/momentmag.o
	$(CC) $(CFLAGS) gutenberg.c  $(ODIR)/momentmag.o -o $(BDIR)/gutenberg $(LDFLAGS)

$(BDIR)/make_pdf: make_pdf.c 
	$(CC) $(CFLAGS) make_pdf.c  -o $(BDIR)/make_pdf $(LDFLAGS) -llapack

$(BDIR)/close2point: close2point.c
	$(CC) $(CFLAGS) close2point.c -o $(BDIR)/close2point $(LDFLAGS)

$(BDIR)/close2pointz: close2pointz.c
	$(CC) $(CFLAGS) close2pointz.c -o $(BDIR)/close2pointz $(LDFLAGS)


$(BDIR)/close2points: close2points.c
	$(CC) $(CFLAGS) close2points.c -o $(BDIR)/close2points $(LDFLAGS)

$(BDIR)/sumfloatbin: sumfloatbin.c
	$(CC) $(CFLAGS) sumfloatbin.c -o $(BDIR)/sumfloatbin $(LDFLAGS)

$(BDIR)/plume_detect: $(ODIR)/plume_detect.o
	$(CC) $(CFLAGS)  $(GGRD_INC) $(GMT_INC) \
	$(ODIR)/plume_detect.o -o $(BDIR)/plume_detect \
	$(GGRD_LIB) $(LDFLAGS) 

$(BDIR)/grds2profile: grds2profile.c
	$(CC) $(CFLAGS)  $(GGRD_INC) $(GMT_INC) grds2profile.c -o $(BDIR)/grds2profile \
		$(GGRD_LIB) $(LDFLAGS)

$(BDIR)/plume_detect_d: $(ODIR)/plume_detect_d.o
	$(CC) $(CFLAGS)  $(GGRD_INC) $(GMT_INC) \
	$(ODIR)/plume_detect_d.o -o $(BDIR)/plume_detect_d \
	$(GGRD_LIB) $(LDFLAGS) 

$(BDIR)/gtr1_code: $(ODIR)/gtr1_code.o $(ODIR)/evgtr1.o
	$(CC) $(CFLAGS) $(ODIR)/gtr1_code.o $(ODIR)/evgtr1.o -o $(BDIR)/gtr1_code $(FTRN_LIB) $(LDFLAGS) 

$(BDIR)/gtr1_avg: $(ODIR)/gtr1_avg.o $(ODIR)/evgtr1.o
	$(CC) $(CFLAGS) $(ODIR)/gtr1_avg.o $(ODIR)/evgtr1.o -o $(BDIR)/gtr1_avg $(FTRN_LIB) $(LDFLAGS)

$(BDIR)/polyfit: polyfit.c $(ODIR)/libnrtools.a \
	$(ODIR)/misc.o 
	$(CC) $(CFLAGS) polyfit.c \
	-o $(BDIR)/polyfit -L$(ODIR)/ $(ODIR)/misc.o \
	-lnrtools $(LDFLAGS)

$(BDIR)/uudecode: uudecode.c
	$(CC) $(CFLAGS) uudecode.c -o $(BDIR)/uudecode $(LDFLAGS)

$(BDIR)/randomize_simple_xyz: randomize_simple_xyz.c
	$(CC) $(CFLAGS) randomize_simple_xyz.c -o $(BDIR)/randomize_simple_xyz $(LDFLAGS)

$(BDIR)/cholesky_lsq: $(ODIR)/cholesky_lsq.o
	$(CC) $(CFLAGS) $(ODIR)/cholesky_lsq.o -o $(BDIR)/cholesky_lsq $(MATHLIBS) $(LDFLAGS)

$(BDIR)/svd_lsq: svd_lsq.c
	$(CC) $(CFLAGS) svd_lsq.c -DOUBLE_PRECISION -o $(BDIR)/svd_lsq $(MATHLIBS) $(LDFLAGS)

$(BDIR)/svd_lsq_s: svd_lsq.c
	$(CC) $(CFLAGS) svd_lsq.c -o $(BDIR)/svd_lsq_s $(MATHLIBS) $(LDFLAGS)

$(BDIR)/erfc: erfc.c
	$(CC) $(CFLAGS) erfc.c -o $(BDIR)/erfc $(LDFLAGS)
$(BDIR)/erf: erf.c
	$(CC) $(CFLAGS) erf.c -o $(BDIR)/erf $(LDFLAGS)

$(BDIR)/half_space_T: half_space_T.c
	$(CC) $(CFLAGS) half_space_T.c -o $(BDIR)/half_space_T $(LDFLAGS)

$(BDIR)/ftest: ftest.c
	$(CC) $(CFLAGS) ftest.c -o $(BDIR)/ftest $(LDFLAGS)

$(BDIR)/ps2raster: ps2raster.c
	$(CC) $(CFLAGS) ps2raster.c -o $(BDIR)/ps2raster $(LDFLAGS)

$(BDIR)/pointingmtpoly: pointingmtpoly.c $(GMT_HDRS) $(ODIR)/polysort.o 
	$(CC) $(CFLAGS) pointingmtpoly.c $(GMT_INC)  $(ODIR)/polysort.o   \
	-o $(BDIR)/pointingmtpoly $(GMT_LIBS) $(LDFLAGS)

$(BDIR)/poly_avg_splitting: poly_avg_splitting.c $(GMT_HDRS)  $(ODIR)/polysort.o 
	$(CC) $(CFLAGS) poly_avg_splitting.c $(GMT_INC)  \
	-o $(BDIR)/poly_avg_splitting   $(ODIR)/polysort.o $(GMT_LIBS) $(LDFLAGS)
#
# GMT6
#
$(BDIR)/grd2cumuhist: grd2cumuhist.c 
	$(CC) $(CFLAGS)  grd2cumuhist.c $(GMT6_INC)  -o $(BDIR)/grd2cumuhist $(GMT6_LIB)

$(ODIR)/grd2poisson.o: grd2poisson.c 
	$(CC) $(CFLAGS) -c grd2poisson.c $(GMT6_INC) -I../fishpak/ \
	-o $(ODIR)/grd2poisson.o   

#
# GMT 4
#

$(ODIR)/grd2poisson_mud.o: grd2poisson_mud.c  ../mudpack/mudpack.h
	$(CC) $(CFLAGS) -c grd2poisson_mud.c $(GMT_INC) -I../mudpack/ \
	-o $(ODIR)/grd2poisson_mud.o


$(ODIR)/grd2poisson_my.o: grd2poisson_my.c   
	$(CC) $(CFLAGS) -c grd2poisson_my.c $(GMT_INC) -I. -o $(ODIR)/grd2poisson_my.o   

# this didn't work?!
$(BDIR)/pointinpoly: pointinpoly.c $(ODIR)/polysort.o  polysort.h
	$(CC) $(CFLAGS) pointinpoly.c $(GMT_INC)  $(ODIR)/polysort.o \
	-o $(BDIR)/pointinpoly $(GMT_LIBS)  $(LDFLAGS)



$(BDIR)/make_two_layer_split_model: make_two_layer_split_model.c  $(ODIR)/calc_two_layer_split.o 
	$(CC) $(CFLAGS) make_two_layer_split_model.c  $(ODIR)/calc_two_layer_split.o \
	-o $(BDIR)/make_two_layer_split_model  $(LDFLAGS)

$(BDIR)/read_itl_sol: read_itl_sol.c   $(ODIR)/calc_two_layer_split.o 
	$(CC) $(CFLAGS) read_itl_sol.c  -o $(BDIR)/read_itl_sol  $(ODIR)/calc_two_layer_split.o \
	$(LDFLAGS)

$(BDIR)/eval_two_layer_split_model: eval_two_layer_split_model.c  $(ODIR)/calc_two_layer_split.o 
	$(CC) $(CFLAGS) eval_two_layer_split_model.c  $(ODIR)/calc_two_layer_split.o \
	-o $(BDIR)/eval_two_layer_split_model  $(LDFLAGS)

$(BDIR)/invert_two_layer_split: invert_two_layer_split.c  $(ODIR)/calc_two_layer_split.o 
	$(CC) $(CFLAGS) invert_two_layer_split.c  $(ODIR)/calc_two_layer_split.o \
	-o $(BDIR)/invert_two_layer_split  $(LDFLAGS)

$(BDIR)/invert_two_layer_split_field: invert_two_layer_split_field.c  $(ODIR)/calc_two_layer_split.o 
	$(CC) $(CFLAGS) invert_two_layer_split_field.c  $(ODIR)/calc_two_layer_split.o \
	-o $(BDIR)/invert_two_layer_split_field  $(LDFLAGS)

$(BDIR)/spherical_spline_invert: $(ODIR)/spherical_spline_invert.o \
	$(ODIR)/spherical_spline.o  spherical_spline.h
	$(CC) $(CFLAGS) $(ODIR)/spherical_spline_invert.o $(ODIR)/spherical_spline.o \
	-o $(BDIR)/spherical_spline_invert $(LDFLAGS) $(MATHLIBS)

$(BDIR)/spherical_spline_evaluate: $(ODIR)/spherical_spline_evaluate.o \
	$(ODIR)/spherical_spline.o  spherical_spline.h
	$(CC) $(CFLAGS) $(ODIR)/spherical_spline_evaluate.o $(ODIR)/spherical_spline.o \
	-o $(BDIR)/spherical_spline_evaluate $(LDFLAGS) $(MATHLIBS)

$(BDIR)/spherical_spline_construct: $(ODIR)/spherical_spline_construct.o \
	$(ODIR)/spherical_spline.o  spherical_spline.h
	$(CC) $(CFLAGS) $(ODIR)/spherical_spline_construct.o $(ODIR)/spherical_spline.o \
	-o $(BDIR)/spherical_spline_construct $(LDFLAGS) $(MATHLIBS)


$(BDIR)/aki2cmt: aki2cmt.f
	$(F77) $(FFLAGS) aki2cmt.f -o $(BDIR)/aki2cmt 

# GMT 6
$(BDIR)/grd2poisson: $(ODIR)/grd2poisson.o
	$(CC) $(CFLAGS) -I$(GMTHOME)/include/ -I$(NETCDFHOME)/include/ \
	-L$(GMTHOME)/lib/ -L$(NETCDFHOME)/lib/ \
	$(ODIR)/grd2poisson.o -o $(BDIR)/grd2poisson \
	-I../fishpak/ -L../fishpak/lib/$(ARCH)/ -lfishpak $(GMT6_LIB) $(LDFLAGS) 

# GMT4
$(BDIR)/grd2poisson_mud: $(ODIR)/grd2poisson_mud.o
	$(CC) $(CFLAGS) -L$(GMT4HOME)/lib/ -L$(NETCDFHOME)/lib/ -L../mudpack/lib/$(ARCH)/ \
	$(ODIR)/grd2poisson_mud.o -o $(BDIR)/grd2poisson_mud \
	-lmudpack $(GMT_LIBS) -lmudpack $(LDFLAGS) 

$(BDIR)/grd2poisson_my: $(ODIR)/grd2poisson_my.o $(ODIR)/nr_util.o $(ODIR)/mglin.o
	$(CC) $(CFLAGS) -I$(GMT4HOME)/include/ -I$(NETCDFHOME)/include/ -L$(GMT4HOME)/lib/ -L$(NETCDFHOME)/lib/ \
	$(ODIR)/grd2poisson_my.o -o $(BDIR)/grd2poisson_my $(ODIR)/nr_util.o $(ODIR)/mglin.o \
	$(GMT_LIBS) $(LDFLAGS) $(LDFLAGS)

$(BDIR)/pca: pca.c $(ODIR)/pca_util.o $(ODIR)/misc.o
	$(CC) $(CFLAGS) pca.c $(ODIR)/pca_util.o $(ODIR)/misc.o -o $(BDIR)/pca \
	$(EISLIB) -lmyeis $(LDFLAGS)

$(BDIR)/median: median.c  $(ODIR)/misc.o
	$(CC) $(CFLAGS) median.c $(ODIR)/misc.o -o $(BDIR)/median $(LDFLAGS)

$(BDIR)/gifmerge: gifmerge.c
	$(CC) $(CFLAGS) gifmerge.c -o $(BDIR)/gifmerge $(LDFLAGS)

$(BDIR)/extract_bin: extract_bin.c $(ODIR)/flip_byte_order.o
	$(CC) $(CFLAGS) extract_bin.c $(ODIR)/flip_byte_order.o -o $(BDIR)/extract_bin $(LDFLAGS)

$(BDIR)/extract_bin_column: extract_bin_column.c $(ODIR)/flip_byte_order.o
	$(CC) $(CFLAGS) extract_bin_column.c $(ODIR)/flip_byte_order.o -o $(BDIR)/extract_bin_column $(LDFLAGS)


$(BDIR)/extract_bin_file: extract_bin_file.c $(ODIR)/flip_byte_order.o
	$(CC) $(CFLAGS) extract_bin_file.c  $(ODIR)/flip_byte_order.o -o $(BDIR)/extract_bin_file $(LDFLAGS)

$(BDIR)/extract_bin_int_file: extract_bin_int_file.c 
	$(CC) $(CFLAGS) extract_bin_int_file.c  \
	-o $(BDIR)/extract_bin_int_file $(LDFLAGS)

$(BDIR)/extract_bin_int: extract_bin_int.c 
	$(CC) $(CFLAGS) extract_bin_int.c   -o $(BDIR)/extract_bin_int $(LDFLAGS)

$(BDIR)/paste_bin: paste_bin.c
	$(CC) $(CFLAGS) paste_bin.c -o $(BDIR)/paste_bin $(LDFLAGS)

$(BDIR)/xyz2lonlatz: xyz2lonlatz.c
	$(CC) $(CFLAGS) xyz2lonlatz.c -o $(BDIR)/xyz2lonlatz $(LDFLAGS)

$(BDIR)/svdfit: $(ODIR)/svdfit.o $(ODIR)/nr_util.o $(ODIR)/svd_util.o
	$(CC) $(CFLAGS) $(ODIR)/svdfit.o $(ODIR)/nr_util.o $(ODIR)/svd_util.o \
	 -o $(BDIR)/svdfit $(LDFLAGS)

$(BDIR)/lonlat2gcpoints: lonlat2gcpoints.c
	$(CC) $(CFLAGS) lonlat2gcpoints.c -o $(BDIR)/lonlat2gcpoints $(LDFLAGS)

$(BDIR)/lonlatzcv2pv: lonlatzcv2pv.c
	$(CC) $(CFLAGS) lonlatzcv2pv.c -o $(BDIR)/lonlatzcv2pv $(LDFLAGS)

$(BDIR)/dist_on_sphere_and_azi: dist_on_sphere_and_azi.c $(ODIR)/stadis.o
	$(CC) $(CFLAGS) dist_on_sphere_and_azi.c $(ODIR)/stadis.o \
	-o $(BDIR)/dist_on_sphere_and_azi $(LDFLAGS)

$(BDIR)/banner: banner.c
	$(CC) $(CFLAGS) banner.c -o $(BDIR)/banner $(LDFLAGS)

$(BDIR)/momt_project: momt_project.c
	$(CC) $(CFLAGS) momt_project.c -o $(BDIR)/momt_project $(LDFLAGS)



$(BDIR)/averageirrdata_spline: averageirrdata.c
	$(CC) $(CFLAGS) -DINTERPOLATION=0 \
		averageirrdata.c \
		-o $(BDIR)/averageirrdata_spline $(LDFLAGS)

$(BDIR)/averageirrdata_linear: averageirrdata.c
	$(CC) $(CFLAGS) -DINTERPOLATION=1 \
		averageirrdata.c \
		-o $(BDIR)/averageirrdata_linear $(LDFLAGS)

$(BDIR)/averageirrdata_quad: averageirrdata.c
	$(CC) $(CFLAGS) -DINTERPOLATION=2 \
		averageirrdata.c \
		-o $(BDIR)/averageirrdata_quad $(LDFLAGS)

$(BDIR)/averageirrdata_cubic: averageirrdata.c
	$(CC) $(CFLAGS) -DINTERPOLATION=3 \
		averageirrdata.c \
		-o $(BDIR)/averageirrdata_cubic $(LDFLAGS)

$(BDIR)/rotatesmat: rotatemat.c $(ODIR)/rotate.o
	$(CC) $(CFLAGS) rotatemat.c  $(DEFS) $(ODIR)/rotate.o -DSPHERICAL -o $(BDIR)/rotatesmat $(LDFLAGS)

$(BDIR)/rotatecmat: rotatemat.c $(ODIR)/rotate.o
	$(CC) $(CFLAGS) rotatemat.c  $(DEFS) $(ODIR)/rotate.o -o $(BDIR)/rotatecmat $(LDFLAGS)

$(BDIR)/rotatesvec: rotatevec.c $(ODIR)/rotate.o
	$(CC) $(CFLAGS) rotatevec.c  $(DEFS) $(ODIR)/rotate.o -DSPHERICAL -o $(BDIR)/rotatesvec $(LDFLAGS)

$(BDIR)/solve_hom3d: $(ODIR)/solve_hom3d.o $(ODIR)/numrec_svd.o
	$(F77) $(LDFLAGS) $(ODIR)/numrec_svd.o \
		$(ODIR)/solve_hom3d.o -o $(BDIR)/solve_hom3d  \
		$(MATHLIBS)  -lm

$(BDIR)/rotatecvec: rotatevec.c  $(ODIR)/rotate.o
	$(CC) $(CFLAGS) rotatevec.c $(DEFS)  $(ODIR)/rotate.o -o $(BDIR)/rotatecvec $(LDFLAGS)

$(BDIR)/unhqx: unhqx.c
	$(CC) $(CFLAGS) unhqx.c -o $(BDIR)/unhqx $(LDFLAGS)

$(BDIR)/correl: $(ODIR)/correl_nr_util.o $(ODIR)/correl_driver.o
	$(CC) $(CFLAGS) $(ODIR)/correl_nr_util.o $(ODIR)/correl_driver.o  \
	-o $(BDIR)/correl $(LDFLAGS)

$(BDIR)/crosscorrelate_delay: $(ODIR)/correl_nr_util.o $(ODIR)/crosscorrelate_delay.o $(ODIR)/correl_util.o 
	$(CC) $(CFLAGS) $(ODIR)/correl_nr_util.o $(ODIR)/crosscorrelate_delay.o $(ODIR)/correl_util.o  \
	-o $(BDIR)/crosscorrelate_delay $(LDFLAGS)

$(BDIR)/crosscorrelate_window: $(ODIR)/crosscorrelate_window.o $(ODIR)/correl_util.o
	$(CC) $(CFLAGS)  $(ODIR)/crosscorrelate_window.o $(ODIR)/correl_util.o  \
	-o $(BDIR)/crosscorrelate_window $(LDFLAGS)

$(BDIR)/gpoints: gpoints.c
	$(CC) $(CFLAGS) gpoints.c -o $(BDIR)/gpoints $(LDFLAGS)

$(BDIR)/invert3x3mat: invert3x3mat.c
	$(CC) $(CFLAGS) invert3x3mat.c -o $(BDIR)/invert3x3mat $(LDFLAGS)



$(BDIR)/extract_dens: extract_dens.c
	$(CC) $(CFLAGS) extract_dens.c -o $(BDIR)/extract_dens $(LDFLAGS)

$(BDIR)/calcms: calcms.c
	$(CC) $(CFLAGS) calcms.c -o $(BDIR)/calcms $(LDFLAGS)

$(BDIR)/dbl2ascii: dbl2ascii.c
	$(CC) $(CFLAGS) dbl2ascii.c -o $(BDIR)/dbl2ascii $(LDFLAGS)

$(BDIR)/resample: resample.c
	$(CC) $(CFLAGS) resample.c -o $(BDIR)/resample $(LDFLAGS)

$(BDIR)/rtp2xyz: rtp2xyz.c
	$(CC) $(CFLAGS) rtp2xyz.c -o $(BDIR)/rtp2xyz $(LDFLAGS)
$(BDIR)/cov_from_events: cov_from_events.c
	$(CC) $(CFLAGS) cov_from_events.c -o $(BDIR)/cov_from_events $(LDFLAGS)

$(BDIR)/stddev: stddev.c
	$(CC) $(CFLAGS) stddev.c -o $(BDIR)/stddev $(LDFLAGS)

$(BDIR)/crosscorrelate: crosscorrelate.c
	$(CC) $(CFLAGS) crosscorrelate.c -o $(BDIR)/crosscorrelate $(LDFLAGS)


$(BDIR)/assign_neighbor_code: assign_neighbor_code.c
	$(CC) $(CFLAGS) assign_neighbor_code.c -o $(BDIR)/assign_neighbor_code $(LDFLAGS)

$(BDIR)/combineoff: combineoff.c
	$(CC) $(CFLAGS) combineoff.c -o $(BDIR)/combineoff $(LDFLAGS)

$(BDIR)/matrix2vector: matrix2vector.c
	$(CC) $(CFLAGS) matrix2vector.c -o $(BDIR)/matrix2vector $(LDFLAGS)

$(BDIR)/simple_spline: simple_spline.c
	$(CC) $(CFLAGS) simple_spline.c -o $(BDIR)/simple_spline $(LDFLAGS)

$(BDIR)/lonlatpv2cv: lonlatpv2cv.c
	$(CC) $(CFLAGS) lonlatpv2cv.c -o $(BDIR)/lonlatpv2cv $(LDFLAGS)

$(BDIR)/multiply_with_transpose: multiply_with_transpose.c
	$(CC) $(CFLAGS) multiply_with_transpose.c -o $(BDIR)/multiply_with_transpose $(LDFLAGS)

$(BDIR)/determine_net_rotation: $(ODIR)/determine_net_rotation.o \
	$(ODIR)/determine_netr_F.o $(ODIR)/read_velocities.o
	$(CC) $(CFLAGS) $(LDFLAGS) $(ODIR)/determine_netr_F.o \
	$(ODIR)/read_velocities.o  $(ODIR)/determine_net_rotation.o  \
	-o $(BDIR)/determine_net_rotation $(FTRN_LIB) $(MATHLIBS) $(LDFLAGS)

$(BDIR)/remove_net_rotation: $(ODIR)/remove_net_rotation.o \
	$(ODIR)/determine_netr_F.o $(ODIR)/read_velocities.o
	$(CC) $(CFLAGS) $(LDFLAGS) $(ODIR)/determine_netr_F.o $(ODIR)/read_velocities.o \
	$(ODIR)/remove_net_rotation.o  \
	-o $(BDIR)/remove_net_rotation $(FTRN_LIB) $(MATHLIBS) $(LDFLAGS)

$(BDIR)/geomean: geomean.c
	$(CC) $(CFLAGS) geomean.c -o $(BDIR)/geomean $(LDFLAGS)

$(BDIR)/misfit: misfit.c
	$(CC) $(CFLAGS) misfit.c -o $(BDIR)/misfit $(LDFLAGS)

$(BDIR)/random: $(ODIR)/random.o $(ODIR)/rand.o $(ODIR)/rotate.o
	$(CC) $(CFLAGS) $(ODIR)/rand.o $(ODIR)/rotate.o $(ODIR)/random.o -o $(BDIR)/random $(LDFLAGS)

$(BDIR)/random_lonlat: $(ODIR)/random_lonlat.o $(ODIR)/rand.o $(ODIR)/rotate.o
	$(CC) $(CFLAGS) $(ODIR)/rand.o $(ODIR)/rotate.o $(ODIR)/random_lonlat.o \
		-o $(BDIR)/random_lonlat $(LDFLAGS)

$(BDIR)/ran_angle: ran_angle.c $(ODIR)/rand.o $(ODIR)/rotate.o
	$(CC) $(CFLAGS) $(ODIR)/rand.o $(ODIR)/rotate.o ran_angle.c -o $(BDIR)/ran_angle $(LDFLAGS)

$(BDIR)/random_gauss: $(ODIR)/random_gauss.o $(ODIR)/rand.o $(ODIR)/rotate.o
	$(CC) $(CFLAGS) $(ODIR)/rand.o $(ODIR)/rotate.o \
		$(ODIR)/random_gauss.o -o $(BDIR)/random_gauss $(LDFLAGS)

$(BDIR)/random_power: $(ODIR)/random_power.o $(ODIR)/rand.o $(ODIR)/rotate.o
	$(CC) $(CFLAGS) $(ODIR)/rand.o $(ODIR)/rotate.o \
		$(ODIR)/random_power.o -o $(BDIR)/random_power $(LDFLAGS) -lm

$(BDIR)/random_tgr: $(ODIR)/random_tgr.o $(ODIR)/rand.o $(ODIR)/rotate.o $(ODIR)/momentmag.o
	$(CC) $(CFLAGS) $(ODIR)/rand.o $(ODIR)/rotate.o  $(ODIR)/momentmag.o \
		$(ODIR)/random_tgr.o -o $(BDIR)/random_tgr $(LDFLAGS) -lm

$(BDIR)/fit_tgr: $(ODIR)/momentmag.o $(ODIR)/fit_tgr.o
	$(CC) $(CFLAGS) $(ODIR)/momentmag.o $(ODIR)/fit_tgr.o -o $(BDIR)/fit_tgr $(LDFLAGS) -lm

$(BDIR)/generate_tgr: $(ODIR)/momentmag.o $(ODIR)/rand.o $(ODIR)/generate_tgr.o
	$(CC) $(CFLAGS) $(ODIR)/momentmag.o $(ODIR)/rand.o \
	$(ODIR)/generate_tgr.o -o $(BDIR)/generate_tgr $(LDFLAGS) -lm

$(BDIR)/random_exp: $(ODIR)/random_exp.o $(ODIR)/rand.o $(ODIR)/rotate.o
	$(CC) $(CFLAGS) $(ODIR)/rand.o $(ODIR)/rotate.o \
		$(ODIR)/random_exp.o -o $(BDIR)/random_exp $(LDFLAGS)

$(BDIR)/random_gvel: $(ODIR)/random_gvel.o $(ODIR)/rand.o $(ODIR)/rotate.o
	$(CC) $(CFLAGS) $(ODIR)/rand.o $(ODIR)/rotate.o \
		$(ODIR)/random_gvel.o -o $(BDIR)/random_gvel $(LDFLAGS)

$(BDIR)/random_tensor: random_tensor.c $(ODIR)/rand.o $(ODIR)/rotate.o $(ODIR)/eigen.o
	$(CC) $(CFLAGS) $(ODIR)/rand.o $(ODIR)/eigen.o $(ODIR)/rotate.o \
	random_tensor.c -o $(BDIR)/random_tensor  $(EISLIB)  $(LDFLAGS) 


$(BDIR)/eigen: $(ODIR)/eigen.main.o $(ODIR)/eigen.o
	$(CC) $(CFLAGS) $(ODIR)/eigen.main.o $(ODIR)/eigen.o \
	-o $(BDIR)/eigen $(EISLIB) $(LDFLAGS) 

$(BDIR)/eigenvalues: $(ODIR)/eigen.ov.o $(ODIR)/eigen.o
	$(CC) $(CFLAGS) $(ODIR)/eigen.ov.o $(ODIR)/eigen.o -o $(BDIR)/eigenvalues \
	$(EISLIB) $(LDFLAGS) 

$(BDIR)/eigen3ds: $(ODIR)/eigen.tds.o $(ODIR)/eigen.o
	$(CC) $(CFLAGS) $(ODIR)/eigen.tds.o  $(ODIR)/eigen.o \
	-o $(BDIR)/eigen3ds $(EISLIB) $(LDFLAGS) 

$(BDIR)/eigenvalues3ds: $(ODIR)/eigen.tds.ov.o $(ODIR)/eigen.o
	$(CC) $(CFLAGS) $(ODIR)/eigen.tds.ov.o $(ODIR)/eigen.o \
	-o $(BDIR)/eigenvalues3ds \
	$(EISLIB) $(LDFLAGS) 

$(BDIR)/determine_coeff: $(ODIR)/chebyshev.o $(ODIR)/numrec_svd.o $(ODIR)/splinesc.o \
		$(ODIR)/splinesf.o $(ODIR)/determine_coeff.o
	$(CC) $(CFLAGS) $(ODIR)/chebyshev.o $(ODIR)/splinesc.o $(ODIR)/splinesf.o \
		$(ODIR)/numrec_svd.o $(ODIR)/determine_coeff.o -o $(BDIR)/determine_coeff \
		$(FTRN_LIB) $(LDFLAGS)  

$(BDIR)/fft: fft.c
	$(CC) $(CFLAGS) fft.c -o $(BDIR)/fft $(LDFLAGS)


$(BDIR)/iaspei91: $(ODIR)/iaspei91.o $(ODIR)/emiask.o
	$(CC) $(CFLAGS)  $(ODIR)/iaspei91.o $(ODIR)/emiask.o \
	-o $(BDIR)/iaspei91 $(LDFLAGS)

$(BDIR)/spat_power_grid: spat_power_grid.c
	$(CC) $(CFLAGS) spat_power_grid.c -o $(BDIR)/spat_power_grid $(LDFLAGS)


$(BDIR)/asciiint2bin: asciiint2bin.c
	$(CC) $(CFLAGS) asciiint2bin.c -o $(BDIR)/asciiint2bin $(LDFLAGS)

$(BDIR)/asciiint2bebin: asciiint2bebin.c $(ODIR)/flip_byte_order.o
	$(CC) $(CFLAGS) asciiint2bebin.c  $(ODIR)/flip_byte_order.o -o $(BDIR)/asciiint2bebin $(LDFLAGS)

$(BDIR)/asciiuint2bin: asciiuint2bin.c
	$(CC) $(CFLAGS) asciiuint2bin.c -o $(BDIR)/asciiuint2bin $(LDFLAGS)

$(BDIR)/asciifloat2bin: asciifloat2bin.c
	$(CC) $(CFLAGS) asciifloat2bin.c -o $(BDIR)/asciifloat2bin $(LDFLAGS)

$(BDIR)/asciifloat2scaled_int: asciifloat2scaled_int.c
	$(CC) $(CFLAGS) asciifloat2scaled_int.c -o $(BDIR)/asciifloat2scaled_int $(LDFLAGS)

$(BDIR)/asciifloat2bebin: asciifloat2bebin.c $(ODIR)/flip_byte_order.o
	$(CC) $(CFLAGS) asciifloat2bebin.c  $(ODIR)/flip_byte_order.o -o $(BDIR)/asciifloat2bebin $(LDFLAGS)

$(BDIR)/comb2col: comb2col.c
	$(CC) $(CFLAGS) -I${HOME}/mylibs/numrec/ comb2col.c -o $(BDIR)/comb2col $(LDFLAGS)

$(BDIR)/fitxy: fitxy.c
	$(CC) $(CFLAGS) -I${HOME}/mylibs/numrec/ fitxy.c -o $(BDIR)/fitxy $(LDFLAGS)

$(BDIR)/fitxyee: fitxyee.c $(ODIR)/fitxyee_util.o 
	$(CC) $(CFLAGS) fitxyee.c $(ODIR)/fitxyee_util.o -o $(BDIR)/fitxyee $(LDFLAGS)

$(BDIR)/p_fisher: p_fisher.c
	$(CC) $(CFLAGS) p_fisher.c -I${HOME}/mylibs/numrec/ -o $(BDIR)/p_fisher $(LDFLAGS)


$(BDIR)/p_student: p_student.c 
	$(CC) $(CFLAGS) p_student.c  -I${HOME}/mylibs/numrec/ -o $(BDIR)/p_student $(LDFLAGS)

$(BDIR)/pearson_corr: pearson_corr.c
	$(CC) $(CFLAGS) pearson_corr.c -I${HOME}/mylibs/numrec/ \
	-o $(BDIR)/pearson_corr $(LDFLAGS)

$(BDIR)/spear: spear.c spear_util.c
	$(CC) $(CFLAGS) spear.c spear_util.c -I${HOME}/mylibs/numrec/ \
	-o $(BDIR)/spear $(LDFLAGS)

$(BDIR)/bootstrap_spear: bootstrap_corr.c spear_util.c
	$(CC) $(CFLAGS)  \
	bootstrap_corr.c spear_util.c -I${HOME}/mylibs/numrec/ \
	-o $(BDIR)/bootstrap_spear $(LDFLAGS)

$(BDIR)/bootstrap_corr: bootstrap_corr.c spear_util.c
	$(CC) $(CFLAGS)  -DUSE_PEARSON \
	bootstrap_corr.c spear_util.c -I${HOME}/mylibs/numrec/ \
	-o $(BDIR)/bootstrap_corr $(LDFLAGS)

$(BDIR)/period: period.c $(MT_LIBS)
	$(CC) $(CFLAGS) period.c -o $(BDIR)/period $(MT_LIBS) $(LDFLAGS)

$(BDIR)/spline_sample: spline_sample.c $(ODIR)/nr_spline.o
	$(CC) $(CFLAGS) spline_sample.c $(ODIR)/nr_spline.o -o $(BDIR)/spline_sample $(LDFLAGS)


$(BDIR)/student: student.c
	$(CC) $(CFLAGS) -I${HOME}/mylibs/numrec/ student.c -o $(BDIR)/student $(LDFLAGS)





# objects from here
$(ODIR)/polysort.o: polysort.c
	$(CC) $(CFLAGS) $(GMT_INC) -c polysort.c -DUNIFORM -o $(ODIR)/polysort.o


$(ODIR)/random.o: my_rand_routines.c
	$(CC) $(CFLAGS) -c my_rand_routines.c -DUNIFORM -o $(ODIR)/random.o

$(ODIR)/random_lonlat.o: my_rand_routines.c
	$(CC) $(CFLAGS) -c my_rand_routines.c -DUNIFORM_LONLAT -o $(ODIR)/random_lonlat.o


$(ODIR)/random_gauss.o: my_rand_routines.c
	$(CC) $(CFLAGS) -c my_rand_routines.c -DGAUSSIAN -o $(ODIR)/random_gauss.o

$(ODIR)/random_power.o: my_rand_routines.c
	$(CC) $(CFLAGS) -c my_rand_routines.c -DPOWERLAW -o $(ODIR)/random_power.o

$(ODIR)/random_tgr.o: my_rand_routines.c
	$(CC) $(CFLAGS) -c my_rand_routines.c -DTGR -o $(ODIR)/random_tgr.o

$(ODIR)/random_exp.o: my_rand_routines.c
	$(CC) $(CFLAGS) -c my_rand_routines.c -DEXPONENTIAL -o $(ODIR)/random_exp.o

$(ODIR)/random_gvel.o: my_rand_routines.c
	$(CC) $(CFLAGS) -c my_rand_routines.c -DGAUSS_VEL -o $(ODIR)/random_gvel.o


$(ODIR)/determine_netr_F.o: determine_netr.F
	$(F77) $(FFLAGS) -c determine_netr.F -o $(ODIR)/determine_netr_F.o

$(ODIR)/determine_net_rotation.o: determine_net_rotation.c
	$(CC) $(CFLAGS) -c determine_net_rotation.c -o $(ODIR)/determine_net_rotation.o

$(ODIR)/remove_net_rotation.o: remove_net_rotation.c
	$(CC) $(CFLAGS) -c remove_net_rotation.c -o $(ODIR)/remove_net_rotation.o

$(ODIR)/chebyshev.o: chebyshev.c
	$(CC) $(CFLAGS) -c chebyshev.c -o $(ODIR)/chebyshev.o

$(ODIR)/splinesc.o: splinesc.c
	$(CC) $(CFLAGS) -c splinesc.c -o $(ODIR)/splinesc.o

$(ODIR)/splinesf.o: splinesf.F
	$(F77) $(FFLAGS) -c splinesf.F -o $(ODIR)/splinesf.o


$(ODIR)/numrec_svd.o: numrec_svd.F
	$(F77) $(FFLAGS) -DUSE_DOUBLE_PRECISION -c numrec_svd.F -o $(ODIR)/numrec_svd.o

$(ODIR)/rotate.o: rotate.c
	$(CC) $(CFLAGS) $(DEFS) -c rotate.c -o $(ODIR)/rotate.o

$(ODIR)/spherical_spline.o: spherical_spline.c spherical_spline.h
	$(CC) $(CFLAGS) $(DEFS) -c spherical_spline.c -o $(ODIR)/spherical_spline.o

$(ODIR)/fit_plane.o: fit_plane.c
	$(CC) $(CFLAGS) -I$(HOME)/progs/src/interact/ $(GMT_INC) \
	-c fit_plane.c -o $(ODIR)/fit_plane.o 

$(ODIR)/svdfit.o: svdfit.c
	$(CC) $(CFLAGS) -c svdfit.c -o $(ODIR)/svdfit.o

$(BDIR)/mpi_example: mpi_example.c
	$(MPICC) mpi_example.c 	-o $(BDIR)/mpi_example

$(ODIR)/determine_coeff.o: determine_coeff.c
	$(CC) $(CFLAGS) -c determine_coeff.c \
	-o $(ODIR)/determine_coeff.o 


$(ODIR)/plume_detect_d.o: plume_detect.c  $(HDR_FLS) 
	$(CC) $(CFLAGS)  -DPD_DOWN -c plume_detect.c $(GMT_INC) $(GGRD_INC) \
	-o $(ODIR)/plume_detect_d.o

$(ODIR)/plume_detect.o: plume_detect.c  $(HDR_FLS) 
	$(CC) $(CFLAGS)   -c plume_detect.c $(GMT_INC) $(GGRD_INC) -o $(ODIR)/plume_detect.o

$(ODIR)/read_velocities.o: read_velocities.c
	$(CC) $(CFLAGS) -c read_velocities.c $(GMT_INC) -o $(ODIR)/read_velocities.o 


$(ODIR)/eigen.main.o: eigen_driver.c
	$(CC) $(CFLAGS) -c eigen_driver.c -o $(ODIR)/eigen.main.o 

$(ODIR)/eigen.ov.o: eigen_driver.c
	$(CC) $(CFLAGS) -c eigen_driver.c -DONLY_VALUES -o $(ODIR)/eigen.ov.o

$(ODIR)/eigen.tds.o: eigen_driver.c
	$(CC) $(CFLAGS) -c eigen_driver.c -DTHREED_SYMMETRIC -o $(ODIR)/eigen.tds.o

$(ODIR)/eigen.tds.ov.o: eigen_driver.c
	$(CC) $(CFLAGS) -c eigen_driver.c -DTHREED_SYMMETRIC \
		-DONLY_VALUES -o $(ODIR)/eigen.tds.ov.o

$(ODIR)/evgtr1.o: evgtr1.f  $(HDR_FLS)
	$(F77) $(FFLAGS) $(F_EXT_SOURCE_FLAG) -c $< -o $(ODIR)/evgtr1.o

$(ODIR)/misc.o: misc.c  $(HDR_FLS)
	$(CC) $(CFLAGS) -c $< -o $(ODIR)/misc.o $(GMT_INC)

$(ODIR)/spherical_spline_invert.o: spherical_spline_invert.c spherical_spline.h
	$(CC) -c $(CFLAGS) spherical_spline_invert.c \
	-o $(ODIR)/spherical_spline_invert.o 

$(ODIR)/spherical_spline_evaluate.o: spherical_spline_evaluate.c spherical_spline.h
	$(CC) -c $(CFLAGS) spherical_spline_evaluate.c \
	-o $(ODIR)/spherical_spline_evaluate.o 

$(ODIR)/spherical_spline_construct.o: spherical_spline_construct.c spherical_spline.h
	$(CC) -c $(CFLAGS) spherical_spline_construct.c \
	-o $(ODIR)/spherical_spline_construct.o 


$(BDIR)/fit_exp_scan: fit_exp_scan.c 
	$(CC) $(CFLAGS) fit_exp_scan.c $(LMFIT_INC) -o $(BDIR)/fit_exp_scan \
		$(LDFLAGS) $(LMFIT_LIBS)  -lm

$(BDIR)/prune_gps_data: prune_gps_data.c 
	$(CC) $(CFLAGS) prune_gps_data.c  -o $(BDIR)/prune_gps_data  $(LDFLAGS) 



$(ODIR)/%.o: %.F   $(HDR_FLS)
	$(F77) $(FFLAGS) -c $< -o $(ODIR)/$*.o


$(ODIR)/%.o: %.f90  $(HDR_FLS)
	$(F90) $(FFLAGS) -c $< -o $(ODIR)/$*.o


$(ODIR)/%.o: %.f  $(HDR_FLS)
	$(F77) $(FFLAGS) -c $< -o $(ODIR)/$*.o

$(ODIR)/%.o: %.c  $(HDR_FLS)
	$(CC) $(CFLAGS) -c $< -o $(ODIR)/$*.o

