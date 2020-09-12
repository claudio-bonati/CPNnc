bin_PROGRAMS = cpn_nc

conf = ./include/conf.h \
       ./lib/conf_def.c \
       ./lib/conf_meas.c \
       ./lib/conf_upd.c 

dSFMT = ./include/dSFMT-common.h \
        ./include/dSFMT.h \
        ./include/dSFMT-params11213.h \
        ./include/dSFMT-params1279.h \
        ./include/dSFMT-params132049.h \
        ./include/dSFMT-params19937.h \
        ./include/dSFMT-params216091.h \
        ./include/dSFMT-params2203.h \
        ./include/dSFMT-params4253.h \
        ./include/dSFMT-params44497.h \
        ./include/dSFMT-params521.h \
        ./include/dSFMT-params86243.h \
        ./include/dSFMT-params.h \
        ./lib/dSFMT.c

endian = ./include/endianness.h  \
         ./lib/endianness.c 

flavourmatrix = ./include/flavour_matrix.h \
                ./lib/flavour_matrix.c

geometry = ./include/geometry.h \
           ./lib/geometry.c

gparam = ./include/gparam.h \
         ./lib/gparam.c

random = ./include/myrngparam.h \
         ./include/random.h \
         ./lib/random.c 

vec = ./include/vec.h \
      ./lib/vec.c 

cpn_nc_SOURCES=./src/cpn_nc.c \
               $(conf) \
               $(dSFMT) \
               $(endian) \
               $(flavourmatrix) \
               $(geometry) \
               $(gparam) \
               ./include/macro.h \
               $(random) \
               $(vec) 

AM_CFLAGS = -std=c99 -O2 -Wall -Wextra -Werror -pedantic -Wconversion

CLEANFILES = dati.dat log.dat dati.log conf.dat* template_input.in 
