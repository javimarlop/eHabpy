# sh conf_grass7eHabplus.sh
# make distclean
# make -j2
# sudo make install


CFLAGS="-Wall" ./configure \
 --prefix=/home/javier/hierba706 \
 --enable-largefile \
 --with-sqlite \
 --with-mysql=yes --with-mysql-includes="/usr/include/mysql" \
 --with-freetype=yes \
 --with-freetype-includes=/usr/include/freetype2 \
 --with-python= \
 --with-cxx \
 --with-tcltk-includes=/usr/include/tcl8.5 \
 --enable-64bit \
 --with-gdal=/usr/bin/gdal-config \
 --with-proj --with-proj-share=/usr/share/proj/ \
 --with-geos=/usr/bin/geos-config \
 --with-nls \
 --with-opengl-includes=/usr/include/GL/ \
 --with-x \
 --with-fftw \
 --with-wx=/usr/bin/wx-config \
 --with-wxwidgets=/usr/bin/wx-config \
 --with-cairo \
 --enable-shared \
 --with-pthread \
 --with-lapack \
 --with-blas
