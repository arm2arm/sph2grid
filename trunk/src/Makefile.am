check_PROGRAMS =  sph2grid 

bin_PROGRAMS = sph2grid 
sph2grid_SOURCES = main.cpp sphvolume.cpp sphvolume.h \
		    data_readers.cpp data_readers.h\
		    utils.h utils.cpp\
                    Render.cpp Render.h
# what flags you want to pass to the C compiler & linker

AM_CPPFLAGS = -I/usr/X11R6/include -I/home/arm2arm/Projects/CImg
              
AM_LDFLAGS = -lpthread -lX11 -lXext -lXrandr



