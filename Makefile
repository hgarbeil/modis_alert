CC = g++ -std=c++11 -g

hdflibs = -L/hbeta/harold/lhome/external/hdf-4.2.13/lib -lm -lmfhdf -ldf -ljpeg -lz -lsz
hdfinc = -I/hbeta/harold/lhome/external/hdf-4.2.13/include 

OBJS = main.o  modis_hdf.o  surftemp.o modis_process.o

modis_process : ${OBJS}
	${CC} ${OBJS} ${hdflibs} -o modis_process

.cpp.o :
	${CC} ${hdfinc} -c $*.cpp
