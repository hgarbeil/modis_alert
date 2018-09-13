CC = g++ -std=c++1y -g

hdflibs = -L/hbeta/harold/lhome/external/hdf-4.2.13/lib -lm -lmfhdf -ldf -ljpeg -lz -lsz -lboost_system -lboost_filesystem
hdfinc = -I/hbeta/harold/lhome/external/hdf-4.2.13/include 

OBJS = main.o  modis_hdf.o  surftemp.o modis_process.o
OBJS_alert = main_alert.o  modis_hdf.o  surftemp.o modis_process.o sinuProjection.o

modis_alert : ${OBJS_alert}
	${CC} ${OBJS_alert} ${hdflibs} -o modis_alert

modis_process : ${OBJS}
	${CC} ${OBJS_alert} ${hdflibs} -o modis_process

.cpp.o :
	${CC} ${hdfinc} -c $*.cpp
