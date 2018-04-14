

include ../../../GDALmake.opt

OBJ	=  ogrflatgeobufdataset.o ogrflatgeobuflayer.o packedrtree.o

CPPFLAGS	:=	-I. -I.. -I../..  $(CPPFLAGS)

default:	$(O_OBJ:.o=.$(OBJ_EXT))

clean:
	rm -f *.o $(O_OBJ)

$(O_OBJ):	ogr_flatgeobuf.h

