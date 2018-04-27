
SHELL = /bin/sh
DEBUG = false
#Adapt to host system
ifeq (,$(findstring Windows,$(OS)))
  RM = rm -f
else
  RM = del
endif



all: doxygen

DOXYGEN=doxygen

doxygen:
	$(DOXYGEN) doxy_conf
	cd latex && make
	mv ./latex/refman.pdf ./latex/Reference_manual.pdf
	cd latex && make clean

remove:
	make clean
	$(RM) -r html latex

clean:
	$(RM) *.bak


.PHONY: remove clean all 
