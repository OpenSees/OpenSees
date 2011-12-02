include ../../../Makefile.def


OBJS       = Subdomain.o SubdomainNodIter.o ShadowSubdomain.o ActorSubdomain.o

# ShadowSubdomain.o ShadowSubdomainActor.o ActorSubdomain.o

# Compilation control

all:	$(OBJS) 

actor:
	$(LINKER) $(LINKFLAGS) ShadowSubdomainActor.o ActorSubdomain.o \
	$(OO_LIBRARY) $(LINKLIBS)  -o ShadowSubdomainActor 

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o ShadowSubdomainActor

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
