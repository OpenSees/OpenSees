include ../../Makefile.def

HOST            = Host
NODE		= Node

PROGRAM         = $(Host) $(Node)

HOST_OBJS       = Host.o 
NODE_OBJS       = Node.o 

# Compilation control

PROGRAM = /tmp/go

OBJS = 

all:        
	@$(CD) $(FE)/actor/channel; $(MAKE);
	@$(CD) $(FE)/actor/machineBroker; $(MAKE);
	@$(CD) $(FE)/actor/message; $(MAKE);
	@$(CD) $(FE)/actor/address; $(MAKE);
	@$(CD) $(FE)/actor/actor; $(MAKE);
	@$(CD) $(FE)/actor/shadow; $(MAKE);	
	@$(CD) $(FE)/actor/objectBroker; $(MAKE);

$(HOST):     $(HOST_OBJS) 
	     @echo "Linking $(PROGRAM)..."
	     $(LINKER) $(LINKFLAGS) $(HOST_OBJS) $(FE_LIBRARY) \
	     $(MACHINE_LINKLIBS) $(BLAS_LIBS) \
	     $(MACHINE_NUMERICAL_LIBS) $(MACHINE_SPECIFIC_LIBS) \
	     -o $(HOST)

$(NODE):     $(NODE_OBJS) 
	     @echo "Linking $(PROGRAM)..."
	     $(LINKER) $(LINKFLAGS) $(NODE_OBJS) $(FE_LIBRARY) \
	     $(MACHINE_LINKLIBS) $(BLAS_LIBS) \
	     $(MACHINE_NUMERICAL_LIBS) $(MACHINE_SPECIFIC_LIBS) \
	     -o $(NODE)

tidy:
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean:  tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean
	@$(RM) $(RMFLAGS) fake $(HOST) $(NODE)

wipe: spotless
	@$(CD) $(FE)/actor/channel; $(MAKE) wipe;
	@$(CD) $(FE)/actor/objectBroker; $(MAKE) wipe;
	@$(CD) $(FE)/actor/machineBroker; $(MAKE) wipe;
	@$(CD) $(FE)/actor/message; $(MAKE) wipe;
	@$(CD) $(FE)/actor/address; $(MAKE) wipe;
	@$(CD) $(FE)/actor/shadow; $(MAKE) wipe;
	@$(CD) $(FE)/actor/actor; $(MAKE) wipe;

# DO NOT DELETE THIS LINE -- make depend depends on it.



