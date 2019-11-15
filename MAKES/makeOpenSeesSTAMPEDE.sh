#!/bin/sh
module load petsc
make wipe
make PROGRAMMING_MODE=SEQUENTIAL -j 2
make wipe
make PROGRAMMING_MODE=PARALLEL -j 2
make wipe
make PROGRAMMING_MODE=PARALLEL_INTERPRETERS -j 2
chmod 'a+rx' $HOME/bin/OpenSees
chmod 'a+rx' $HOME/bin/OpenSeesSP
chmod 'a+rx' $HOME/bin/OpenSeesMP
