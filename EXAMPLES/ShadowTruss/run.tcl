#!/home/fmk/bin/OpenSees

set port 2222;
set machine 127.0.0.1;  # local machine

exec ./serverTruss $port &
exec ./example1 $port $machine &


