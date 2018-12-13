FROM ubuntu:18.04

RUN apt-get update 
RUN apt-get -y install sudo g++ gcc gfortran make cmake git vim \
	python3 python3-pip python3-dev \
	openssh-server
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y tcl8.6 tcl8.6-dev
RUN apt-get clean

RUN useradd -m docker && echo "docker:docker" | chpasswd && adduser docker sudo

RUN sed -i 's/PermitRootLogin prohibit-password/PermitRootLogin yes/' /etc/ssh/sshd_config
RUN sed 's@session\s*required\s*pam_loginuid.so@session optional pam_loginuid.so@g' -i /etc/pam.d/sshd
RUN mkdir /var/run/sshd

RUN cd /home/docker && git clone https://github.com/frqc/OpenSees.git
RUN mkdir /home/docker/bin /home/docker/lib
RUN cd /home/docker/OpenSees && cp MAKES/Makefile.def.DOCKER Makefile.def && make
RUN cd /home/docker/OpenSees/SRC/interpreter/ && make 
RUN cp /home/docker/OpenSees/SRC/interpreter/opensees.so /usr/lib/python3/dist-packages


ENV NOTVISIBLE "in users profile"
RUN echo "export VISIBLE=now" >> /etc/profile

RUN chsh -s /bin/bash

EXPOSE 22
CMD ["/usr/sbin/sshd", "-D"]

# USER docker
# CMD /bin/bash

