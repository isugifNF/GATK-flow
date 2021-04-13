FROM broadinstitute/gatk:4.1.3.0

RUN apt-get update
RUN apt-get install -y build-essential wget curl git autoconf automake
RUN apt-get install -y gcc g++ bison make
RUN apt-get install -y bwa datamash vcftools

RUN cd /opt; git clone https://github.com/lh3/bioawk.git; cd bioawk; make; mv bioawk maketab /usr/bin/
RUN cd /opt; curl -L https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.0pre2/bwa-mem2-2.0pre2_x64-linux.tar.bz2 | tar jxf - 

ENV PATH=$PATH:/opt/bwa-mem2-2.0pre2_x64-linux

CMD ["ls"]

LABEL author="Jennifer Chang"
LABEL last-update="2021-04-12"
