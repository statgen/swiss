FROM alpine:3.8

# Install required packages for swiss to install. Many of swiss' dependencies
# require compiling C/C++ code.
RUN apk update && \
  apk add python2 python2-dev py2-pip py2-virtualenv git gcc g++ make \
  lzo lzo-dev zlib zlib-dev bzip2 bzip2-dev xz xz-dev curl curl-dev bash

# Grab the latest htslib for tabix.
RUN cd /root && \
  mkdir -p install/htslib && \
  cd install/htslib && \
  curl -OJL https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2 && \
  tar xf htslib-1.9.tar.bz2 && \
  cd htslib-1.9 && \
  make install

# Grab PLINK.
RUN cd /root && \
  mkdir -p install/plink && \
  cd install/plink && \
  curl -OJ https://www.cog-genomics.org/static/bin/plink180807/plink_linux_x86_64.zip && \
  unzip plink_linux_x86_64.zip && \
  cp plink /usr/local/bin

# Upgrade pip to latest.
RUN pip install --upgrade pip

# Install swiss.
RUN pip install git+https://github.com/statgen/swiss.git@v1.0.0

# Create a swiss group and user to execute as.
RUN addgroup swiss && adduser -G swiss -s /bin/bash -D swiss
WORKDIR /home/swiss
USER swiss

# Download necessary swiss data (LD, GWAS catalogs, etc.)
#RUN swiss --download-data

# Run the swiss executable.
ENTRYPOINT ["bash"]
