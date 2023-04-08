FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive

ARG star_version=2.7.10a_alpha_220818
ARG samtools_version=1.15.1
ARG bbmap_version=38.97
ARG rsem_version=1.3.3
ARG subread_version=2.0.2
ARG salmon_version=1.10.0

#Install OS packages
RUN apt-get update && apt-get -y --no-install-recommends -qq install \
    wget gcc build-essential software-properties-common libz-dev \
    git libncurses5-dev libbz2-dev liblzma-dev default-jre bsdmainutils

#Install STAR
RUN wget --no-check-certificate https://github.com/alexdobin/STAR/archive/${star_version}.tar.gz && \
    tar -xzf ${star_version}.tar.gz -C /opt && \
    cd /opt/STAR-${star_version}/source && \
    make STAR CXXFLAGS_SIMD="-msse4.2" && \
    cd / && rm ${star_version}.tar.gz 

#Install seqtk
RUN git clone https://github.com/lh3/seqtk.git && \
    mv seqtk /opt && \
    cd /opt/seqtk && \
    make

#Install samtools
RUN wget https://github.com/samtools/samtools/releases/download/${samtools_version}/samtools-${samtools_version}.tar.bz2 && \
    tar -xvf samtools-${samtools_version}.tar.bz2 -C /opt && \
    cd /opt/samtools-${samtools_version} && \
    ./configure && \
    make && \
    make install && \
    cd / && rm samtools-${samtools_version}.tar.bz2  

#Install BBMap
RUN wget https://sourceforge.net/projects/bbmap/files/BBMap_${bbmap_version}.tar.gz && \
    tar -xzf BBMap_${bbmap_version}.tar.gz -C /opt && \
    cd /opt/bbmap && \
    ./stats.sh in=resources/phix174_ill.ref.fa.gz && \
    cd / && rm BBMap_${bbmap_version}.tar.gz

#Install RSEM
RUN wget https://github.com/deweylab/RSEM/archive/refs/tags/v${rsem_version}.tar.gz && \
    tar -xzf v${rsem_version}.tar.gz -C /opt && \
    cd /opt/RSEM-${rsem_version} && \
    make && \
    cd / && rm v${rsem_version}.tar.gz

#Install Subread (for featureCounts)
RUN wget https://github.com/ShiLab-Bioinformatics/subread/releases/download/${subread_version}/subread-${subread_version}-Linux-x86_64.tar.gz && \ 
    tar -xvf subread-${subread_version}-Linux-x86_64.tar.gz -C /opt && \
    cd / && rm subread-${subread_version}-Linux-x86_64.tar.gz

#Install Salmon
RUN wget https://github.com/COMBINE-lab/salmon/releases/download/v${salmon_version}/salmon-${salmon_version}_linux_x86_64.tar.gz && \
    tar -xvf salmon-${salmon_version}_linux_x86_64.tar.gz -C /opt && \
    mv /opt/salmon-latest_linux_x86_64 /opt/salmon-${salmon_version}_linux_x86_64 && \ 
    cd / && rm salmon-${salmon_version}_linux_x86_64.tar.gz

ENV PATH="${PATH}:/opt/STAR-${star_version}/source:/opt/seqtk:/opt/bbmap:/opt/RSEM-${rsem_version}:/opt/salmon-${salmon_version}_linux_x86_64/bin:/opt/subread-${subread_version}-Linux-x86_64/bin"     

#Saving Software Versions to a file
RUN echo "STAR version: ${star_version}" >> versions.txt && \
    echo "Samtools version: ${samtools_version}" >> versions.txt && \
    echo "BBMap version: ${bbmap_version}" >> versions.txt && \
    echo "RSEM version: ${rsem_version}" >> versions.txt && \
    echo "Salmon version: ${salmon_version}" >> versions.txt && \
    echo "Subread version: ${subread_version}" >> versions.txt && \
    seqtk_version=`strings $(which seqtk) | grep 'Version:' | cut -f 2 -d " "` && \
    echo "seqtk version: ${seqtk_version}" >> versions.txt 
