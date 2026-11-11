FROM ubuntu:22.04
RUN apt-get update \
 && apt-get install -y --no-install-recommends ca-certificates wget \
 && rm -rf /var/lib/apt/lists/*
RUN wget -q -O /tmp/subread.tar.gz "https://sourceforge.net/projects/subread/files/subread-1.4.6-p3/subread-1.4.6-p3-Linux-x86_64.tar.gz/download" \
 && tar -xzf /tmp/subread.tar.gz -C /opt \
 && rm /tmp/subread.tar.gz
ENV PATH=/opt/subread-1.4.6-p3-Linux-x86_64/bin:$PATH
