FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    gfortran python3 make binutils \
 && rm -rf /var/lib/apt/lists/*

WORKDIR /WAVE
COPY . /WAVE

RUN cd WAVE && python3 python/make_wave.py

CMD ["bash"]
