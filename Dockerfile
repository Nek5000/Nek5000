FROM debian:jessie

# Setup build environment
RUN apt-get clean && apt-get update && apt-get install -y --fix-missing \
  bzip2 \
  gfortran \
  git \
  liblapack-dev \
  libmpich2-dev \
  make \
  mpich2 \
  wget \
&& rm -rf /var/lib/apt/lists/*

ENV HOME /home/nek
ENV PYTHONUNBUFFERED 1
ENV VERBOSE_TESTS "True"

# More stuff for testing
WORKDIR /home/nek
RUN wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh -O miniconda.sh
RUN bash miniconda.sh -b -p $HOME/miniconda
ENV PATH $HOME/miniconda/bin:$PATH
RUN hash -r
RUN conda config --set always_yes yes --set changeps1 no
RUN conda update -q conda
RUN conda info -a
RUN conda install pytest
#RUN pip install unittest

# The tests and tools
RUN git clone https://github.com/RonRahaman/NekTests-deprecated.git -b unittests_plus_examples --recursive nek_tests

# Pull in the sources
WORKDIR /home/nek/nek5000
ADD . . 
WORKDIR /home/nek/nek_tests

# Set env 
ENV SOURCE_ROOT /home/nek/nek5000
ENV TOOLS_ROOT /home/nek/nek5000/tools
ENV TOOLS_BIN /home/nek/nek5000/tools/bin
ENV SCRIPTS_ROOT /home/nek/nek5000/bin
ENV EXAMPLES_ROOT /home/nek/nek_tests/examples
ENV LOG_ROOT /home/nek/nek_tests/examples/TestLogs
ENV F77 mpif77
ENV CC mpicc
ENV IFMPI true
ENV PYTHONPATH "/home/nek/nek_tests/:${PYTHONPATH}"

# Run some tests
ENTRYPOINT ["py.test", "-s", "-v"] 
