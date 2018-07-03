FROM reskyner/ccp4
RUN apt-get update
RUN apt-get install -y r-base ncbi-blast+
RUN apt-get install -y tcsh
RUN wget https://pdb-redo.eu/software/pdb-redo.tar.bz2 
RUN mv pdb-redo.tar.bz2 /usr/local/pdb-redo.tar.bz2
RUN bunzip2 /usr/local/pdb-redo.tar.bz2
RUN tar -zvf  /usr/local/pdb-redo.tar
ADD ccp4.setup-csh /usr/local/pdb-redo/ccp4.setup-csh
ADD pdb_redo.csh /usr/local/pdb-redo/pdb_redo.csh
