#!/bin/sh

RESULTDIR=result/

if [ ! -d ${RESULTDIR} ];
then
    mkdir ${RESULTDIR}
fi

P=${PBS_NP}

. ../params.sh

TIMEFILE=${RESULTDIR}/${P}

mpirun ${MPIRUN_PARAMS} ./mpi_master_worker 1 0 10 100000000 1000 2> ${TIMEFILE}

process_time_file ${TIMEFILE}
