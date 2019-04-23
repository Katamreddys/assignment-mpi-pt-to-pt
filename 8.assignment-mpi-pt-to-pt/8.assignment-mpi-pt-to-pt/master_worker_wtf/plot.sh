#!/bin/sh


RESULTDIR=result/
PLOTDIR=plot/
#h=`hostname`


if [ ! -d ${RESULTDIR} ];
then
    mkdir ${RESULTDIR}
fi 

if [ ! -d ${PLOTDIR} ];
then
    mkdir ${PLOTDIR}
fi 



#preparing files for plotting

#for P in 1 2 4 8 16 32 36; 
#for P in 2 4 8 16 32 36; #There is no P=1 on master worker
for P in 2 4 8 16 24 32; #There is no P=1 on master worker
do
    echo ${P} \
	 $(cat ${RESULTDIR}/1) \
	 $(cat ${RESULTDIR}/${P})
done > ${RESULTDIR}/speedup_master_worker



# create plots

gnuplot <<EOF
set terminal pdf
set output '${PLOTDIR}/master_worker_speedup.pdf'
set style data linespoints
set key top left;
set xlabel 'P';
set ylabel 'speedup';
set xrange [0:40];
set yrange [0:40];
set ytics 2;
set xtics 2;

plot '${RESULTDIR}/speedup_master_worker' u 1:(\$2/\$3) t 'speedup' lc 1 lw 3;

EOF
