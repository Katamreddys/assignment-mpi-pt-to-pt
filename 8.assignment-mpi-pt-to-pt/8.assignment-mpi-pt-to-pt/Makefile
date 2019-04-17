
assignment-mpi-pt-to-pt.pdf: assignment-mpi-pt-to-pt.tex
	pdflatex assignment-mpi-pt-to-pt.tex

assignment-mpi-pt-to-pt.tgz: approx.cpp params.sh assignment-mpi-pt-to-pt.pdf \
                             pingpong/mpi_ping_pong.cpp pingpong/Makefile pingpong/run.sh \
                             heat/libgen.a heat/Makefile heat/mpi_heat.cpp heat/plot.sh heat/run-strong.sh heat/run-weak.sh heat/table-weak.sh heat/table-strong.sh heat/test.sh \
                             master_worker/cases.txt master_worker/libfunctions.a master_worker/Makefile master_worker/mpi_master_worker.cpp master_worker/plot.sh master_worker/queue.sh master_worker/run_numint.sh master_worker/test.sh
	tar zcvf assignment-mpi-pt-to-pt.tgz \
                             approx.cpp params.sh assignment-mpi-pt-to-pt.pdf \
                             pingpong/mpi_ping_pong.cpp pingpong/Makefile pingpong/run.sh \
                             heat/libgen.a heat/Makefile heat/mpi_heat.cpp heat/plot.sh heat/run-strong.sh heat/run-weak.sh heat/table-weak.sh heat/table-strong.sh heat/test.sh \
                             master_worker/cases.txt master_worker/libfunctions.a master_worker/Makefile master_worker/mpi_master_worker.cpp master_worker/plot.sh master_worker/queue.sh master_worker/run_numint.sh master_worker/test.sh
