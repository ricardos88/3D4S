CODENAME	= main
LOCDIR		= /home/ricardo/Dropbox/Codigo
NP		= 1
SIM_NUM         = 16
EXTA		= .pbs
EXTB		= .sh

# -ksp_type gmres -pc_type jacobi
# -ksp_type preonly -pc_type ilu
# -ksp_type gmres -pc_type ilu
# -ksp_rtol 1.e-14 -snes_rtol 1.e-14
# -ts_adapt_type none -snes_max_it 200
# -ts_monitor -snes_monitor -ksp_monitor
# -ts_ssp_type rks2

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

all: delete create compile run

all1: delete create compile run1

allbarkla: delete create compile barkla

delete: 
	rm -r results/${SIM_NUM};

create: 
	mkdir results/${SIM_NUM};

compile: ${CODENAME}.o  chkopts
	-${CLINKER} -o ${CODENAME} ${CODENAME}.o ${PETSC_KSP_LIB}
	${RM} ${CODENAME}.o

run:
	-@${MPIEXEC} -n ${NP} ./${CODENAME} ${SIM_NUM} -ts_adapt_type none -log_view > ${CODENAME}_${SIM_NUM}.tmp 2>&1 &

run1:
	./${CODENAME} ${SIM_NUM} -ts_adapt_type none -log_view > ${CODENAME}_${SIM_NUM}.tmp 2>&1 &

barkla:
	sbatch job${NP}${EXTB} ${SIM_NUM};

sbarkla:
	sbatch sjob${NP}${EXTB} ${SIM_NUM};

carcher: ${CODENAME}.o  chkopts
	-${CLINKER} -o ${CODENAME}${NP} ${CODENAME}.o  ${PETSC_KSP_LIB}
	${RM} ${CODENAME}.o

archer:
	qsub ${CODENAME}${NP}${EXTA}

include ${PETSC_DIR}/lib/petsc/conf/test
