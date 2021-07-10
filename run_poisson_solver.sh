#!/bin/sh

export FF_HOME="./fastflow"

SOLVER=jacobi_seq # jacobi_seq, gauss_seq, jacobi_par, gauss_par
make "$SOLVER"

HEIGHT=1000
WIDTH=1000
EPS=1e-12
ITERS=1000

# required when running par methods, otherwise keep defaults sequential solvers are unaffected 
NUM_WORKERS=2
SCHED=0 # static scheduling, that is, ~(#iteration_space/num_workers) (see ff docs)
BLOCK_HEIGHT=1
BLOCK_WIDTH=1

echo "\n\n"
echo "* Running solver ${SOLVER} *\n"
echo "grid size = ${HEIGHT} x ${WIDTH}"
echo "iters = ${ITERS} / no. workers = ${NUM_WORKERS} / scheduling = ${SCHED} / block_height = ${BLOCK_HEIGHT} / block_width = ${BLOCK_WIDTH}"

./$SOLVER \
	$HEIGHT $WIDTH $EPS $NUM_WORKERS $SCHED $BLOCK_HEIGHT $BLOCK_WIDTH $ITERS

echo "\n"
