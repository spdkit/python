#! /usr/bin/env bash
#BSUB -J scheduler
#BSUB -q normal
#BSUB -n 72
#BSUB -o out.%J.txt
#BSUB -e error.%J.txt
#BSUB -R span[ptile=72]

# step0: load env vars
source ~/apps/env.rc
apps load python

# step 1: start remote execution services
(
# start a scheduler for incoming computation requests
gosh-remote -v bootstrap as-scheduler &
sleep 1
# also start a worker for computation in the same node
gosh-remote -v bootstrap as-worker
) 2>&1 | tee gosh-remote.log &

# step 2: run job
sleep 1
./scripts/run-jobs.py

# step 3: when job done, kill background services
sleep 1
pkill gosh-remote
