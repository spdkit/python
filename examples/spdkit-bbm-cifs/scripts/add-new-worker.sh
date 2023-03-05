#! /usr/bin/env bash

bsub <<EOF
#BSUB -J worker
#BSUB -q normal
#BSUB -n 72
#BSUB -o out.%J.txt
#BSUB -e error.%J.txt
#BSUB -R span[ptile=72]

# load env vars for apps
source ~/apps/env.rc

# start a worker for computation
gosh-remote -v bootstrap as-worker
EOF
