#! /usr/bin/env bash

bsub -J shell -q normal -n 72 -R span[ptile=72] -Is /bin/bash
