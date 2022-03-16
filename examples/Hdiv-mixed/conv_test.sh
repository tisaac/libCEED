#!/bin/bash

# Copyright (c) 2017, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory. LLNL-CODE-734707.
# All Rights reserved. See files LICENSE and NOTICE for details.
#
# This file is part of CEED, a collection of benchmarks, miniapps, software
# libraries and APIs for efficient high-order finite element and spectral
# element discretizations for exascale applications. For more information and
# source code availability see http://github.com/ceed.
#
# The CEED research is supported by the Exascale Computing Project 17-SC-20-SC,
# a collaborative effort of two U.S. Department of Energy organizations (Office
# of Science and the National Nuclear Security Administration) responsible for
# the planning and preparation of a capable exascale ecosystem, including
# software, applications, hardware, advanced system engineering and early
# testbed platforms, in support of the nation's exascale computing imperative.

declare -A run_flags
    run_flags[pc_type]=svd
    run_flags[problem]=darcy3d
    run_flags[dm_plex_dim]=3
    run_flags[dm_plex_box_faces]=2,2,2

declare -A test_flags
    test_flags[res_start]=2
    test_flags[res_stride]=1
    test_flags[res_end]=6

file_name=conv_test_result.csv

echo "run,mesh_res,error_u,error_p" > $file_name

i=0

for ((res=${test_flags[res_start]}; res<=${test_flags[res_end]}; res+=${test_flags[res_stride]})); do
    run_flags[dm_plex_box_faces]=$res,$res,$res
    args=''
    for arg in "${!run_flags[@]}"; do
        if ! [[ -z ${run_flags[$arg]} ]]; then
            args="$args -$arg ${run_flags[$arg]}"
        fi
    done
    ./main $args | grep "L2 Error of u and p" | awk -v i="$i" -v res="$res" '{ printf "%d,%d,%.5f,%.5f\n", i, res, $8, $9}' >> $file_name
    i=$((i+1))
done

