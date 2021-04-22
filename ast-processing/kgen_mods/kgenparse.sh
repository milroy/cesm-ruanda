#!/bin/bash 
# 
# kgen_run.sh 

CAM5_HOME=/glade/u/home/milroy/cesm/cesm1_3_beta11/models/atm/cam/
KGEN=/glade/u/home/milroy/git/KGen/bin/kgen
CASE_DIR=/glade/u/home/milroy/cesm/cesm1_3_beta11/scripts/ufect.kgenparse.000/
OUTPUT_DIR=/glade/scratch/milroy/ufect.kgenparse.000/

while read i
do
    python ${KGEN} --outdir ${OUTPUT_DIR} \
     --cmd-build "cd ${CASE_DIR}; ./ufect.kgenparse.000.clean_build; ./ufect.kgenparse.000.build" \
     --cmd-run "cd ${CASE_DIR}; ./ufect.kgenparse.000.submit" \
     $i 2>/dev/null
done < mods.ini
# /glade/u/home/milroy/cesm/cesm1_3_beta11/models/atm/cam/src/physics/cam/micro_mg_cam.F90:

