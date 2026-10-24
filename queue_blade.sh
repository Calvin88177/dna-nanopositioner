#!/bin/bash
# three-point bending of the 128-bp blade: +70 pN then -70 pN at mid-span
# supports at base-pair planes 34 and 93 (span 59 bp = 20.1 nm), load on planes 63 and 64
set -e
OX=/home/claude/nanopositioner/oxDNA/build/bin/oxDNA
SIM=/home/claude/nanopositioner/sim
cd "$SIM/blade128"
for run in "70 Pp70" "-70 Pm70"; do
    set -- $run
    python3 "$SIM/setup_bend.py" blade128.json 128 34 93 63 64 "$1" "$2" start.dat 120000 500 > "setup_$2.log"
    "$OX" "$2_input" > "$2.out" 2>&1
    echo "$2 done $(date +%H:%M:%S)" >> prod_status.txt
done
