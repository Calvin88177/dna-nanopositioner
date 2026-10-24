#!/bin/bash
# two coaxial duplexes pressed end to end: a force series at 25 C
OX=/home/claude/nanopositioner/oxDNA/build/bin/oxDNA
SIM=/home/claude/nanopositioner/sim
cd "$SIM/duplex"
for F in 70 140 35 0 210; do
    python3 "$SIM/setup_duplex.py" inputs "F$F" dup.json.oxdna "$F" 2000000 2000 25C 3 > "setup_F$F.log"
    "$OX" "F${F}_input" > "F$F.out" 2>&1
    echo "F$F done $(date +%H:%M:%S)" >> status.txt
done
