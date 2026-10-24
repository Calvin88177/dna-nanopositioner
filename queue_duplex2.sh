#!/bin/bash
# second force series for the two-duplex contact: narrows the threshold, repeats 70 pN, and adds a cold run
OX=/home/claude/nanopositioner/oxDNA/build/bin/oxDNA
SIM=/home/claude/nanopositioner/sim
cd "$SIM/duplex"
while ! grep -q "F210 done" status.txt 2>/dev/null; do sleep 20; done
for run in "50 25C F50" "18 25C F18" "70 25C F70b" "70 -40C F70cold" "50 25C F50b"; do
    set -- $run
    python3 "$SIM/setup_duplex.py" inputs "$3" dup.json.oxdna "$1" 2000000 2000 "$2" 3 > "setup_$3.log"
    "$OX" "$3_input" > "$3.out" 2>&1
    echo "$3 done $(date +%H:%M:%S)" >> status.txt
done
