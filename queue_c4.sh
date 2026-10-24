#!/bin/bash
# the four-helix pad, 2 bp proud, under the same 70 pN and ideal guidance as c2 and c0
SIM=/home/claude/nanopositioner/sim
"$SIM/run_contact_ideal.sh" c4 2 70 100000 6 i > "$SIM/contact_c4/driver_i.out" 2>&1
echo "c4 done $(date +%H:%M:%S)" >> "$SIM/queue_status.txt"
