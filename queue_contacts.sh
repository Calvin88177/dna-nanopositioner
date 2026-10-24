#!/bin/bash
# run the two ideal-guide contact tests one after the other
SIM=/home/claude/nanopositioner/sim
"$SIM/run_contact_ideal.sh" c2 2 70 100000 6 i > "$SIM/contact_c2/driver_i.out" 2>&1
"$SIM/run_contact_ideal.sh" c0 0 70 100000 6 i > "$SIM/contact_c0/driver_i.out" 2>&1
echo "queue done $(date +%H:%M:%S)" >> "$SIM/queue_status.txt"
