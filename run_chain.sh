#!/bin/bash
# run_chain.sh <input1> <input2> ...   - run oxDNA inputs one after another in the current directory
OX=/home/claude/nanopositioner/oxDNA/build/bin/oxDNA
for inp in "$@"; do
    "$OX" "$inp" > "${inp}.out" 2>&1 || { echo "FAILED $inp" >> chain_status.txt; exit 1; }
    echo "finished $inp $(date +%H:%M:%S)" >> chain_status.txt
done
echo "all done" >> chain_status.txt
