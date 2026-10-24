#!/bin/bash
# run_contact_ideal.sh <name> <h> <F_pN> <prod_steps> [<free_bp>] [<tagprefix>]
# Contact test in the ideal-guide limit: every nucleotide further than <free_bp> base pairs
# from the contact is held in x and y at its lattice position, so nothing can slip sideways.
set -e
name=$1; h=$2; F=$3; steps=$4; free=${5:-6}; pre=${6:-i}
OX=/home/claude/nanopositioner/oxDNA/build/bin/oxDNA
SIM=/home/claude/nanopositioner/sim
cd "$SIM/contact_$name"

if [ ! -f relax_mc_last.dat ]; then
cat > relax_mc <<EOF
backend = CPU
backend_precision = double
sim_type = MC
ensemble = NVT
interaction_type = DNA2
salt_concentration = 1.0
use_average_seq = 1
T = 25C
delta_translation = 0.1
delta_rotation = 0.1
verlet_skin = 0.5
steps = 1000
max_backbone_force = 5
max_backbone_force_far = 10
topology = $name.json.top
conf_file = $name.json.oxdna
trajectory_file = relax_mc_traj.dat
lastconf_file = relax_mc_last.dat
energy_file = relax_mc_energy.dat
print_conf_interval = 1000
print_energy_every = 500
time_scale = linear
refresh_vel = 1
restart_step_counter = 1
no_stdout_energy = 1
log_file = relax_mc.log
EOF
"$OX" relax_mc > relax_mc.out 2>&1
echo "relax_mc done $(date +%H:%M:%S)" >> status_$pre.txt
fi

python3 "$SIM/setup_contact.py" inputs $name $h ${pre}relax relax_mc_last.dat 0 2.0 5000 5000 25C ideal $free > setup_$pre.log
cat > ${pre}relax_md <<EOF
backend = CPU
backend_precision = double
sim_type = MD
interaction_type = DNA2
salt_concentration = 1.0
use_average_seq = 1
T = 25C
dt = 0.003
verlet_skin = 0.05
thermostat = brownian
newtonian_steps = 103
diff_coeff = 2.5
steps = 5000
max_backbone_force = 5
max_backbone_force_far = 10
topology = $name.json.top
conf_file = relax_mc_last.dat
trajectory_file = ${pre}relax_md_traj.dat
lastconf_file = ${pre}relax_md_last.dat
energy_file = ${pre}relax_md_energy.dat
print_conf_interval = 5000
print_energy_every = 1000
time_scale = linear
refresh_vel = 1
restart_step_counter = 1
no_stdout_energy = 1
external_forces = 1
external_forces_file = ${pre}relax_forces.txt
log_file = ${pre}relax_md.log
EOF
"$OX" ${pre}relax_md > ${pre}relax_md.out 2>&1
echo "relax_md done $(date +%H:%M:%S)" >> status_$pre.txt

# settle and production keep the lattice references of the MC-relaxed structure
python3 "$SIM/setup_contact.py" inputs $name $h ${pre}settle relax_mc_last.dat $F 2.0 15000 1000 25C ideal $free >> setup_$pre.log
sed -i "s/^conf_file = .*/conf_file = ${pre}relax_md_last.dat/" ${pre}settle_input
"$OX" ${pre}settle_input > ${pre}settle.out 2>&1
echo "settle done $(date +%H:%M:%S)" >> status_$pre.txt

python3 "$SIM/setup_contact.py" inputs $name $h ${pre}prod relax_mc_last.dat $F 2.0 $steps 200 25C ideal $free >> setup_$pre.log
sed -i "s/^conf_file = .*/conf_file = ${pre}settle_last.dat/" ${pre}prod_input
"$OX" ${pre}prod_input > ${pre}prod.out 2>&1
echo "prod done $(date +%H:%M:%S)" >> status_$pre.txt
