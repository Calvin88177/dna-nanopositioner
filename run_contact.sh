#!/bin/bash
# run_contact.sh <name> <h> <F_pN> <prod_steps>
# relax -> settle under weak guidance with the preload on -> production with stiff guidance
set -e
name=$1; h=$2; F=$3; steps=$4
OX=/home/claude/nanopositioner/oxDNA/build/bin/oxDNA
SIM=/home/claude/nanopositioner/sim
cd "$SIM/contact_$name"

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

cat > relax_md <<EOF
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
steps = 10000
max_backbone_force = 5
max_backbone_force_far = 10
topology = $name.json.top
conf_file = relax_mc_last.dat
trajectory_file = relax_md_traj.dat
lastconf_file = relax_md_last.dat
energy_file = relax_md_energy.dat
print_conf_interval = 10000
print_energy_every = 2000
time_scale = linear
refresh_vel = 1
restart_step_counter = 1
no_stdout_energy = 1
external_forces = 1
external_forces_file = relaxmd_forces.txt
log_file = relax_md.log
EOF

"$OX" relax_mc > relax_mc.out 2>&1
echo "relax_mc done $(date +%H:%M:%S)" >> status.txt
# weak guidance during MD relaxation and settling, references from the MC-relaxed structure
python3 "$SIM/setup_contact.py" inputs $name $h relaxmd relax_mc_last.dat 0 0.05 10000 > setup.log
"$OX" relax_md > relax_md.out 2>&1
echo "relax_md done $(date +%H:%M:%S)" >> status.txt
python3 "$SIM/setup_contact.py" inputs $name $h settle relax_md_last.dat $F 0.05 40000 2000 >> setup.log
"$OX" settle_input > settle.out 2>&1
echo "settle done $(date +%H:%M:%S)" >> status.txt
python3 "$SIM/setup_contact.py" inputs $name $h prod settle_last.dat $F 2.0 $steps 200 >> setup.log
"$OX" prod_input > prod.out 2>&1
echo "prod done $(date +%H:%M:%S)" >> status.txt
