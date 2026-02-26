#!/bin/bash
# ==============================================================================
# BASH SCRIPT
#
# PURPOSE & CONTEXT
#   A cluster submission script designed to run a "quench" simulation using
#   the `mpol` executable. It defines the physical and numerical parameters,
#   constructs a PBS/Torque job script dynamically, and handles the post-run
#   data retrieval and cleanup.
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. SIMULATION PARAMETERS
# Define the physical and numerical parameters for the TDVP multi-polaron run.
# ------------------------------------------------------------------------------
npadd=10           # Max number of coherent states (polarons) to add dynamically
delta=0.2          # Qubit tunneling splitting (Delta)
alpha=0.4          # Dimensionless light-matter coupling strength
tmax=8             # Total simulation time
merr=0.0000003     # McLachlan error threshold to trigger basis expansion
dt=0.001           # RK4 integration time step
verbose=1          # Output verbosity flag
p0=0.0000001       # Initial vacuum amplitude assigned to newly added polarons

# ------------------------------------------------------------------------------
# 2. JOB SCRIPT PREPARATION
# ------------------------------------------------------------------------------
newjobfile="quench"
cp jobfile $newjobfile

# ------------------------------------------------------------------------------
# 3. EXECUTION COMMAND CONSTRUCTION
# ------------------------------------------------------------------------------
# Construct the main executable command with all the defined parameters.

echo "./mpol -npini 1 -npadd $npadd -del $delta -al $alpha -wmax 5 -wc 1 -tmax $tmax -dt $dt -merr $merr -p0 $p0 -verbose $verbose" >> $newjobfile

# ------------------------------------------------------------------------------
# 4. POST-PROCESSING & CLEANUP COMMANDS
# ------------------------------------------------------------------------------
# Append a command to move the output data from the compute node's temporary 
# scratch directory ($tpdir) back to the directory where `qsub` was launched.
echo "mv ../MPOL_\$tpdir/data/* \$PBS_O_WORKDIR/data/" >> $newjobfile

# Append a command to safely delete the temporary scratch directory on the node.
echo "rm -r ../MPOL_\$tpdir" >> $newjobfile

# ------------------------------------------------------------------------------
# 5. CLUSTER SUBMISSION
# ------------------------------------------------------------------------------
# Submit the fully constructed job script to the PBS/Torque queue.
qsub $newjobfile
rm $newjobfile

