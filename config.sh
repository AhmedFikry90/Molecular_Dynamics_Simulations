#!/bin/bash
#
# Global configuration for MD simulation workflow
#
# Author: AhmedFikry90

# Version information
VERSION="1.0.0"

# Runtime options
VERBOSE="${VERBOSE:-false}"
DRY_RUN="${DRY_RUN:-false}"
FORCE_OVERWRITE="${FORCE_OVERWRITE:-false}"

# Logging configuration
LOG_LEVEL=3                  # 0=none, 1=error, 2=warning, 3=info, 4=debug
LOG_FILE="${LOG_FILE:-}"     # If empty, log to stdout
LOG_TIMESTAMP=true           # Include timestamp in log messages

# MD engine selection
# Supported: gromacs, amber, namd
MD_ENGINE="${MD_ENGINE:-gromacs}"

# Default paths
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
INPUT_DIR="${INPUT_DIR:-input}"
OUTPUT_DIR="${OUTPUT_DIR:-output}"
TEMP_DIR="${TEMP_DIR:-temp}"
LOG_DIR="${LOG_DIR:-logs}"

# Resource allocation
NUM_THREADS="${NUM_THREADS:-$(nproc 2>/dev/null || echo 1)}"  # CPU threads for OpenMP
NUM_GPUS="${NUM_GPUS:-0}"                                    # GPUs for GPU acceleration
MPI_PROCESSES="${MPI_PROCESSES:-1}"                          # MPI processes

# System configuration
TIMESTEP="${TIMESTEP:-0.002}"                # Simulation timestep (ps)
TEMPERATURE="${TEMPERATURE:-300}"            # Simulation temperature (K)
PRESSURE="${PRESSURE:-1.0}"                  # Simulation pressure (bar)
CUTOFF="${CUTOFF:-1.0}"                      # Non-bonded cutoff (nm)

# Simulation settings
MINIMIZATION_STEPS="${MINIMIZATION_STEPS:-5000}"            # Steps for energy minimization
EQUILIBRATION_STEPS="${EQUILIBRATION_STEPS:-100000}"        # Steps for equilibration
PRODUCTION_STEPS="${PRODUCTION_STEPS:-5000000}"             # Steps for production run

# Output settings
ENERGY_OUTPUT_FREQ="${ENERGY_OUTPUT_FREQ:-1000}"            # Energy output frequency (steps)
TRAJECTORY_OUTPUT_FREQ="${TRAJECTORY_OUTPUT_FREQ:-1000}"    # Trajectory output frequency (steps)

# Job submission settings for HPC environments
JOB_SCHEDULER="${JOB_SCHEDULER:-none}"     # Options: none, slurm, pbs, sge
WALLTIME="${WALLTIME:-24:00:00}"           # Wall time limit (HH:MM:SS)
QUEUE="${QUEUE:-}"                         # Queue/partition name
ACCOUNT="${ACCOUNT:-}"                     # Account/project name

# File extensions for different MD engines
case "$MD_ENGINE" in
  gromacs)
    CONFIG_EXT="mdp"
    TRAJECTORY_EXT="xtc"
    ;;
  amber)
    CONFIG_EXT="in"
    TRAJECTORY_EXT="nc"
    ;;
  namd)
    CONFIG_EXT="conf"
    TRAJECTORY_EXT="dcd"
    ;;
  *)
    CONFIG_EXT="conf"
    TRAJECTORY_EXT="trj"
    ;;
esac

# Create required directories if they don't exist
mkdir -p "$INPUT_DIR" "$OUTPUT_DIR" "$TEMP_DIR" "$LOG_DIR"

# Load local configuration if exists
if [[ -f ".md_workflow_config" ]]; then
  source ".md_workflow_config"
fi