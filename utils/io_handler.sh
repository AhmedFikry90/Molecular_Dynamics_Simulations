#!/bin/bash
#
# I/O handling utilities for MD simulation workflow
#
# Author: AhmedFikry90

# Create a timestamped backup of a file
# Usage: backup_file "file_path"
backup_file() {
  local file_path="$1"
  
  if [[ ! -f "$file_path" ]]; then
    log_debug "No file to backup: $file_path"
    return 0
  fi
  
  local backup_path="${file_path}.bak.$(date +%Y%m%d%H%M%S)"
  log_debug "Creating backup: $backup_path"
  
  if [[ "$DRY_RUN" == "true" ]]; then
    log_info "[DRY-RUN] Would backup $file_path to $backup_path"
    return 0
  fi
  
  cp "$file_path" "$backup_path" || {
    log_error "Failed to create backup of $file_path"
    return 1
  }
  
  return 0
}

# Create a unique output directory for a simulation stage
# Usage: create_output_directory "stage_name"
create_output_directory() {
  local stage_name="$1"
  local timestamp=$(date +%Y%m%d%H%M%S)
  local output_dir="${OUTPUT_DIR}/${stage_name}_${timestamp}"
  
  log_debug "Creating output directory: $output_dir"
  
  if [[ "$DRY_RUN" == "true" ]]; then
    log_info "[DRY-RUN] Would create directory: $output_dir"
    echo "$output_dir"
    return 0
  fi
  
  mkdir -p "$output_dir" || {
    log_error "Failed to create output directory: $output_dir"
    return 1
  }
  
  echo "$output_dir"
  return 0
}

# Save configuration data to a file
# Usage: save_config "file_path" "config_data"
save_config() {
  local file_path="$1"
  local config_data="$2"
  
  log_debug "Saving configuration to: $file_path"
  
  if [[ "$DRY_RUN" == "true" ]]; then
    log_info "[DRY-RUN] Would save configuration to: $file_path"
    return 0
  fi
  
  echo "$config_data" > "$file_path" || {
    log_error "Failed to save configuration to: $file_path"
    return 1
  }
  
  return 0
}

# Generate a template configuration file based on MD engine
# Usage: generate_config_template "stage" "output_path"
generate_config_template() {
  local stage="$1"
  local output_path="$2"
  
  log_debug "Generating $MD_ENGINE configuration template for $stage: $output_path"
  
  if [[ "$DRY_RUN" == "true" ]]; then
    log_info "[DRY-RUN] Would generate configuration template: $output_path"
    return 0
  fi
  
  local template_content=""
  
  case "${MD_ENGINE,,}" in
    gromacs)
      case "${stage,,}" in
        minimize)
          template_content="; GROMACS energy minimization mdp file
integrator      = steep
emtol           = 1000.0
emstep          = 0.01
nsteps          = $MINIMIZATION_STEPS

cutoff-scheme   = Verlet
nstlist         = 1
rlist           = 1.0
coulombtype     = PME
rcoulomb        = 1.0
vdwtype         = Cut-off
rvdw            = 1.0
pbc             = xyz"
          ;;
        equilibrate_nvt)
          template_content="; GROMACS NVT equilibration mdp file
integrator      = md
dt              = $TIMESTEP
nsteps          = $EQUILIBRATION_STEPS

cutoff-scheme   = Verlet
nstlist         = 10
rlist           = 1.0
coulombtype     = PME
rcoulomb        = 1.0
vdwtype         = Cut-off
rvdw            = 1.0
pbc             = xyz

tcoupl          = V-rescale
tc-grps         = System
tau_t           = 0.1
ref_t           = $TEMPERATURE

constraints     = h-bonds
constraint_algorithm = LINCS

nstxout         = $TRAJECTORY_OUTPUT_FREQ
nstvout         = $TRAJECTORY_OUTPUT_FREQ
nstenergy       = $ENERGY_OUTPUT_FREQ
nstlog          = $ENERGY_OUTPUT_FREQ"
          ;;
        equilibrate_npt)
          template_content="; GROMACS NPT equilibration mdp file
integrator      = md
dt              = $TIMESTEP
nsteps          = $EQUILIBRATION_STEPS

cutoff-scheme   = Verlet
nstlist         = 10
rlist           = 1.0
coulombtype     = PME
rcoulomb        = 1.0
vdwtype         = Cut-off
rvdw            = 1.0
pbc             = xyz

tcoupl          = V-rescale
tc-grps         = System
tau_t           = 0.1
ref_t           = $TEMPERATURE

pcoupl          = Parrinello-Rahman
pcoupltype      = isotropic
tau_p           = 2.0
ref_p           = $PRESSURE
compressibility = 4.5e-5

constraints     = h-bonds
constraint_algorithm = LINCS

nstxout         = $TRAJECTORY_OUTPUT_FREQ
nstvout         = $TRAJECTORY_OUTPUT_FREQ
nstenergy       = $ENERGY_OUTPUT_FREQ
nstlog          = $ENERGY_OUTPUT_FREQ"
          ;;
        production)
          template_content="; GROMACS production run mdp file
integrator      = md
dt              = $TIMESTEP
nsteps          = $PRODUCTION_STEPS

cutoff-scheme   = Verlet
nstlist         = 10
rlist           = 1.0
coulombtype     = PME
rcoulomb        = 1.0
vdwtype         = Cut-off
rvdw            = 1.0
pbc             = xyz

tcoupl          = V-rescale
tc-grps         = System
tau_t           = 0.1
ref_t           = $TEMPERATURE

pcoupl          = Parrinello-Rahman
pcoupltype      = isotropic
tau_p           = 2.0
ref_p           = $PRESSURE
compressibility = 4.5e-5

constraints     = h-bonds
constraint_algorithm = LINCS

nstxout         = $TRAJECTORY_OUTPUT_FREQ
nstvout         = $TRAJECTORY_OUTPUT_FREQ
nstenergy       = $ENERGY_OUTPUT_FREQ
nstlog          = $ENERGY_OUTPUT_FREQ"
          ;;
        *)
          log_error "Unknown stage for GROMACS template: $stage"
          return 1
          ;;
      esac
      ;;
      
    amber)
      case "${stage,,}" in
        minimize)
          template_content="Minimization
&cntrl
  imin=1, maxcyc=$MINIMIZATION_STEPS, 
  ntb=1, 
  cut=$CUTOFF, 
  ntpr=$ENERGY_OUTPUT_FREQ,
  ntxo=1, 
  ntwr=$ENERGY_OUTPUT_FREQ
/"
          ;;
        equilibrate_nvt)
          template_content="NVT equilibration
&cntrl
  imin=0, 
  ntx=1, 
  irest=0, 
  nstlim=$EQUILIBRATION_STEPS, 
  dt=$TIMESTEP, 
  ntf=2, 
  ntc=2, 
  temp0=$TEMPERATURE, 
  ntpr=$ENERGY_OUTPUT_FREQ, 
  ntwx=$TRAJECTORY_OUTPUT_FREQ, 
  cut=$CUTOFF, 
  ntb=1, 
  ntt=3, 
  gamma_ln=2.0, 
  ig=-1
/"
          ;;
        equilibrate_npt)
          template_content="NPT equilibration
&cntrl
  imin=0, 
  ntx=5, 
  irest=1, 
  nstlim=$EQUILIBRATION_STEPS, 
  dt=$TIMESTEP, 
  ntf=2, 
  ntc=2, 
  temp0=$TEMPERATURE, 
  ntpr=$ENERGY_OUTPUT_FREQ, 
  ntwx=$TRAJECTORY_OUTPUT_FREQ, 
  cut=$CUTOFF, 
  ntb=2, 
  ntp=1, 
  pres0=$PRESSURE, 
  ntt=3, 
  gamma_ln=2.0, 
  ig=-1
/"
          ;;
        production)
          template_content="Production
&cntrl
  imin=0, 
  ntx=5, 
  irest=1, 
  nstlim=$PRODUCTION_STEPS, 
  dt=$TIMESTEP, 
  ntf=2, 
  ntc=2, 
  temp0=$TEMPERATURE, 
  ntpr=$ENERGY_OUTPUT_FREQ, 
  ntwx=$TRAJECTORY_OUTPUT_FREQ, 
  cut=$CUTOFF, 
  ntb=2, 
  ntp=1, 
  pres0=$PRESSURE, 
  ntt=3, 
  gamma_ln=2.0, 
  ig=-1
/"
          ;;
        *)
          log_error "Unknown stage for AMBER template: $stage"
          return 1
          ;;
      esac
      ;;
      
    namd)
      case "${stage,,}" in
        minimize)
          template_content="# NAMD configuration file for minimization
# Generated by MD Workflow

structure       input.psf
coordinates     input.pdb
outputname      min

# Force field parameters
paraTypeCharmm  on
parameters      par_all36_prot.prm

# Basic dynamics settings
exclude         scaled1-4
1-4scaling      1.0
switching       on
switchdist      8.0
cutoff          12.0
pairlistdist    13.5
timestep        1.0
nonbondedFreq   1
fullElectFrequency 2

# Periodic Boundary Conditions
cellBasisVector1    50.0    0.0    0.0
cellBasisVector2     0.0   50.0    0.0
cellBasisVector3     0.0    0.0   50.0
cellOrigin           0.0    0.0    0.0
wrapAll             on

# Output
outputEnergies      $ENERGY_OUTPUT_FREQ
dcdfreq             $TRAJECTORY_OUTPUT_FREQ
xstFreq             $TRAJECTORY_OUTPUT_FREQ

# Minimization
minimize            $MINIMIZATION_STEPS
"
          ;;
        equilibrate_nvt)
          template_content="# NAMD configuration file for NVT equilibration
# Generated by MD Workflow

structure       input.psf
coordinates     input.pdb
outputname      nvt

# Force field parameters
paraTypeCharmm  on
parameters      par_all36_prot.prm

# Basic dynamics settings
exclude         scaled1-4
1-4scaling      1.0
switching       on
switchdist      8.0
cutoff          12.0
pairlistdist    13.5
timestep        $TIMESTEP
nonbondedFreq   1
fullElectFrequency 2

# Periodic Boundary Conditions
cellBasisVector1    50.0    0.0    0.0
cellBasisVector2     0.0   50.0    0.0
cellBasisVector3     0.0    0.0   50.0
cellOrigin           0.0    0.0    0.0
wrapAll             on

# Output
outputEnergies      $ENERGY_OUTPUT_FREQ
dcdfreq             $TRAJECTORY_OUTPUT_FREQ
xstFreq             $TRAJECTORY_OUTPUT_FREQ

# Temperature control
langevin            on
langevinDamping     1.0
langevinTemp        $TEMPERATURE
langevinHydrogen    off

# Constraints
rigidBonds          all
rigidTolerance      0.00001
rigidIterations     100

# Run
run                 $EQUILIBRATION_STEPS
"
          ;;
        equilibrate_npt)
          template_content="# NAMD configuration file for NPT equilibration
# Generated by MD Workflow

structure       input.psf
coordinates     input.pdb
velocities      nvt.restart.vel
extendedSystem  nvt.restart.xsc
outputname      npt

# Force field parameters
paraTypeCharmm  on
parameters      par_all36_prot.prm

# Basic dynamics settings
exclude         scaled1-4
1-4scaling      1.0
switching       on
switchdist      8.0
cutoff          12.0
pairlistdist    13.5
timestep        $TIMESTEP
nonbondedFreq   1
fullElectFrequency 2

# Periodic Boundary Conditions
wrapAll             on

# Output
outputEnergies      $ENERGY_OUTPUT_FREQ
dcdfreq             $TRAJECTORY_OUTPUT_FREQ
xstFreq             $TRAJECTORY_OUTPUT_FREQ
restartfreq         $ENERGY_OUTPUT_FREQ

# Temperature control
langevin            on
langevinDamping     1.0
langevinTemp        $TEMPERATURE
langevinHydrogen    off

# Pressure control
useGroupPressure    yes
useFlexibleCell     no
useConstantArea     no
langevinPiston      on
langevinPistonTarget $PRESSURE
langevinPistonPeriod 100.0
langevinPistonDecay  50.0
langevinPistonTemp   $TEMPERATURE

# Constraints
rigidBonds          all
rigidTolerance      0.00001
rigidIterations     100

# Run
run                 $EQUILIBRATION_STEPS
"
          ;;
        production)
          template_content="# NAMD configuration file for production
# Generated by MD Workflow

structure       input.psf
coordinates     input.pdb
velocities      npt.restart.vel
extendedSystem  npt.restart.xsc
outputname      prod

# Force field parameters
paraTypeCharmm  on
parameters      par_all36_prot.prm

# Basic dynamics settings
exclude         scaled1-4
1-4scaling      1.0
switching       on
switchdist      8.0
cutoff          12.0
pairlistdist    13.5
timestep        $TIMESTEP
nonbondedFreq   1
fullElectFrequency 2

# Periodic Boundary Conditions
wrapAll             on

# Output
outputEnergies      $ENERGY_OUTPUT_FREQ
dcdfreq             $TRAJECTORY_OUTPUT_FREQ
xstFreq             $TRAJECTORY_OUTPUT_FREQ
restartfreq         $ENERGY_OUTPUT_FREQ

# Temperature control
langevin            on
langevinDamping     1.0
langevinTemp        $TEMPERATURE
langevinHydrogen    off

# Pressure control
useGroupPressure    yes
useFlexibleCell     no
useConstantArea     no
langevinPiston      on
langevinPistonTarget $PRESSURE
langevinPistonPeriod 100.0
langevinPistonDecay  50.0
langevinPistonTemp   $TEMPERATURE

# Constraints
rigidBonds          all
rigidTolerance      0.00001
rigidIterations     100

# Run
run                 $PRODUCTION_STEPS
"
          ;;
        *)
          log_error "Unknown stage for NAMD template: $stage"
          return 1
          ;;
      esac
      ;;
      
    *)
      log_error "Unsupported MD engine for configuration template: $MD_ENGINE"
      return 1
      ;;
  esac
  
  # Save the template
  save_config "$output_path" "$template_content"
  
  return 0
}

# Create a job submission script for HPC environments
# Usage: create_job_script "job_name" "command" "output_script"
create_job_script() {
  local job_name="$1"
  local command="$2"
  local output_script="$3"
  
  log_debug "Creating job script: $output_script"
  
  if [[ "$DRY_RUN" == "true" ]]; then
    log_info "[DRY-RUN] Would create job script: $output_script"
    return 0
  fi
  
  local script_content=""
  
  case "${JOB_SCHEDULER,,}" in
    slurm)
      script_content="#!/bin/bash
#SBATCH --job-name=${job_name}
#SBATCH --output=${job_name}.%j.out
#SBATCH --error=${job_name}.%j.err
#SBATCH --time=${WALLTIME}
#SBATCH --ntasks=${MPI_PROCESSES}
#SBATCH --cpus-per-task=${NUM_THREADS}
"

      if [[ -n "$QUEUE" ]]; then
        script_content+="#SBATCH --partition=${QUEUE}
"
      fi
      
      if [[ -n "$ACCOUNT" ]]; then
        script_content+="#SBATCH --account=${ACCOUNT}
"
      fi
      
      if [[ "$NUM_GPUS" -gt 0 ]]; then
        script_content+="#SBATCH --gres=gpu:${NUM_GPUS}
"
      fi
      
      script_content+="
# Load modules and set environment variables here

# Run the command
${command}
"
      ;;
      
    pbs)
      script_content="#!/bin/bash
#PBS -N ${job_name}
#PBS -o ${job_name}.out
#PBS -e ${job_name}.err
#PBS -l walltime=${WALLTIME}
#PBS -l select=1:ncpus=${NUM_THREADS}:mpiprocs=${MPI_PROCESSES}
"

      if [[ -n "$QUEUE" ]]; then
        script_content+="#PBS -q ${QUEUE}
"
      fi
      
      if [[ -n "$ACCOUNT" ]]; then
        script_content+="#PBS -A ${ACCOUNT}
"
      fi
      
      if [[ "$NUM_GPUS" -gt 0 ]]; then
        script_content+="#PBS -l ngpus=${NUM_GPUS}
"
      fi
      
      script_content+="
# Load modules and set environment variables here

# Change to the submission directory
cd \$PBS_O_WORKDIR

# Run the command
${command}
"
      ;;
      
    sge)
      script_content="#!/bin/bash
#$ -N ${job_name}
#$ -o ${job_name}.out
#$ -e ${job_name}.err
#$ -l h_rt=${WALLTIME}
#$ -pe mpi ${MPI_PROCESSES}
"

      if [[ -n "$QUEUE" ]]; then
        script_content+="#$ -q ${QUEUE}
"
      fi
      
      if [[ -n "$ACCOUNT" ]]; then
        script_content+="#$ -P ${ACCOUNT}
"
      fi
      
      script_content+="
# Load modules and set environment variables here

# Run the command
${command}
"
      ;;
      
    none)
      script_content="#!/bin/bash
# Direct execution script for ${job_name}
# Created: $(date)

${command}
"
      ;;
      
    *)
      log_error "Unsupported job scheduler: $JOB_SCHEDULER"
      return 1
      ;;
  esac
  
  # Save the script
  save_config "$output_script" "$script_content"
  
  # Make the script executable
  chmod +x "$output_script" || {
    log_error "Failed to make job script executable: $output_script"
    return 1
  }
  
  return 0
}

# Submit a job to the HPC scheduler
# Usage: submit_job "script_path"
submit_job() {
  local script_path="$1"
  
  log_info "Submitting job: $script_path"
  
  if [[ "$DRY_RUN" == "true" ]]; then
    log_info "[DRY-RUN] Would submit job: $script_path"
    return 0
  fi
  
  if [[ ! -x "$script_path" ]]; then
    log_error "Job script is not executable: $script_path"
    return 1
  fi
  
  local job_id=""
  
  case "${JOB_SCHEDULER,,}" in
    slurm)
      job_id=$(sbatch "$script_path" | grep -o -E '[0-9]+')
      ;;
    pbs)
      job_id=$(qsub "$script_path")
      ;;
    sge)
      job_id=$(qsub "$script_path" | grep -o -E '[0-9]+')
      ;;
    none)
      # Direct execution
      log_info "Executing script directly: $script_path"
      "$script_path" &
      job_id=$!
      ;;
    *)
      log_error "Unsupported job scheduler: $JOB_SCHEDULER"
      return 1
      ;;
  esac
  
  if [[ -z "$job_id" ]]; then
    log_error "Failed to submit job: $script_path"
    return 1
  fi
  
  log_success "Job submitted with ID: $job_id"
  echo "$job_id"
  return 0
}