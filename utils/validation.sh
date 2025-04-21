#!/bin/bash
#
# Input validation utilities for MD simulation workflow
#
# Author: AhmedFikry90

# Validate a file exists and is readable
# Usage: validate_file "file_path" "file_description"
validate_file() {
  local file_path="$1"
  local description="${2:-file}"
  
  if [[ -z "$file_path" ]]; then
    log_error "No $description specified"
    return 1
  fi
  
  if [[ ! -f "$file_path" ]]; then
    log_error "Cannot find $description: $file_path"
    return 1
  fi
  
  if [[ ! -r "$file_path" ]]; then
    log_error "Cannot read $description: $file_path (permission denied)"
    return 1
  fi
  
  return 0
}

# Validate a directory exists and is writable
# Usage: validate_directory "dir_path" "dir_description"
validate_directory() {
  local dir_path="$1"
  local description="${2:-directory}"
  
  if [[ -z "$dir_path" ]]; then
    log_error "No $description specified"
    return 1
  fi
  
  if [[ ! -d "$dir_path" ]]; then
    log_info "Creating $description: $dir_path"
    mkdir -p "$dir_path" || {
      log_error "Failed to create $description: $dir_path"
      return 1
    }
  fi
  
  if [[ ! -w "$dir_path" ]]; then
    log_error "Cannot write to $description: $dir_path (permission denied)"
    return 1
  fi
  
  return 0
}

# Validate a numeric value within a range
# Usage: validate_numeric "value" "min" "max" "description"
validate_numeric() {
  local value="$1"
  local min="$2"
  local max="$3"
  local description="${4:-value}"
  
  if ! [[ "$value" =~ ^[0-9]+([.][0-9]+)?$ ]]; then
    log_error "$description must be a number: $value"
    return 1
  fi
  
  if (( $(echo "$value < $min" | bc -l) )); then
    log_error "$description must be at least $min: $value"
    return 1
  fi
  
  if (( $(echo "$value > $max" | bc -l) )); then
    log_error "$description must be at most $max: $value"
    return 1
  fi
  
  return 0
}

# Validate that a required program is in PATH
# Usage: validate_program "program_name" "package_name"
validate_program() {
  local program="$1"
  local package="${2:-$program}"
  
  if ! command -v "$program" &> /dev/null; then
    log_error "Required program '$program' not found in PATH"
    log_info "Please install $package or ensure it's in your PATH"
    return 1
  fi
  
  return 0
}

# Validate an input structure file based on its format
# Usage: validate_structure_file "file_path"
validate_structure_file() {
  local file_path="$1"
  
  # First check if file exists and is readable
  validate_file "$file_path" "structure file" || return 1
  
  # Check file extension
  local ext="${file_path##*.}"
  case "${ext,,}" in
    pdb)
      # Check for ATOM records
      if ! grep -q "^ATOM" "$file_path"; then
        log_error "Invalid PDB file: No ATOM records found in $file_path"
        return 1
      fi
      ;;
    gro)
      # Check for minimum content
      if [[ $(wc -l < "$file_path") -lt 3 ]]; then
        log_error "Invalid GRO file: Too few lines in $file_path"
        return 1
      fi
      ;;
    xyz)
      # Check for minimum content
      if [[ $(wc -l < "$file_path") -lt 2 ]]; then
        log_error "Invalid XYZ file: Too few lines in $file_path"
        return 1
      fi
      ;;
    mol2)
      # Check for MOLECULE record
      if ! grep -q "@<TRIPOS>MOLECULE" "$file_path"; then
        log_error "Invalid MOL2 file: No @<TRIPOS>MOLECULE record found in $file_path"
        return 1
      fi
      ;;
    *)
      log_warning "Unknown structure file extension: .$ext"
      log_info "Supported formats: .pdb, .gro, .xyz, .mol2"
      ;;
  esac
  
  log_debug "Structure file validated: $file_path"
  return 0
}

# Validate a topology file based on its format
# Usage: validate_topology_file "file_path" "md_engine"
validate_topology_file() {
  local file_path="$1"
  local md_engine="${2:-$MD_ENGINE}"
  
  # First check if file exists and is readable
  validate_file "$file_path" "topology file" || return 1
  
  # Check based on MD engine
  case "${md_engine,,}" in
    gromacs)
      # Check for topology section
      if ! grep -q "^\[ \(moleculetype\|atoms\) \]" "$file_path"; then
        log_warning "Possible invalid GROMACS topology: missing required sections in $file_path"
      fi
      ;;
    amber)
      # Check for unit section
      if ! grep -q -i "^\(ATOM\|BOND\)" "$file_path"; then
        log_warning "Possible invalid AMBER topology: missing required sections in $file_path"
      fi
      ;;
    namd)
      # For NAMD, topology usually comes from PSF files
      if [[ "${file_path##*.}" != "psf" ]]; then
        log_warning "NAMD typically uses PSF files for topology, but got: $file_path"
      fi
      ;;
    *)
      log_warning "Unknown MD engine for topology validation: $md_engine"
      ;;
  esac
  
  log_debug "Topology file validated: $file_path"
  return 0
}

# Validate simulation parameters
# Usage: validate_simulation_parameters
validate_simulation_parameters() {
  local errors=0
  
  # Validate timestep
  if ! validate_numeric "$TIMESTEP" 0.0001 0.01 "Timestep"; then
    errors=$((errors + 1))
    log_warning "Recommended timestep range: 0.0005-0.005 ps"
  fi
  
  # Validate temperature
  if ! validate_numeric "$TEMPERATURE" 0 1000 "Temperature"; then
    errors=$((errors + 1))
  fi
  
  # Validate pressure
  if ! validate_numeric "$PRESSURE" 0 100 "Pressure"; then
    errors=$((errors + 1))
  fi
  
  # Validate steps
  if ! validate_numeric "$MINIMIZATION_STEPS" 100 1000000 "Minimization steps"; then
    errors=$((errors + 1))
  fi
  
  if ! validate_numeric "$EQUILIBRATION_STEPS" 1000 10000000 "Equilibration steps"; then
    errors=$((errors + 1))
  fi
  
  if ! validate_numeric "$PRODUCTION_STEPS" 1000 1000000000 "Production steps"; then
    errors=$((errors + 1))
  fi
  
  # Resource allocation
  if ! validate_numeric "$NUM_THREADS" 1 128 "Number of threads"; then
    errors=$((errors + 1))
  fi
  
  if ! validate_numeric "$NUM_GPUS" 0 16 "Number of GPUs"; then
    errors=$((errors + 1))
  fi
  
  if ! validate_numeric "$MPI_PROCESSES" 1 1024 "Number of MPI processes"; then
    errors=$((errors + 1))
  fi
  
  # Output frequency
  if ! validate_numeric "$ENERGY_OUTPUT_FREQ" 1 1000000 "Energy output frequency"; then
    errors=$((errors + 1))
  fi
  
  if ! validate_numeric "$TRAJECTORY_OUTPUT_FREQ" 1 1000000 "Trajectory output frequency"; then
    errors=$((errors + 1))
  fi
  
  if [[ $errors -gt 0 ]]; then
    log_warning "Found $errors issues with simulation parameters"
    return 1
  fi
  
  log_debug "Simulation parameters validated"
  return 0
}