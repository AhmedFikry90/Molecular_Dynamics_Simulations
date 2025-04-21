#!/bin/bash
#
# Energy minimization module for MD simulation workflow
#
# Author: AhmedFikry90

# Run energy minimization
minimize_energy() {
  log_info "Starting energy minimization..."
  
  # Create output directory
  local output_dir=$(create_output_directory "minimize")
  if [[ $? -ne 0 ]]; then
    log_error "Failed to create output directory"
    return 1
  fi
  
  # Perform MD engine-specific energy minimization
  case "${MD_ENGINE,,}" in
    gromacs)
      minimize_energy_gromacs "$output_dir"
      ;;
    amber)
      minimize_energy_amber "$output_dir"
      ;;
    namd)
      minimize_energy_namd "$output_dir"
      ;;
    *)
      log_error "Unsupported MD engine for energy minimization: $MD_ENGINE"
      return 1
      ;;
  esac
  
  # Return the result of the minimization
  local result=$?
  if [[ $result -eq 0 ]]; then
    log_success "Energy minimization completed successfully"
  else
    log_error "Energy minimization failed"
  fi
  
  return $result
}

# Run energy minimization for GROMACS
minimize_energy_gromacs() {
  local output_dir="$1"
  
  log_info "Running energy minimization with GROMACS..."
  
  # Find the most recent setup files
  local setup_dir=$(find "$OUTPUT_DIR" -type d -name "setup_*" | sort -r | head -n 1)
  
  if [[ -z "$setup_dir" || ! -d "$setup_dir" ]]; then
    log_warning "No setup files found. Running system setup first..."
    setup_system
    setup_dir=$(find "$OUTPUT_DIR" -type d -name "setup_*" | sort -r | head -n 1)
    
    if [[ -z "$setup_dir" || ! -d "$setup_dir" ]]; then
      log_error "Failed to set up system. Cannot continue minimization."
      return 1
    }
  fi
  
  # Check for required files
  local struct_file="${setup_dir}/solvated_ions.gro"
  if [[ ! -f "$struct_file" ]]; then
    struct_file="${setup_dir}/solvated.gro"
    if [[ ! -f "$struct_file" ]]; then
      struct_file="${setup_dir}/system.gro"
      if [[ ! -f "$struct_file" ]]; then
        log_error "No structure file found in setup directory"
        return 1
      fi
    fi
  fi
  
  local topol_file="${setup_dir}/topol.top"
  if [[ ! -f "$topol_file" ]]; then
    log_error "No topology file found in setup directory"
    return 1
  fi
  
  # Copy files to minimization directory
  if [[ "$DRY_RUN" == "true" ]]; then
    log_info "[DRY-RUN] Would copy setup files to minimization directory"
  else
    cp "$struct_file" "${output_dir}/system.gro"
    cp "$topol_file" "${output_dir}/topol.top"
    
    # Also copy position restraint file if it exists
    if [[ -f "${setup_dir}/posre.itp" ]]; then
      cp "${setup_dir}/posre.itp" "${output_dir}/"
    fi
  fi
  
  # Generate minimization configuration
  local mdp_file="${output_dir}/em.mdp"
  generate_config_template "minimize" "$mdp_file"
  
  # Run energy minimization
  if command -v gmx &> /dev/null; then
    # Generate tpr file
    local grompp_cmd="gmx grompp -f $mdp_file -c ${output_dir}/system.gro -p ${output_dir}/topol.top -o ${output_dir}/em.tpr"
    
    if [[ "$DRY_RUN" == "true" ]]; then
      log_info "[DRY-RUN] Would prepare minimization with: $grompp_cmd"
    else
      log_command "$grompp_cmd"
      eval "$grompp_cmd" > "${output_dir}/grompp.log" 2>&1 || {
        log_error "Failed to prepare minimization with grompp"
        return 1
      }
    fi
    
    # Run minimization
    local mdrun_cmd="gmx mdrun -v -deffnm ${output_dir}/em -ntomp $NUM_THREADS"
    
    if [[ "$DRY_RUN" == "true" ]]; then
      log_info "[DRY-RUN] Would run energy minimization with: $mdrun_cmd"
    else
      log_command "$mdrun_cmd"
      eval "$mdrun_cmd" > "${output_dir}/mdrun_em.log" 2>&1 || {
        log_error "Energy minimization failed"
        return 1
      }
    fi
  else
    log_error "GROMACS (gmx) command not found in PATH"
    return 1
  fi
  
  log_success "GROMACS energy minimization complete"
  log_info "Minimization results in: $output_dir"
  return 0
}

# Run energy minimization for AMBER
minimize_energy_amber() {
  local output_dir="$1"
  
  log_info "Running energy minimization with AMBER..."
  
  # Find the most recent setup files
  local setup_dir=$(find "$OUTPUT_DIR" -type d -name "setup_*" | sort -r | head -n 1)
  
  if [[ -z "$setup_dir" || ! -d "$setup_dir" ]]; then
    log_warning "No setup files found. Running system setup first..."
    setup_system
    setup_dir=$(find "$OUTPUT_DIR" -type d -name "setup_*" | sort -r | head -n 1)
    
    if [[ -z "$setup_dir" || ! -d "$setup_dir" ]]; then
      log_error "Failed to set up system. Cannot continue minimization."
      return 1
    }
  fi
  
  # Check for required files
  local prmtop_file="${setup_dir}/solvated.prmtop"
  local inpcrd_file="${setup_dir}/solvated.inpcrd"
  
  if [[ ! -f "$prmtop_file" || ! -f "$inpcrd_file" ]]; then
    prmtop_file="${setup_dir}/system.prmtop"
    inpcrd_file="${setup_dir}/system.inpcrd"
    
    if [[ ! -f "$prmtop_file" || ! -f "$inpcrd_file" ]]; then
      log_error "Required AMBER files not found in setup directory"
      return 1
    fi
  fi
  
  # Copy files to minimization directory
  if [[ "$DRY_RUN" == "true" ]]; then
    log_info "[DRY-RUN] Would copy setup files to minimization directory"
  else
    cp "$prmtop_file" "${output_dir}/system.prmtop"
    cp "$inpcrd_file" "${output_dir}/system.inpcrd"
  fi
  
  # Generate minimization configuration
  local min_file="${output_dir}/min.in"
  generate_config_template "minimize" "$min_file"
  
  # Run energy minimization
  if command -v pmemd &> /dev/null; then
    local amber_cmd="pmemd -O -i ${min_file} -p ${output_dir}/system.prmtop -c ${output_dir}/system.inpcrd -r ${output_dir}/min.rst -o ${output_dir}/min.out -inf ${output_dir}/min.info"
    
    # Use GPU if available
    if [[ "$NUM_GPUS" -gt 0 ]] && command -v pmemd.cuda &> /dev/null; then
      amber_cmd="pmemd.cuda -O -i ${min_file} -p ${output_dir}/system.prmtop -c ${output_dir}/system.inpcrd -r ${output_dir}/min.rst -o ${output_dir}/min.out -inf ${output_dir}/min.info"
    fi
    
    if [[ "$DRY_RUN" == "true" ]]; then
      log_info "[DRY-RUN] Would run energy minimization with: $amber_cmd"
    else
      log_command "$amber_cmd"
      eval "$amber_cmd" || {
        log_error "Energy minimization failed"
        return 1
      }
    fi
  else
    log_error "AMBER (pmemd) command not found in PATH"
    return 1
  fi
  
  log_success "AMBER energy minimization complete"
  log_info "Minimization results in: $output_dir"
  return 0
}

# Run energy minimization for NAMD
minimize_energy_namd() {
  local output_dir="$1"
  
  log_info "Running energy minimization with NAMD..."
  
  # Find the most recent setup files
  local setup_dir=$(find "$OUTPUT_DIR" -type d -name "setup_*" | sort -r | head -n 1)
  
  if [[ -z "$setup_dir" || ! -d "$setup_dir" ]]; then
    log_warning "No setup files found. Running system setup first..."
    setup_system
    setup_dir=$(find "$OUTPUT_DIR" -type d -name "setup_*" | sort -r | head -n 1)
    
    if [[ -z "$setup_dir" || ! -d "$setup_dir" ]]; then
      log_error "Failed to set up system. Cannot continue minimization."
      return 1
    }
  fi
  
  # Check for required files
  local psf_file="${setup_dir}/ionized.psf"
  local pdb_file="${setup_dir}/ionized.pdb"
  
  if [[ ! -f "$psf_file" || ! -f "$pdb_file" ]]; then
    psf_file="${setup_dir}/solvated.psf"
    pdb_file="${setup_dir}/solvated.pdb"
    
    if [[ ! -f "$psf_file" || ! -f "$pdb_file" ]]; then
      psf_file="${setup_dir}/system.psf"
      pdb_file="${setup_dir}/system.pdb"
      
      if [[ ! -f "$psf_file" || ! -f "$pdb_file" ]]; then
        log_error "Required NAMD files not found in setup directory"
        return 1
      fi
    fi
  fi
  
  # Copy files to minimization directory
  if [[ "$DRY_RUN" == "true" ]]; then
    log_info "[DRY-RUN] Would copy setup files to minimization directory"
  else
    cp "$psf_file" "${output_dir}/input.psf"
    cp "$pdb_file" "${output_dir}/input.pdb"
    
    # Also copy parameter files if they exist
    if [[ -d "${setup_dir}/parameters" ]]; then
      mkdir -p "${output_dir}/parameters"
      cp "${setup_dir}/parameters/"* "${output_dir}/parameters/" 2>/dev/null
    fi
  fi
  
  # Generate minimization configuration
  local min_file="${output_dir}/min.conf"
  generate_config_template "minimize" "$min_file"
  
  # Update the paths in the configuration file
  if [[ "$DRY_RUN" != "true" ]]; then
    if [[ -d "${output_dir}/parameters" ]]; then
      # Find parameter files in the directory
      local param_files=("${output_dir}/parameters/"*.prm)
      if [[ ${#param_files[@]} -gt 0 ]]; then
        log_info "Updating parameter file paths in configuration..."
        sed -i "s|par_all36_prot.prm|${param_files[0]}|g" "$min_file"
      fi
    fi
  fi
  
  # Run energy minimization
  if command -v namd2 &> /dev/null; then
    local namd_cmd="namd2 +p$NUM_THREADS $min_file > ${output_dir}/min.log"
    
    if [[ "$DRY_RUN" == "true" ]]; then
      log_info "[DRY-RUN] Would run energy minimization with: $namd_cmd"
    else
      log_command "$namd_cmd"
      eval "$namd_cmd" || {
        log_error "Energy minimization failed"
        return 1
      }
    fi
  else
    log_error "NAMD (namd2) command not found in PATH"
    return 1
  fi
  
  log_success "NAMD energy minimization complete"
  log_info "Minimization results in: $output_dir"
  return 0
}