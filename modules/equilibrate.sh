#!/bin/bash
#
# System equilibration module for MD simulation workflow
#
# Author: AhmedFikry90

# Run system equilibration
equilibrate_system() {
  log_info "Starting system equilibration..."
  
  # Create output directory
  local output_dir=$(create_output_directory "equilibrate")
  if [[ $? -ne 0 ]]; then
    log_error "Failed to create output directory"
    return 1
  fi
  
  # Perform MD engine-specific equilibration
  case "${MD_ENGINE,,}" in
    gromacs)
      equilibrate_system_gromacs "$output_dir"
      ;;
    amber)
      equilibrate_system_amber "$output_dir"
      ;;
    namd)
      equilibrate_system_namd "$output_dir"
      ;;
    *)
      log_error "Unsupported MD engine for equilibration: $MD_ENGINE"
      return 1
      ;;
  esac
  
  # Return the result of the equilibration
  local result=$?
  if [[ $result -eq 0 ]]; then
    log_success "System equilibration completed successfully"
  else
    log_error "System equilibration failed"
  fi
  
  return $result
}

# Run equilibration for GROMACS
equilibrate_system_gromacs() {
  local output_dir="$1"
  
  log_info "Running system equilibration with GROMACS..."
  
  # Find the most recent minimization files
  local min_dir=$(find "$OUTPUT_DIR" -type d -name "minimize_*" | sort -r | head -n 1)
  
  if [[ -z "$min_dir" || ! -d "$min_dir" ]]; then
    log_warning "No minimization results found. Running energy minimization first..."
    minimize_energy
    min_dir=$(find "$OUTPUT_DIR" -type d -name "minimize_*" | sort -r | head -n 1)
    
    if [[ -z "$min_dir" || ! -d "$min_dir" ]]; then
      log_error "Failed to minimize energy. Cannot continue equilibration."
      return 1
    }
  fi
  
  # Check for required files
  local struct_file="${min_dir}/em.gro"
  if [[ ! -f "$struct_file" ]]; then
    log_error "No minimized structure file found in minimization directory"
    return 1
  fi
  
  local topol_file="${min_dir}/topol.top"
  if [[ ! -f "$topol_file" ]]; then
    log_error "No topology file found in minimization directory"
    return 1
  fi
  
  # Copy files to equilibration directory
  if [[ "$DRY_RUN" == "true" ]]; then
    log_info "[DRY-RUN] Would copy minimization files to equilibration directory"
  else
    cp "$struct_file" "${output_dir}/em.gro"
    cp "$topol_file" "${output_dir}/topol.top"
    
    # Also copy position restraint file if it exists
    if [[ -f "${min_dir}/posre.itp" ]]; then
      cp "${min_dir}/posre.itp" "${output_dir}/"
    fi
  fi
  
  # Two-step equilibration: NVT followed by NPT
  if command -v gmx &> /dev/null; then
    # Generate NVT equilibration configuration
    local nvt_mdp="${output_dir}/nvt.mdp"
    generate_config_template "equilibrate_nvt" "$nvt_mdp"
    
    # Generate NPT equilibration configuration
    local npt_mdp="${output_dir}/npt.mdp"
    generate_config_template "equilibrate_npt" "$npt_mdp"
    
    # Run NVT equilibration
    if [[ "$DRY_RUN" != "true" ]]; then
      # Prepare NVT run
      local nvt_grompp_cmd="gmx grompp -f $nvt_mdp -c ${output_dir}/em.gro -p ${output_dir}/topol.top -o ${output_dir}/nvt.tpr"
      log_command "$nvt_grompp_cmd"
      eval "$nvt_grompp_cmd" > "${output_dir}/grompp_nvt.log" 2>&1 || {
        log_error "Failed to prepare NVT equilibration with grompp"
        return 1
      }
      
      # Run NVT equilibration
      local nvt_mdrun_cmd="gmx mdrun -v -deffnm ${output_dir}/nvt -ntomp $NUM_THREADS"
      log_command "$nvt_mdrun_cmd"
      eval "$nvt_mdrun_cmd" > "${output_dir}/mdrun_nvt.log" 2>&1 || {
        log_error "NVT equilibration failed"
        return 1
      }
      
      # Prepare NPT run
      local npt_grompp_cmd="gmx grompp -f $npt_mdp -c ${output_dir}/nvt.gro -t ${output_dir}/nvt.cpt -p ${output_dir}/topol.top -o ${output_dir}/npt.tpr"
      log_command "$npt_grompp_cmd"
      eval "$npt_grompp_cmd" > "${output_dir}/grompp_npt.log" 2>&1 || {
        log_error "Failed to prepare NPT equilibration with grompp"
        return 1
      }
      
      # Run NPT equilibration
      local npt_mdrun_cmd="gmx mdrun -v -deffnm ${output_dir}/npt -ntomp $NUM_THREADS"
      log_command "$npt_mdrun_cmd"
      eval "$npt_mdrun_cmd" > "${output_dir}/mdrun_npt.log" 2>&1 || {
        log_error "NPT equilibration failed"
        return 1
      }
    else
      log_info "[DRY-RUN] Would run NVT and NPT equilibration"
    fi
  else
    log_error "GROMACS (gmx) command not found in PATH"
    return 1
  fi
  
  log_success "GROMACS equilibration complete"
  log_info "Equilibration results in: $output_dir"
  return 0
}

# Run equilibration for AMBER
equilibrate_system_amber() {
  local output_dir="$1"
  
  log_info "Running system equilibration with AMBER..."
  
  # Find the most recent minimization files
  local min_dir=$(find "$OUTPUT_DIR" -type d -name "minimize_*" | sort -r | head -n 1)
  
  if [[ -z "$min_dir" || ! -d "$min_dir" ]]; then
    log_warning "No minimization results found. Running energy minimization first..."
    minimize_energy
    min_dir=$(find "$OUTPUT_DIR" -type d -name "minimize_*" | sort -r | head -n 1)
    
    if [[ -z "$min_dir" || ! -d "$min_dir" ]]; then
      log_error "Failed to minimize energy. Cannot continue equilibration."
      return 1
    }
  fi
  
  # Check for required files
  local prmtop_file="${min_dir}/system.prmtop"
  local rst_file="${min_dir}/min.rst"
  
  if [[ ! -f "$prmtop_file" || ! -f "$rst_file" ]]; then
    log_error "Required AMBER files not found in minimization directory"
    return 1
  fi
  
  # Copy files to equilibration directory
  if [[ "$DRY_RUN" == "true" ]]; then
    log_info "[DRY-RUN] Would copy minimization files to equilibration directory"
  else
    cp "$prmtop_file" "${output_dir}/system.prmtop"
    cp "$rst_file" "${output_dir}/min.rst"
  fi
  
  # Two-step equilibration: NVT followed by NPT
  if command -v pmemd &> /dev/null; then
    # Generate NVT equilibration configuration
    local nvt_file="${output_dir}/nvt.in"
    generate_config_template "equilibrate_nvt" "$nvt_file"
    
    # Generate NPT equilibration configuration
    local npt_file="${output_dir}/npt.in"
    generate_config_template "equilibrate_npt" "$npt_file"
    
    # Set up GPU or CPU commands
    local nvt_cmd="pmemd -O -i ${nvt_file} -p ${output_dir}/system.prmtop -c ${output_dir}/min.rst -r ${output_dir}/nvt.rst -o ${output_dir}/nvt.out -inf ${output_dir}/nvt.info -x ${output_dir}/nvt.nc"
    local npt_cmd="pmemd -O -i ${npt_file} -p ${output_dir}/system.prmtop -c ${output_dir}/nvt.rst -r ${output_dir}/npt.rst -o ${output_dir}/npt.out -inf ${output_dir}/npt.info -x ${output_dir}/npt.nc"
    
    # Use GPU if available
    if [[ "$NUM_GPUS" -gt 0 ]] && command -v pmemd.cuda &> /dev/null; then
      nvt_cmd="pmemd.cuda -O -i ${nvt_file} -p ${output_dir}/system.prmtop -c ${output_dir}/min.rst -r ${output_dir}/nvt.rst -o ${output_dir}/nvt.out -inf ${output_dir}/nvt.info -x ${output_dir}/nvt.nc"
      npt_cmd="pmemd.cuda -O -i ${npt_file} -p ${output_dir}/system.prmtop -c ${output_dir}/nvt.rst -r ${output_dir}/npt.rst -o ${output_dir}/npt.out -inf ${output_dir}/npt.info -x ${output_dir}/npt.nc"
    fi
    
    # Run NVT equilibration
    if [[ "$DRY_RUN" == "true" ]]; then
      log_info "[DRY-RUN] Would run NVT equilibration with: $nvt_cmd"
      log_info "[DRY-RUN] Would run NPT equilibration with: $npt_cmd"
    else
      log_command "$nvt_cmd"
      eval "$nvt_cmd" || {
        log_error "NVT equilibration failed"
        return 1
      }
      
      log_command "$npt_cmd"
      eval "$npt_cmd" || {
        log_error "NPT equilibration failed"
        return 1
      }
    fi
  else
    log_error "AMBER (pmemd) command not found in PATH"
    return 1
  fi
  
  log_success "AMBER equilibration complete"
  log_info "Equilibration results in: $output_dir"
  return 0
}

# Run equilibration for NAMD
equilibrate_system_namd() {
  local output_dir="$1"
  
  log_info "Running system equilibration with NAMD..."
  
  # Find the most recent minimization files
  local min_dir=$(find "$OUTPUT_DIR" -type d -name "minimize_*" | sort -r | head -n 1)
  
  if [[ -z "$min_dir" || ! -d "$min_dir" ]]; then
    log_warning "No minimization results found. Running energy minimization first..."
    minimize_energy
    min_dir=$(find "$OUTPUT_DIR" -type d -name "minimize_*" | sort -r | head -n 1)
    
    if [[ -z "$min_dir" || ! -d "$min_dir" ]]; then
      log_error "Failed to minimize energy. Cannot continue equilibration."
      return 1
    }
  fi
  
  # Check for required files
  local psf_file="${min_dir}/input.psf"
  local min_output="${min_dir}/min.coor"
  
  if [[ ! -f "$psf_file" ]]; then
    log_error "PSF file not found in minimization directory"
    return 1
  }
  
  # NAMD might not have produced the expected output file
  if [[ ! -f "$min_output" ]]; then
    min_output="${min_dir}/min.restart.coor"
    if [[ ! -f "$min_output" ]]; then
      # Fall back to the input PDB file
      min_output="${min_dir}/input.pdb"
      if [[ ! -f "$min_output" ]]; then
        log_error "No minimized structure found in minimization directory"
        return 1
      fi
    fi
  fi
  
  # Copy files to equilibration directory
  if [[ "$DRY_RUN" == "true" ]]; then
    log_info "[DRY-RUN] Would copy minimization files to equilibration directory"
  else
    cp "$psf_file" "${output_dir}/input.psf"
    cp "$min_output" "${output_dir}/input.pdb"
    
    # Also copy parameter files if they exist
    if [[ -d "${min_dir}/parameters" ]]; then
      mkdir -p "${output_dir}/parameters"
      cp "${min_dir}/parameters/"* "${output_dir}/parameters/" 2>/dev/null
    fi
  fi
  
  # Two-step equilibration: NVT followed by NPT
  if command -v namd2 &> /dev/null; then
    # Generate NVT equilibration configuration
    local nvt_file="${output_dir}/nvt.conf"
    generate_config_template "equilibrate_nvt" "$nvt_file"
    
    # Generate NPT equilibration configuration
    local npt_file="${output_dir}/npt.conf"
    generate_config_template "equilibrate_npt" "$npt_file"
    
    # Update the paths in the configuration files
    if [[ "$DRY_RUN" != "true" ]]; then
      if [[ -d "${output_dir}/parameters" ]]; then
        # Find parameter files in the directory
        local param_files=("${output_dir}/parameters/"*.prm)
        if [[ ${#param_files[@]} -gt 0 ]]; then
          log_info "Updating parameter file paths in configuration..."
          sed -i "s|par_all36_prot.prm|${param_files[0]}|g" "$nvt_file"
          sed -i "s|par_all36_prot.prm|${param_files[0]}|g" "$npt_file"
        fi
      fi
    fi
    
    # Run NVT equilibration
    if [[ "$DRY_RUN" == "true" ]]; then
      log_info "[DRY-RUN] Would run NVT and NPT equilibration with NAMD"
    else
      local nvt_cmd="namd2 +p$NUM_THREADS $nvt_file > ${output_dir}/nvt.log"
      log_command "$nvt_cmd"
      eval "$nvt_cmd" || {
        log_error "NVT equilibration failed"
        return 1
      }
      
      local npt_cmd="namd2 +p$NUM_THREADS $npt_file > ${output_dir}/npt.log"
      log_command "$npt_cmd"
      eval "$npt_cmd" || {
        log_error "NPT equilibration failed"
        return 1
      }
    fi
  else
    log_error "NAMD (namd2) command not found in PATH"
    return 1
  fi
  
  log_success "NAMD equilibration complete"
  log_info "Equilibration results in: $output_dir"
  return 0
}