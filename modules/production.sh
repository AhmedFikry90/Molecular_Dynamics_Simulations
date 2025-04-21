#!/bin/bash
#
# Production run module for MD simulation workflow
#
# Author: AhmedFikry90

# Run production MD simulation
run_production() {
  log_info "Starting production MD run..."
  
  # Create output directory
  local output_dir=$(create_output_directory "production")
  if [[ $? -ne 0 ]]; then
    log_error "Failed to create output directory"
    return 1
  fi
  
  # Perform MD engine-specific production run
  case "${MD_ENGINE,,}" in
    gromacs)
      run_production_gromacs "$output_dir"
      ;;
    amber)
      run_production_amber "$output_dir"
      ;;
    namd)
      run_production_namd "$output_dir"
      ;;
    *)
      log_error "Unsupported MD engine for production run: $MD_ENGINE"
      return 1
      ;;
  esac
  
  # Return the result of the production run
  local result=$?
  if [[ $result -eq 0 ]]; then
    log_success "Production MD run completed successfully"
  else
    log_error "Production MD run failed"
  fi
  
  return $result
}

# Run production simulation for GROMACS
run_production_gromacs() {
  local output_dir="$1"
  
  log_info "Running production simulation with GROMACS..."
  
  # Find the most recent equilibration files
  local equil_dir=$(find "$OUTPUT_DIR" -type d -name "equilibrate_*" | sort -r | head -n 1)
  
  if [[ -z "$equil_dir" || ! -d "$equil_dir" ]]; then
    log_warning "No equilibration results found. Running equilibration first..."
    equilibrate_system
    equil_dir=$(find "$OUTPUT_DIR" -type d -name "equilibrate_*" | sort -r | head -n 1)
    
    if [[ -z "$equil_dir" || ! -d "$equil_dir" ]]; then
      log_error "Failed to equilibrate system. Cannot continue production run."
      return 1
    }
  fi
  
  # Check for required files
  local struct_file="${equil_dir}/npt.gro"
  local cpt_file="${equil_dir}/npt.cpt"
  local topol_file="${equil_dir}/topol.top"
  
  if [[ ! -f "$struct_file" || ! -f "$cpt_file" || ! -f "$topol_file" ]]; then
    log_error "Required GROMACS equilibration files not found"
    return 1
  fi
  
  # Copy files to production directory
  if [[ "$DRY_RUN" == "true" ]]; then
    log_info "[DRY-RUN] Would copy equilibration files to production directory"
  else
    cp "$struct_file" "${output_dir}/npt.gro"
    cp "$cpt_file" "${output_dir}/npt.cpt"
    cp "$topol_file" "${output_dir}/topol.top"
    
    # Also copy position restraint file if it exists
    if [[ -f "${equil_dir}/posre.itp" ]]; then
      cp "${equil_dir}/posre.itp" "${output_dir}/"
    fi
  fi
  
  # Generate production run configuration
  local md_mdp="${output_dir}/md.mdp"
  generate_config_template "production" "$md_mdp"
  
  # Run production simulation
  if command -v gmx &> /dev/null; then
    # Prepare production run
    local grompp_cmd="gmx grompp -f $md_mdp -c ${output_dir}/npt.gro -t ${output_dir}/npt.cpt -p ${output_dir}/topol.top -o ${output_dir}/md.tpr"
    
    if [[ "$DRY_RUN" == "true" ]]; then
      log_info "[DRY-RUN] Would prepare production run with: $grompp_cmd"
    else
      log_command "$grompp_cmd"
      eval "$grompp_cmd" > "${output_dir}/grompp_md.log" 2>&1 || {
        log_error "Failed to prepare production run with grompp"
        return 1
      }
    fi
    
    # Run production simulation
    local mdrun_cmd="gmx mdrun -v -deffnm ${output_dir}/md -ntomp $NUM_THREADS"
    
    # Create a job script for HPC if using a job scheduler
    if [[ "$JOB_SCHEDULER" != "none" ]]; then
      local job_script="${output_dir}/submit_md.sh"
      create_job_script "md_production" "$mdrun_cmd" "$job_script"
      
      if [[ "$DRY_RUN" == "true" ]]; then
        log_info "[DRY-RUN] Would submit job: $job_script"
      else
        log_info "Submitting production job to $JOB_SCHEDULER"
        submit_job "$job_script"
      fi
    else
      if [[ "$DRY_RUN" == "true" ]]; then
        log_info "[DRY-RUN] Would run production simulation with: $mdrun_cmd"
      else
        log_command "$mdrun_cmd"
        eval "$mdrun_cmd" > "${output_dir}/mdrun_md.log" 2>&1 || {
          log_error "Production simulation failed"
          return 1
        }
      fi
    fi
  else
    log_error "GROMACS (gmx) command not found in PATH"
    return 1
  fi
  
  log_success "GROMACS production simulation complete"
  log_info "Production results in: $output_dir"
  return 0
}

# Run production simulation for AMBER
run_production_amber() {
  local output_dir="$1"
  
  log_info "Running production simulation with AMBER..."
  
  # Find the most recent equilibration files
  local equil_dir=$(find "$OUTPUT_DIR" -type d -name "equilibrate_*" | sort -r | head -n 1)
  
  if [[ -z "$equil_dir" || ! -d "$equil_dir" ]]; then
    log_warning "No equilibration results found. Running equilibration first..."
    equilibrate_system
    equil_dir=$(find "$OUTPUT_DIR" -type d -name "equilibrate_*" | sort -r | head -n 1)
    
    if [[ -z "$equil_dir" || ! -d "$equil_dir" ]]; then
      log_error "Failed to equilibrate system. Cannot continue production run."
      return 1
    }
  fi
  
  # Check for required files
  local prmtop_file="${equil_dir}/system.prmtop"
  local rst_file="${equil_dir}/npt.rst"
  
  if [[ ! -f "$prmtop_file" || ! -f "$rst_file" ]]; then
    log_error "Required AMBER equilibration files not found"
    return 1
  fi
  
  # Copy files to production directory
  if [[ "$DRY_RUN" == "true" ]]; then
    log_info "[DRY-RUN] Would copy equilibration files to production directory"
  else
    cp "$prmtop_file" "${output_dir}/system.prmtop"
    cp "$rst_file" "${output_dir}/npt.rst"
  fi
  
  # Generate production run configuration
  local md_file="${output_dir}/md.in"
  generate_config_template "production" "$md_file"
  
  # Run production simulation
  if command -v pmemd &> /dev/null; then
    # Set up command for CPU or GPU
    local md_cmd="pmemd -O -i ${md_file} -p ${output_dir}/system.prmtop -c ${output_dir}/npt.rst -r ${output_dir}/md.rst -o ${output_dir}/md.out -inf ${output_dir}/md.info -x ${output_dir}/md.nc"
    
    # Use GPU if available
    if [[ "$NUM_GPUS" -gt 0 ]] && command -v pmemd.cuda &> /dev/null; then
      md_cmd="pmemd.cuda -O -i ${md_file} -p ${output_dir}/system.prmtop -c ${output_dir}/npt.rst -r ${output_dir}/md.rst -o ${output_dir}/md.out -inf ${output_dir}/md.info -x ${output_dir}/md.nc"
    fi
    
    # Create a job script for HPC if using a job scheduler
    if [[ "$JOB_SCHEDULER" != "none" ]]; then
      local job_script="${output_dir}/submit_md.sh"
      create_job_script "md_production" "$md_cmd" "$job_script"
      
      if [[ "$DRY_RUN" == "true" ]]; then
        log_info "[DRY-RUN] Would submit job: $job_script"
      else
        log_info "Submitting production job to $JOB_SCHEDULER"
        submit_job "$job_script"
      fi
    else
      if [[ "$DRY_RUN" == "true" ]]; then
        log_info "[DRY-RUN] Would run production simulation with: $md_cmd"
      else
        log_command "$md_cmd"
        eval "$md_cmd" || {
          log_error "Production simulation failed"
          return 1
        }
      fi
    fi
  else
    log_error "AMBER (pmemd) command not found in PATH"
    return 1
  fi
  
  log_success "AMBER production simulation complete"
  log_info "Production results in: $output_dir"
  return 0
}

# Run production simulation for NAMD
run_production_namd() {
  local output_dir="$1"
  
  log_info "Running production simulation with NAMD..."
  
  # Find the most recent equilibration files
  local equil_dir=$(find "$OUTPUT_DIR" -type d -name "equilibrate_*" | sort -r | head -n 1)
  
  if [[ -z "$equil_dir" || ! -d "$equil_dir" ]]; then
    log_warning "No equilibration results found. Running equilibration first..."
    equilibrate_system
    equil_dir=$(find "$OUTPUT_DIR" -type d -name "equilibrate_*" | sort -r | head -n 1)
    
    if [[ -z "$equil_dir" || ! -d "$equil_dir" ]]; then
      log_error "Failed to equilibrate system. Cannot continue production run."
      return 1
    }
  fi
  
  # Check for required files
  local psf_file="${equil_dir}/input.psf"
  local pdb_file="${equil_dir}/input.pdb"
  local npt_output="${equil_dir}/npt.restart.coor"
  local npt_vel="${equil_dir}/npt.restart.vel"
  local npt_xsc="${equil_dir}/npt.restart.xsc"
  
  if [[ ! -f "$psf_file" ]]; then
    log_error "PSF file not found in equilibration directory"
    return 1
  fi
  
  # Fall back to original PDB if restart files don't exist
  if [[ ! -f "$npt_output" || ! -f "$npt_vel" || ! -f "$npt_xsc" ]]; then
    log_warning "NPT restart files not found, using original input files"
    npt_output="$pdb_file"
  fi
  
  # Copy files to production directory
  if [[ "$DRY_RUN" == "true" ]]; then
    log_info "[DRY-RUN] Would copy equilibration files to production directory"
  else
    cp "$psf_file" "${output_dir}/input.psf"
    cp "$pdb_file" "${output_dir}/input.pdb"
    
    if [[ -f "$npt_output" ]]; then
      cp "$npt_output" "${output_dir}/npt.restart.coor"
    fi
    
    if [[ -f "$npt_vel" ]]; then
      cp "$npt_vel" "${output_dir}/npt.restart.vel"
    fi
    
    if [[ -f "$npt_xsc" ]]; then
      cp "$npt_xsc" "${output_dir}/npt.restart.xsc"
    fi
    
    # Also copy parameter files if they exist
    if [[ -d "${equil_dir}/parameters" ]]; then
      mkdir -p "${output_dir}/parameters"
      cp "${equil_dir}/parameters/"* "${output_dir}/parameters/" 2>/dev/null
    fi
  fi
  
  # Generate production run configuration
  local md_file="${output_dir}/md.conf"
  generate_config_template "production" "$md_file"
  
  # Update the paths in the configuration file
  if [[ "$DRY_RUN" != "true" ]]; then
    if [[ -d "${output_dir}/parameters" ]]; then
      # Find parameter files in the directory
      local param_files=("${output_dir}/parameters/"*.prm)
      if [[ ${#param_files[@]} -gt 0 ]]; then
        log_info "Updating parameter file paths in configuration..."
        sed -i "s|par_all36_prot.prm|${param_files[0]}|g" "$md_file"
      fi
    fi
  fi
  
  # Run production simulation
  if command -v namd2 &> /dev/null; then
    local md_cmd="namd2 +p$NUM_THREADS $md_file > ${output_dir}/md.log"
    
    # Create a job script for HPC if using a job scheduler
    if [[ "$JOB_SCHEDULER" != "none" ]]; then
      local job_script="${output_dir}/submit_md.sh"
      create_job_script "md_production" "$md_cmd" "$job_script"
      
      if [[ "$DRY_RUN" == "true" ]]; then
        log_info "[DRY-RUN] Would submit job: $job_script"
      else
        log_info "Submitting production job to $JOB_SCHEDULER"
        submit_job "$job_script"
      fi
    else
      if [[ "$DRY_RUN" == "true" ]]; then
        log_info "[DRY-RUN] Would run production simulation with: $md_cmd"
      else
        log_command "$md_cmd"
        eval "$md_cmd" || {
          log_error "Production simulation failed"
          return 1
        }
      fi
    fi
  else
    log_error "NAMD (namd2) command not found in PATH"
    return 1
  fi
  
  log_success "NAMD production simulation complete"
  log_info "Production results in: $output_dir"
  return 0
}