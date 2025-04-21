#!/bin/bash
#
# System setup module for MD simulation workflow
#
# Author: AhmedFikry90

# Set up the simulation system
setup_system() {
  log_info "Starting simulation system setup..."
  
  # Create output directory
  local output_dir=$(create_output_directory "setup")
  if [[ $? -ne 0 ]]; then
    log_error "Failed to create output directory"
    return 1
  fi
  
  # Perform MD engine-specific system setup
  case "${MD_ENGINE,,}" in
    gromacs)
      setup_system_gromacs "$output_dir"
      ;;
    amber)
      setup_system_amber "$output_dir"
      ;;
    namd)
      setup_system_namd "$output_dir"
      ;;
    *)
      log_error "Unsupported MD engine for system setup: $MD_ENGINE"
      return 1
      ;;
  esac
  
  # Return the result of the setup
  local result=$?
  if [[ $result -eq 0 ]]; then
    log_success "System setup completed successfully"
  else
    log_error "System setup failed"
  fi
  
  return $result
}

# Set up system for GROMACS simulations
setup_system_gromacs() {
  local output_dir="$1"
  
  log_info "Setting up system for GROMACS..."
  
  # Find the most recent parameterized files
  local param_dir=$(find "$OUTPUT_DIR" -type d -name "parameterize_*" | sort -r | head -n 1)
  
  if [[ -z "$param_dir" || ! -d "$param_dir" ]]; then
    log_warning "No parameterized files found. Running parameterization first..."
    parameterize_structure
    param_dir=$(find "$OUTPUT_DIR" -type d -name "parameterize_*" | sort -r | head -n 1)
    
    if [[ -z "$param_dir" || ! -d "$param_dir" ]]; then
      log_error "Failed to parameterize structure. Cannot continue setup."
      return 1
    }
  fi
  
  # Check for required files
  local struct_file="${param_dir}/system.gro"
  local topol_file="${param_dir}/topol.top"
  
  if [[ ! -f "$struct_file" || ! -f "$topol_file" ]]; then
    log_error "Required GROMACS files not found in parameterized directory"
    return 1
  fi
  
  # Copy files to setup directory
  if [[ "$DRY_RUN" == "true" ]]; then
    log_info "[DRY-RUN] Would copy parameterized files to setup directory"
  else
    cp "$struct_file" "${output_dir}/system.gro"
    cp "$topol_file" "${output_dir}/topol.top"
    
    # Also copy position restraint file if it exists
    if [[ -f "${param_dir}/posre.itp" ]]; then
      cp "${param_dir}/posre.itp" "${output_dir}/"
    fi
  fi
  
  # Add solvent to the system
  if command -v gmx &> /dev/null; then
    local solvate_cmd="gmx solvate -cp ${output_dir}/system.gro -cs spc216.gro -o ${output_dir}/solvated.gro -p ${output_dir}/topol.top"
    
    if [[ "$DRY_RUN" == "true" ]]; then
      log_info "[DRY-RUN] Would solvate the system"
    else
      log_command "$solvate_cmd"
      eval "$solvate_cmd" > "${output_dir}/solvate.log" 2>&1 || {
        log_error "Failed to solvate the system"
        return 1
      }
    fi
    
    # Add ions to neutralize the system
    if [[ "$DRY_RUN" != "true" ]]; then
      # Generate a minimal mdp file for genion
      local mdp_file="${output_dir}/ions.mdp"
      echo "integrator = steep" > "$mdp_file"
      echo "nsteps = 0" >> "$mdp_file"
      echo "emtol = 1000.0" >> "$mdp_file"
      
      # Generate a tpr file for genion
      local grompp_cmd="gmx grompp -f $mdp_file -c ${output_dir}/solvated.gro -p ${output_dir}/topol.top -o ${output_dir}/ions.tpr"
      log_command "$grompp_cmd"
      eval "$grompp_cmd" > "${output_dir}/grompp_ions.log" 2>&1 || {
        log_error "Failed to generate tpr file for adding ions"
        return 1
      }
      
      # Add ions
      local genion_cmd="gmx genion -s ${output_dir}/ions.tpr -o ${output_dir}/solvated_ions.gro -p ${output_dir}/topol.top -neutral"
      log_command "$genion_cmd"
      echo "SOL" | eval "$genion_cmd" > "${output_dir}/genion.log" 2>&1 || {
        log_warning "Failed to add ions. Continuing with solvated system only."
        cp "${output_dir}/solvated.gro" "${output_dir}/solvated_ions.gro"
      }
    else
      log_info "[DRY-RUN] Would add ions to neutralize the system"
    fi
  else
    log_error "GROMACS (gmx) command not found in PATH"
    return 1
  fi
  
  log_success "GROMACS system setup complete"
  log_info "Setup files in: $output_dir"
  return 0
}

# Set up system for AMBER simulations
setup_system_amber() {
  local output_dir="$1"
  
  log_info "Setting up system for AMBER..."
  
  # Find the most recent parameterized files
  local param_dir=$(find "$OUTPUT_DIR" -type d -name "parameterize_*" | sort -r | head -n 1)
  
  if [[ -z "$param_dir" || ! -d "$param_dir" ]]; then
    log_warning "No parameterized files found. Running parameterization first..."
    parameterize_structure
    param_dir=$(find "$OUTPUT_DIR" -type d -name "parameterize_*" | sort -r | head -n 1)
    
    if [[ -z "$param_dir" || ! -d "$param_dir" ]]; then
      log_error "Failed to parameterize structure. Cannot continue setup."
      return 1
    }
  fi
  
  # Check for required files
  local prmtop_file="${param_dir}/system.prmtop"
  local inpcrd_file="${param_dir}/system.inpcrd"
  
  if [[ ! -f "$prmtop_file" || ! -f "$inpcrd_file" ]]; then
    log_error "Required AMBER files not found in parameterized directory"
    return 1
  fi
  
  # Copy files to setup directory
  if [[ "$DRY_RUN" == "true" ]]; then
    log_info "[DRY-RUN] Would copy parameterized files to setup directory"
  else
    cp "$prmtop_file" "${output_dir}/system.prmtop"
    cp "$inpcrd_file" "${output_dir}/system.inpcrd"
  fi
  
  # Use tleap to add solvent and ions
  if command -v tleap &> /dev/null; then
    # Create a temporary PDB from the inpcrd file
    if command -v ambpdb &> /dev/null && [[ "$DRY_RUN" != "true" ]]; then
      local ambpdb_cmd="ambpdb -p ${output_dir}/system.prmtop -c ${output_dir}/system.inpcrd > ${output_dir}/system.pdb"
      log_command "$ambpdb_cmd"
      eval "$ambpdb_cmd" || {
        log_warning "Failed to convert AMBER structure to PDB. Continuing anyway."
      }
    fi
    
    # Create tleap input for solvation
    local tleap_input="${TEMP_DIR}/tleap_solvate.in"
    cat > "$tleap_input" << EOF
source leaprc.protein.ff14SB
source leaprc.water.tip3p
mol = loadprmtop ${output_dir}/system.prmtop
loadamberparm mol ${output_dir}/system.prmtop ${output_dir}/system.inpcrd
solvatebox mol TIP3PBOX 10.0
addions mol Na+ 0
addions mol Cl- 0
saveamberparm mol ${output_dir}/solvated.prmtop ${output_dir}/solvated.inpcrd
savepdb mol ${output_dir}/solvated.pdb
quit
EOF
    
    local tleap_cmd="tleap -f $tleap_input"
    
    if [[ "$DRY_RUN" == "true" ]]; then
      log_info "[DRY-RUN] Would solvate the system using tleap"
    else
      log_command "$tleap_cmd"
      eval "$tleap_cmd" > "${output_dir}/tleap_solvate.log" 2>&1 || {
        log_error "Failed to solvate the system using tleap"
        return 1
      }
    fi
  else
    log_error "AMBER tleap command not found in PATH"
    return 1
  fi
  
  log_success "AMBER system setup complete"
  log_info "Setup files in: $output_dir"
  return 0
}

# Set up system for NAMD simulations
setup_system_namd() {
  local output_dir="$1"
  
  log_info "Setting up system for NAMD..."
  
  # Find the most recent parameterized files
  local param_dir=$(find "$OUTPUT_DIR" -type d -name "parameterize_*" | sort -r | head -n 1)
  
  if [[ -z "$param_dir" || ! -d "$param_dir" ]]; then
    log_warning "No parameterized files found. Running parameterization first..."
    parameterize_structure
    param_dir=$(find "$OUTPUT_DIR" -type d -name "parameterize_*" | sort -r | head -n 1)
    
    if [[ -z "$param_dir" || ! -d "$param_dir" ]]; then
      log_error "Failed to parameterize structure. Cannot continue setup."
      return 1
    }
  fi
  
  # Check for required files
  local psf_file="${param_dir}/system.psf"
  local pdb_file="${param_dir}/system.pdb"
  
  if [[ ! -f "$psf_file" || ! -f "$pdb_file" ]]; then
    log_error "Required NAMD files not found in parameterized directory"
    return 1
  fi
  
  # Copy files to setup directory
  if [[ "$DRY_RUN" == "true" ]]; then
    log_info "[DRY-RUN] Would copy parameterized files to setup directory"
  else
    cp "$psf_file" "${output_dir}/system.psf"
    cp "$pdb_file" "${output_dir}/system.pdb"
    
    # Also copy parameter files if they exist
    if [[ -d "${param_dir}/parameters" ]]; then
      mkdir -p "${output_dir}/parameters"
      cp "${param_dir}/parameters/"* "${output_dir}/parameters/" 2>/dev/null
    fi
  fi
  
  # Use VMD to solvate the system if available
  if command -v vmd &> /dev/null; then
    # Create VMD script for solvation
    local vmd_script="${TEMP_DIR}/solvate.tcl"
    cat > "$vmd_script" << EOF
package require solvate
package require autoionize

# Load the molecule
mol new ${output_dir}/system.psf
mol addfile ${output_dir}/system.pdb

# Determine the size of the molecule
set sel [atomselect top all]
set minmax [measure minmax \$sel]
set vecDiff [vecsub [lindex \$minmax 1] [lindex \$minmax 0]]
set padding 10.0

# Create a water box that extends 10 Angstroms beyond the molecule
solvate ${output_dir}/system.psf ${output_dir}/system.pdb -t 10 -o ${output_dir}/solvated

# Add ions to neutralize the system
autoionize -psf ${output_dir}/solvated.psf -pdb ${output_dir}/solvated.pdb -o ${output_dir}/ionized -neutralize

exit
EOF
    
    local vmd_cmd="vmd -dispdev text -e $vmd_script"
    
    if [[ "$DRY_RUN" == "true" ]]; then
      log_info "[DRY-RUN] Would solvate and ionize the system using VMD"
    else
      log_command "$vmd_cmd"
      eval "$vmd_cmd" > "${output_dir}/vmd_solvate.log" 2>&1 || {
        log_error "Failed to solvate the system using VMD"
        return 1
      }
    fi
  else
    log_warning "VMD not found. Cannot solvate the system automatically."
    log_info "You may need to manually solvate the system using VMD or other tools."
  fi
  
  log_success "NAMD system setup complete"
  log_info "Setup files in: $output_dir"
  return 0
}