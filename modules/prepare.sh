#!/bin/bash
#
# Structure preparation module for MD simulation workflow
#
# Author: AhmedFikry90

# Prepare molecule structure for simulation
prepare_structure() {
  log_info "Starting structure preparation..."
  
  # Create output directory
  local output_dir=$(create_output_directory "prepare")
  if [[ $? -ne 0 ]]; then
    log_error "Failed to create output directory"
    return 1
  fi
  
  # Validate input structure file
  if [[ -z "$STRUCTURE_FILE" ]]; then
    log_error "No structure file specified. Please set STRUCTURE_FILE in your configuration."
    return 1
  fi
  
  if ! validate_file "$STRUCTURE_FILE" "structure file"; then
    log_error "Invalid or missing structure file: $STRUCTURE_FILE"
    return 1
  fi
  
  log_info "Using structure file: $STRUCTURE_FILE"
  
  # Perform MD engine-specific preparation
  case "${MD_ENGINE,,}" in
    gromacs)
      prepare_structure_gromacs "$STRUCTURE_FILE" "$output_dir"
      ;;
    amber)
      prepare_structure_amber "$STRUCTURE_FILE" "$output_dir"
      ;;
    namd)
      prepare_structure_namd "$STRUCTURE_FILE" "$output_dir"
      ;;
    *)
      log_error "Unsupported MD engine for structure preparation: $MD_ENGINE"
      return 1
      ;;
  esac
  
  # Return the result of the preparation
  local result=$?
  if [[ $result -eq 0 ]]; then
    log_success "Structure preparation completed successfully"
  else
    log_error "Structure preparation failed"
  fi
  
  return $result
}

# Prepare structure for GROMACS simulations
prepare_structure_gromacs() {
  local input_file="$1"
  local output_dir="$2"
  
  log_info "Preparing structure for GROMACS..."
  
  # Determine file format
  local file_ext="${input_file##*.}"
  local output_file="${output_dir}/processed.gro"
  local clean_pdb="${TEMP_DIR}/clean.pdb"
  
  if [[ "${file_ext,,}" != "pdb" ]]; then
    log_info "Converting input file to PDB format..."
    if command -v obabel &> /dev/null; then
      if [[ "$DRY_RUN" == "true" ]]; then
        log_info "[DRY-RUN] Would convert $input_file to PDB format using Open Babel"
      else
        obabel "$input_file" -opdb -O "$clean_pdb" -h 2>/dev/null || {
          log_error "Failed to convert structure to PDB format"
          return 1
        }
      fi
    else
      log_warning "Open Babel not found. Using original file."
      clean_pdb="$input_file"
    fi
  else
    # Copy and clean the PDB file
    if [[ "$DRY_RUN" == "true" ]]; then
      log_info "[DRY-RUN] Would clean PDB file $input_file"
    else
      # Remove non-ATOM records (simple cleaning)
      grep "^ATOM" "$input_file" > "$clean_pdb" || {
        log_warning "No ATOM records found, using original file"
        cp "$input_file" "$clean_pdb"
      }
    fi
  fi
  
  # Create a topology using pdb2gmx
  if command -v gmx &> /dev/null; then
    local pdb2gmx_cmd="gmx pdb2gmx -f $clean_pdb -o $output_file -p ${output_dir}/topol.top -ff oplsaa -water spce -ignh"
    
    if [[ "$DRY_RUN" == "true" ]]; then
      log_info "[DRY-RUN] Would run: $pdb2gmx_cmd"
    else
      log_command "$pdb2gmx_cmd"
      # Auto-select the first option for all prompts
      echo "1" | eval "$pdb2gmx_cmd" > "${output_dir}/pdb2gmx.log" 2>&1 || {
        log_error "Failed to create GROMACS topology with pdb2gmx"
        return 1
      }
    fi
  else
    log_error "GROMACS (gmx) command not found in PATH"
    return 1
  fi
  
  # Create a box
  if [[ "$DRY_RUN" != "true" && -f "$output_file" ]]; then
    local editconf_cmd="gmx editconf -f $output_file -o ${output_dir}/box.gro -c -d 1.0 -bt cubic"
    
    log_command "$editconf_cmd"
    eval "$editconf_cmd" > "${output_dir}/editconf.log" 2>&1 || {
      log_error "Failed to create simulation box with editconf"
      return 1
    }
  elif [[ "$DRY_RUN" == "true" ]]; then
    log_info "[DRY-RUN] Would create simulation box"
  fi
  
  log_success "GROMACS structure preparation complete"
  log_info "Prepared files in: $output_dir"
  return 0
}

# Prepare structure for AMBER simulations
prepare_structure_amber() {
  local input_file="$1"
  local output_dir="$2"
  
  log_info "Preparing structure for AMBER..."
  
  # Determine file format
  local file_ext="${input_file##*.}"
  local clean_pdb="${TEMP_DIR}/clean.pdb"
  
  if [[ "${file_ext,,}" != "pdb" ]]; then
    log_info "Converting input file to PDB format..."
    if command -v obabel &> /dev/null; then
      if [[ "$DRY_RUN" == "true" ]]; then
        log_info "[DRY-RUN] Would convert $input_file to PDB format using Open Babel"
      else
        obabel "$input_file" -opdb -O "$clean_pdb" -h 2>/dev/null || {
          log_error "Failed to convert structure to PDB format"
          return 1
        }
      fi
    else
      log_warning "Open Babel not found. Using original file."
      cp "$input_file" "$clean_pdb"
    fi
  else
    # Copy the PDB file
    cp "$input_file" "$clean_pdb"
  fi
  
  # Run tleap to prepare the system (if AMBER tools are available)
  if command -v tleap &> /dev/null; then
    local tleap_input="${TEMP_DIR}/tleap.in"
    local tleap_cmd="tleap -f $tleap_input"
    
    # Create tleap input file
    cat > "$tleap_input" << EOF
source leaprc.protein.ff14SB
source leaprc.water.tip3p
mol = loadpdb $clean_pdb
check mol
saveamberparm mol ${output_dir}/system.prmtop ${output_dir}/system.inpcrd
savepdb mol ${output_dir}/system.pdb
quit
EOF
    
    if [[ "$DRY_RUN" == "true" ]]; then
      log_info "[DRY-RUN] Would run tleap to prepare AMBER system"
    else
      log_command "$tleap_cmd"
      eval "$tleap_cmd" > "${output_dir}/tleap.log" 2>&1 || {
        log_error "Failed to prepare AMBER system with tleap"
        return 1
      }
    fi
  else
    log_error "AMBER tleap command not found in PATH"
    return 1
  fi
  
  log_success "AMBER structure preparation complete"
  log_info "Prepared files in: $output_dir"
  return 0
}

# Prepare structure for NAMD simulations
prepare_structure_namd() {
  local input_file="$1"
  local output_dir="$2"
  
  log_info "Preparing structure for NAMD..."
  
  # Determine file format
  local file_ext="${input_file##*.}"
  local clean_pdb="${TEMP_DIR}/clean.pdb"
  
  if [[ "${file_ext,,}" != "pdb" ]]; then
    log_info "Converting input file to PDB format..."
    if command -v obabel &> /dev/null; then
      if [[ "$DRY_RUN" == "true" ]]; then
        log_info "[DRY-RUN] Would convert $input_file to PDB format using Open Babel"
      else
        obabel "$input_file" -opdb -O "$clean_pdb" -h 2>/dev/null || {
          log_error "Failed to convert structure to PDB format"
          return 1
        }
      fi
    else
      log_warning "Open Babel not found. Using original file."
      cp "$input_file" "$clean_pdb"
    fi
  else
    # Copy the PDB file
    cp "$input_file" "$clean_pdb"
  fi
  
  # Run psfgen to prepare the system (if VMD/NAMD tools are available)
  if command -v psfgen &> /dev/null; then
    local psfgen_input="${TEMP_DIR}/psfgen.tcl"
    local psfgen_cmd="psfgen < $psfgen_input"
    
    # Create psfgen input file
    cat > "$psfgen_input" << EOF
package require psfgen
topology toppar/top_all36_prot.rtf
pdbalias residue HIS HSE
pdbalias atom ILE CD1 CD
segment PROT {pdb $clean_pdb}
coordpdb $clean_pdb PROT
guesscoord
writepdb ${output_dir}/system.pdb
writepsf ${output_dir}/system.psf
exit
EOF
    
    if [[ "$DRY_RUN" == "true" ]]; then
      log_info "[DRY-RUN] Would run psfgen to prepare NAMD system"
    else
      log_command "$psfgen_cmd"
      eval "$psfgen_cmd" > "${output_dir}/psfgen.log" 2>&1 || {
        log_error "Failed to prepare NAMD system with psfgen"
        return 1
      }
    fi
  else
    log_error "NAMD psfgen command not found in PATH"
    return 1
  fi
  
  log_success "NAMD structure preparation complete"
  log_info "Prepared files in: $output_dir"
  return 0
}