#!/bin/bash
#
# Force field parameterization module for MD simulation workflow
#
# Author: AhmedFikry90

# Assign force field parameters to the molecule
parameterize_structure() {
  log_info "Starting force field parameterization..."
  
  # Create output directory
  local output_dir=$(create_output_directory "parameterize")
  if [[ $? -ne 0 ]]; then
    log_error "Failed to create output directory"
    return 1
  fi
  
  # Validate input structure
  if [[ -z "$STRUCTURE_FILE" ]]; then
    log_error "No structure file specified. Please set STRUCTURE_FILE in your configuration."
    return 1
  fi
  
  if ! validate_file "$STRUCTURE_FILE" "structure file"; then
    log_error "Invalid or missing structure file: $STRUCTURE_FILE"
    return 1
  fi
  
  # Perform MD engine-specific parameterization
  case "${MD_ENGINE,,}" in
    gromacs)
      parameterize_gromacs "$STRUCTURE_FILE" "$output_dir"
      ;;
    amber)
      parameterize_amber "$STRUCTURE_FILE" "$output_dir"
      ;;
    namd)
      parameterize_namd "$STRUCTURE_FILE" "$output_dir"
      ;;
    *)
      log_error "Unsupported MD engine for parameterization: $MD_ENGINE"
      return 1
      ;;
  esac
  
  # Return the result of the parameterization
  local result=$?
  if [[ $result -eq 0 ]]; then
    log_success "Force field parameterization completed successfully"
  else
    log_error "Force field parameterization failed"
  fi
  
  return $result
}

# Parameterize for GROMACS simulations
parameterize_gromacs() {
  local input_file="$1"
  local output_dir="$2"
  
  log_info "Parameterizing for GROMACS..."
  
  # Path to preprocessed structure from prepare module
  local prepared_dir=$(find "$OUTPUT_DIR" -type d -name "prepare_*" | sort -r | head -n 1)
  
  if [[ -z "$prepared_dir" || ! -d "$prepared_dir" ]]; then
    log_warning "No prepared structure directory found. Running preparation first..."
    prepare_structure
    prepared_dir=$(find "$OUTPUT_DIR" -type d -name "prepare_*" | sort -r | head -n 1)
    
    if [[ -z "$prepared_dir" || ! -d "$prepared_dir" ]]; then
      log_error "Failed to prepare structure. Cannot continue parameterization."
      return 1
    fi
  fi
  
  # Check if box.gro exists from prepare module
  local struct_file="${prepared_dir}/box.gro"
  if [[ ! -f "$struct_file" ]]; then
    log_error "Expected structure file not found: $struct_file"
    return 1
  }
  
  # Check if topology file exists from prepare module
  local topol_file="${prepared_dir}/topol.top"
  if [[ ! -f "$topol_file" ]]; then
    log_error "Expected topology file not found: $topol_file"
    return 1
  }
  
  # Copy files to output directory
  if [[ "$DRY_RUN" == "true" ]]; then
    log_info "[DRY-RUN] Would copy prepared files to parameterize directory"
  else
    cp "$struct_file" "${output_dir}/system.gro"
    cp "$topol_file" "${output_dir}/topol.top"
    log_info "Copied prepared structure and topology to: $output_dir"
  fi
  
  # Generate position restraint file
  if command -v gmx &> /dev/null; then
    local genrestr_cmd="gmx genrestr -f ${output_dir}/system.gro -o ${output_dir}/posre.itp -fc 1000 1000 1000"
    
    if [[ "$DRY_RUN" == "true" ]]; then
      log_info "[DRY-RUN] Would generate position restraints"
    else
      log_command "$genrestr_cmd"
      echo "Protein" | eval "$genrestr_cmd" > "${output_dir}/genrestr.log" 2>&1 || {
        log_warning "Failed to generate position restraints. Continuing without them."
      }
    fi
  else
    log_error "GROMACS (gmx) command not found in PATH"
    return 1
  fi
  
  log_success "GROMACS parameterization complete"
  log_info "Parameterized files in: $output_dir"
  return 0
}

# Parameterize for AMBER simulations
parameterize_amber() {
  local input_file="$1"
  local output_dir="$2"
  
  log_info "Parameterizing for AMBER..."
  
  # Path to preprocessed structure from prepare module
  local prepared_dir=$(find "$OUTPUT_DIR" -type d -name "prepare_*" | sort -r | head -n 1)
  
  if [[ -z "$prepared_dir" || ! -d "$prepared_dir" ]]; then
    log_warning "No prepared structure directory found. Running preparation first..."
    prepare_structure
    prepared_dir=$(find "$OUTPUT_DIR" -type d -name "prepare_*" | sort -r | head -n 1)
    
    if [[ -z "$prepared_dir" || ! -d "$prepared_dir" ]]; then
      log_error "Failed to prepare structure. Cannot continue parameterization."
      return 1
    fi
  fi
  
  # Check if prmtop and inpcrd files exist from prepare module
  local prmtop_file="${prepared_dir}/system.prmtop"
  local inpcrd_file="${prepared_dir}/system.inpcrd"
  
  if [[ ! -f "$prmtop_file" || ! -f "$inpcrd_file" ]]; then
    log_error "Expected AMBER topology or coordinate files not found"
    return 1
  }
  
  # Copy files to output directory
  if [[ "$DRY_RUN" == "true" ]]; then
    log_info "[DRY-RUN] Would copy prepared files to parameterize directory"
  else
    cp "$prmtop_file" "${output_dir}/system.prmtop"
    cp "$inpcrd_file" "${output_dir}/system.inpcrd"
    log_info "Copied prepared topology and coordinates to: $output_dir"
  fi
  
  log_success "AMBER parameterization complete"
  log_info "Parameterized files in: $output_dir"
  return 0
}

# Parameterize for NAMD simulations
parameterize_namd() {
  local input_file="$1"
  local output_dir="$2"
  
  log_info "Parameterizing for NAMD..."
  
  # Path to preprocessed structure from prepare module
  local prepared_dir=$(find "$OUTPUT_DIR" -type d -name "prepare_*" | sort -r | head -n 1)
  
  if [[ -z "$prepared_dir" || ! -d "$prepared_dir" ]]; then
    log_warning "No prepared structure directory found. Running preparation first..."
    prepare_structure
    prepared_dir=$(find "$OUTPUT_DIR" -type d -name "prepare_*" | sort -r | head -n 1)
    
    if [[ -z "$prepared_dir" || ! -d "$prepared_dir" ]]; then
      log_error "Failed to prepare structure. Cannot continue parameterization."
      return 1
    fi
  fi
  
  # Check if PSF and PDB files exist from prepare module
  local psf_file="${prepared_dir}/system.psf"
  local pdb_file="${prepared_dir}/system.pdb"
  
  if [[ ! -f "$psf_file" || ! -f "$pdb_file" ]]; then
    log_error "Expected NAMD PSF or PDB files not found"
    return 1
  }
  
  # Copy files to output directory
  if [[ "$DRY_RUN" == "true" ]]; then
    log_info "[DRY-RUN] Would copy prepared files to parameterize directory"
  else
    cp "$psf_file" "${output_dir}/system.psf"
    cp "$pdb_file" "${output_dir}/system.pdb"
    log_info "Copied prepared PSF and PDB files to: $output_dir"
  fi
  
  # Create parameter files directory if it doesn't exist
  local param_dir="${output_dir}/parameters"
  mkdir -p "$param_dir"
  
  # Create a simple script to download CHARMM parameters if they don't exist
  local download_script="${output_dir}/download_parameters.sh"
  cat > "$download_script" << 'EOF'
#!/bin/bash

# Directory for parameter files
PARAM_DIR="parameters"
mkdir -p "$PARAM_DIR"

# List of CHARMM parameter files to download
PARAM_FILES=(
  "par_all36_prot.prm"
  "toppar_water_ions.str"
  "top_all36_prot.rtf"
)

# Download URL
BASE_URL="https://www.ks.uiuc.edu/Research/namd/2.12/release/lib/charmmpar/"

# Download each file
for file in "${PARAM_FILES[@]}"; do
  if [ ! -f "${PARAM_DIR}/${file}" ]; then
    echo "Downloading ${file}..."
    curl -s "${BASE_URL}/${file}" -o "${PARAM_DIR}/${file}"
    if [ $? -ne 0 ]; then
      echo "Failed to download ${file}"
    else
      echo "Downloaded ${file} successfully"
    fi
  else
    echo "${file} already exists, skipping download."
  fi
done

echo "Parameter files ready for NAMD simulation"
EOF
  
  # Make the script executable
  chmod +x "$download_script"
  
  # Execute the script if not in dry-run mode
  if [[ "$DRY_RUN" != "true" ]]; then
    log_info "Downloading parameter files for NAMD..."
    "$download_script" || {
      log_warning "Failed to download some parameter files. Manual configuration may be required."
    }
  else
    log_info "[DRY-RUN] Would download parameter files for NAMD"
  fi
  
  log_success "NAMD parameterization complete"
  log_info "Parameterized files in: $output_dir"
  return 0
}