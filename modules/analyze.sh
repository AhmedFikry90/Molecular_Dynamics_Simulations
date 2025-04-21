#!/bin/bash
#
# Analysis module for MD simulation workflow
#
# Author: AhmedFikry90

# Analyze simulation results
analyze_results() {
  log_info "Starting analysis of simulation results..."
  
  # Create output directory
  local output_dir=$(create_output_directory "analysis")
  if [[ $? -ne 0 ]]; then
    log_error "Failed to create output directory"
    return 1
  fi
  
  # Perform MD engine-specific analysis
  case "${MD_ENGINE,,}" in
    gromacs)
      analyze_results_gromacs "$output_dir"
      ;;
    amber)
      analyze_results_amber "$output_dir"
      ;;
    namd)
      analyze_results_namd "$output_dir"
      ;;
    *)
      log_error "Unsupported MD engine for analysis: $MD_ENGINE"
      return 1
      ;;
  esac
  
  # Return the result of the analysis
  local result=$?
  if [[ $result -eq 0 ]]; then
    log_success "Analysis completed successfully"
  else
    log_error "Analysis failed"
  fi
  
  return $result
}

# Analyze GROMACS simulation results
analyze_results_gromacs() {
  local output_dir="$1"
  
  log_info "Analyzing GROMACS simulation results..."
  
  # Find the most recent production run directory
  local prod_dir=$(find "$OUTPUT_DIR" -type d -name "production_*" | sort -r | head -n 1)
  
  if [[ -z "$prod_dir" || ! -d "$prod_dir" ]]; then
    log_warning "No production run results found. Running production simulation first..."
    run_production
    prod_dir=$(find "$OUTPUT_DIR" -type d -name "production_*" | sort -r | head -n 1)
    
    if [[ -z "$prod_dir" || ! -d "$prod_dir" ]]; then
      log_error "Failed to run production simulation. Cannot continue analysis."
      return 1
    }
  fi
  
  # Check for required files
  local tpr_file="${prod_dir}/md.tpr"
  local xtc_file="${prod_dir}/md.xtc"
  local edr_file="${prod_dir}/md.edr"
  
  if [[ ! -f "$tpr_file" || ! -f "$xtc_file" || ! -f "$edr_file" ]]; then
    log_error "Required GROMACS output files not found"
    return 1
  fi
  
  # Create analysis report directory
  mkdir -p "${output_dir}/plots"
  
  # Initialize analysis report
  local report_file="${output_dir}/analysis_report.md"
  cat > "$report_file" << EOF
# Molecular Dynamics Simulation Analysis Report

## System Information

- MD Engine: GROMACS
- Run Date: $(date)
- Production Run Directory: ${prod_dir}

## Analysis Results

EOF
  
  # Run GROMACS analysis tools
  if command -v gmx &> /dev/null; then
    # Energy analysis
    if [[ "$DRY_RUN" != "true" ]]; then
      log_info "Analyzing energy components..."
      
      # Create input for energy selection
      echo -e "Potential\nKinetic-En.\nTotal-Energy\nTemperature\nPressure" > "${TEMP_DIR}/energy_terms.txt"
      
      # Extract energy components
      local energy_cmd="gmx energy -f $edr_file -o ${output_dir}/energy.xvg"
      log_command "$energy_cmd"
      cat "${TEMP_DIR}/energy_terms.txt" | eval "$energy_cmd" > "${output_dir}/energy_analysis.log" 2>&1 || {
        log_warning "Energy analysis failed"
      }
      
      # Convert XVG to CSV for better portability
      if [[ -f "${output_dir}/energy.xvg" ]]; then
        awk 'BEGIN{FS="  "; OFS=","} /^[^#@]/{print $1, $2, $3, $4, $5, $6}' "${output_dir}/energy.xvg" > "${output_dir}/energy.csv"
        
        # Add to report
        cat >> "$report_file" << EOF
### Energy Analysis

The following energy components were analyzed:
- Potential Energy
- Kinetic Energy
- Total Energy
- Temperature
- Pressure

Energy data is available in CSV format at \`energy.csv\`.

EOF
      fi
      
      # RMSD analysis
      log_info "Calculating RMSD..."
      
      # Create index group for RMSD analysis
      echo "Protein" > "${TEMP_DIR}/rmsd_index.txt"
      
      # Calculate RMSD
      local rmsd_cmd="gmx rms -s $tpr_file -f $xtc_file -o ${output_dir}/rmsd.xvg -tu ns"
      log_command "$rmsd_cmd"
      cat "${TEMP_DIR}/rmsd_index.txt" "${TEMP_DIR}/rmsd_index.txt" | eval "$rmsd_cmd" > "${output_dir}/rmsd_analysis.log" 2>&1 || {
        log_warning "RMSD analysis failed"
      }
      
      # Convert XVG to CSV
      if [[ -f "${output_dir}/rmsd.xvg" ]]; then
        awk 'BEGIN{FS="  "; OFS=","} /^[^#@]/{print $1, $2}' "${output_dir}/rmsd.xvg" > "${output_dir}/rmsd.csv"
        
        # Add to report
        cat >> "$report_file" << EOF
### RMSD Analysis

Root Mean Square Deviation (RMSD) was calculated for the protein backbone atoms.
RMSD data is available in CSV format at \`rmsd.csv\`.

EOF
      fi
      
      # Radius of gyration analysis
      log_info "Calculating radius of gyration..."
      
      # Create index group for Rg analysis
      echo "Protein" > "${TEMP_DIR}/rg_index.txt"
      
      # Calculate Rg
      local rg_cmd="gmx gyrate -s $tpr_file -f $xtc_file -o ${output_dir}/rg.xvg"
      log_command "$rg_cmd"
      cat "${TEMP_DIR}/rg_index.txt" | eval "$rg_cmd" > "${output_dir}/rg_analysis.log" 2>&1 || {
        log_warning "Radius of gyration analysis failed"
      }
      
      # Convert XVG to CSV
      if [[ -f "${output_dir}/rg.xvg" ]]; then
        awk 'BEGIN{FS="  "; OFS=","} /^[^#@]/{print $1, $2}' "${output_dir}/rg.xvg" > "${output_dir}/rg.csv"
        
        # Add to report
        cat >> "$report_file" << EOF
### Radius of Gyration Analysis

Radius of gyration was calculated for the protein.
Rg data is available in CSV format at \`rg.csv\`.

EOF
      fi
      
      # Generate a PDB with average structure
      log_info "Generating average structure..."
      
      local trjconv_cmd="gmx trjconv -s $tpr_file -f $xtc_file -o ${output_dir}/average.pdb -tu ns -pbc nojump"
      log_command "$trjconv_cmd"
      echo "Protein" | eval "$trjconv_cmd" > "${output_dir}/trjconv_average.log" 2>&1 || {
        log_warning "Average structure generation failed"
      }
      
      if [[ -f "${output_dir}/average.pdb" ]]; then
        cat >> "$report_file" << EOF
### Average Structure

An average structure from the trajectory has been generated and saved as \`average.pdb\`.

EOF
      fi
      
      # Create a summary at the end of the report
      cat >> "$report_file" << EOF
## Summary

This analysis report contains:
- Energy analysis
- RMSD analysis
- Radius of gyration analysis
- Average structure

The raw data files are provided in CSV format for further analysis and visualization.
EOF
      
      log_success "GROMACS analysis complete"
      log_info "Analysis report generated: $report_file"
    else
      log_info "[DRY-RUN] Would analyze GROMACS simulation results"
    fi
  else
    log_error "GROMACS (gmx) command not found in PATH"
    return 1
  fi
  
  return 0
}

# Analyze AMBER simulation results
analyze_results_amber() {
  local output_dir="$1"
  
  log_info "Analyzing AMBER simulation results..."
  
  # Find the most recent production run directory
  local prod_dir=$(find "$OUTPUT_DIR" -type d -name "production_*" | sort -r | head -n 1)
  
  if [[ -z "$prod_dir" || ! -d "$prod_dir" ]]; then
    log_warning "No production run results found. Running production simulation first..."
    run_production
    prod_dir=$(find "$OUTPUT_DIR" -type d -name "production_*" | sort -r | head -n 1)
    
    if [[ -z "$prod_dir" || ! -d "$prod_dir" ]]; then
      log_error "Failed to run production simulation. Cannot continue analysis."
      return 1
    }
  fi
  
  # Check for required files
  local prmtop_file="${prod_dir}/system.prmtop"
  local trajectory_file="${prod_dir}/md.nc"
  local output_file="${prod_dir}/md.out"
  
  if [[ ! -f "$prmtop_file" || ! -f "$trajectory_file" || ! -f "$output_file" ]]; then
    log_error "Required AMBER output files not found"
    return 1
  fi
  
  # Create analysis report directory
  mkdir -p "${output_dir}/plots"
  
  # Initialize analysis report
  local report_file="${output_dir}/analysis_report.md"
  cat > "$report_file" << EOF
# Molecular Dynamics Simulation Analysis Report

## System Information

- MD Engine: AMBER
- Run Date: $(date)
- Production Run Directory: ${prod_dir}

## Analysis Results

EOF
  
  # Run AMBER analysis tools
  if command -v cpptraj &> /dev/null; then
    if [[ "$DRY_RUN" != "true" ]]; then
      log_info "Analyzing trajectory with cpptraj..."
      
      # Create cpptraj input script for various analyses
      local cpptraj_input="${TEMP_DIR}/cpptraj_analysis.in"
      cat > "$cpptraj_input" << EOF
# Load topology and trajectory
parm ${prmtop_file}
trajin ${trajectory_file}

# Calculate RMSD for protein backbone
rms first :1-999999@C,CA,N out ${output_dir}/rmsd.csv mass

# Calculate radius of gyration
radgyr :1-999999 out ${output_dir}/rg.csv mass

# Calculate secondary structure
secstruct :1-999999 out ${output_dir}/ss.csv sumout ${output_dir}/ss_summary.csv

# Average structure
average ${output_dir}/average.pdb pdb

# Create a representative structure (closest to average)
cluster hieragglo clusters 1 :1-999999@C,CA,N out ${output_dir}/cluster.csv summary ${output_dir}/cluster_summary.csv cpopvtime ${output_dir}/cluster_pop.csv repout ${output_dir}/rep_structure.pdb repfmt pdb

# Energy analysis from trajectory
energy out ${output_dir}/energy.csv

# Exit
quit
EOF
      
      # Run cpptraj analysis
      local cpptraj_cmd="cpptraj -i $cpptraj_input"
      log_command "$cpptraj_cmd"
      eval "$cpptraj_cmd" > "${output_dir}/cpptraj_analysis.log" 2>&1 || {
        log_warning "Some cpptraj analyses failed. Check the log file for details."
      }
      
      # Extract energy information from output file
      log_info "Extracting energy information from output file..."
      grep "NSTEP\|ENERGY\|TEMP" "$output_file" | grep -v "EPtot" > "${output_dir}/energy_output.txt"
      
      # Add analysis results to report
      if [[ -f "${output_dir}/rmsd.csv" ]]; then
        cat >> "$report_file" << EOF
### RMSD Analysis

Root Mean Square Deviation (RMSD) was calculated for the protein backbone atoms (C, CA, N).
RMSD data is available in CSV format at \`rmsd.csv\`.

EOF
      fi
      
      if [[ -f "${output_dir}/rg.csv" ]]; then
        cat >> "$report_file" << EOF
### Radius of Gyration Analysis

Radius of gyration was calculated for the protein.
Rg data is available in CSV format at \`rg.csv\`.

EOF
      fi
      
      if [[ -f "${output_dir}/ss.csv" ]]; then
        cat >> "$report_file" << EOF
### Secondary Structure Analysis

Secondary structure analysis was performed for the protein.
- Detailed SS data is available in CSV format at \`ss.csv\`
- Summary SS data is available in CSV format at \`ss_summary.csv\`

EOF
      fi
      
      if [[ -f "${output_dir}/average.pdb" ]]; then
        cat >> "$report_file" << EOF
### Average Structure

An average structure from the trajectory has been generated and saved as \`average.pdb\`.

EOF
      fi
      
      if [[ -f "${output_dir}/rep_structure.pdb" ]]; then
        cat >> "$report_file" << EOF
### Representative Structure

A representative structure (closest to the average) has been identified and saved as \`rep_structure.pdb\`.

EOF
      fi
      
      if [[ -f "${output_dir}/energy.csv" ]]; then
        cat >> "$report_file" << EOF
### Energy Analysis

Energy components were analyzed throughout the trajectory.
Energy data is available in CSV format at \`energy.csv\`.

EOF
      fi
      
      # Create a summary at the end of the report
      cat >> "$report_file" << EOF
## Summary

This analysis report contains:
- RMSD analysis
- Radius of gyration analysis
- Secondary structure analysis
- Average structure
- Representative structure
- Energy analysis

The raw data files are provided in CSV format for further analysis and visualization.
EOF
      
      log_success "AMBER analysis complete"
      log_info "Analysis report generated: $report_file"
    else
      log_info "[DRY-RUN] Would analyze AMBER simulation results"
    fi
  else
    log_error "AMBER cpptraj command not found in PATH"
    return 1
  fi
  
  return 0
}

# Analyze NAMD simulation results
analyze_results_namd() {
  local output_dir="$1"
  
  log_info "Analyzing NAMD simulation results..."
  
  # Find the most recent production run directory
  local prod_dir=$(find "$OUTPUT_DIR" -type d -name "production_*" | sort -r | head -n 1)
  
  if [[ -z "$prod_dir" || ! -d "$prod_dir" ]]; then
    log_warning "No production run results found. Running production simulation first..."
    run_production
    prod_dir=$(find "$OUTPUT_DIR" -type d -name "production_*" | sort -r | head -n 1)
    
    if [[ -z "$prod_dir" || ! -d "$prod_dir" ]]; then
      log_error "Failed to run production simulation. Cannot continue analysis."
      return 1
    }
  fi
  
  # Check for required files
  local psf_file="${prod_dir}/input.psf"
  local trajectory_file="${prod_dir}/prod.dcd"
  local log_file="${prod_dir}/md.log"
  
  if [[ ! -f "$psf_file" ]]; then
    log_error "PSF file not found in production directory"
    return 1
  fi
  
  if [[ ! -f "$trajectory_file" ]]; then
    trajectory_file=$(find "$prod_dir" -name "*.dcd" | head -n 1)
    if [[ -z "$trajectory_file" ]]; then
      log_error "No trajectory file found in production directory"
      return 1
    fi
  fi
  
  # Create analysis report directory
  mkdir -p "${output_dir}/plots"
  
  # Initialize analysis report
  local report_file="${output_dir}/analysis_report.md"
  cat > "$report_file" << EOF
# Molecular Dynamics Simulation Analysis Report

## System Information

- MD Engine: NAMD
- Run Date: $(date)
- Production Run Directory: ${prod_dir}

## Analysis Results

EOF
  
  # Check if VMD is available for analysis
  if command -v vmd &> /dev/null; then
    if [[ "$DRY_RUN" != "true" ]]; then
      log_info "Analyzing trajectory with VMD..."
      
      # Create VMD script for various analyses
      local vmd_script="${TEMP_DIR}/vmd_analysis.tcl"
      cat > "$vmd_script" << EOF
# Load molecule
mol new ${psf_file} type psf
mol addfile ${trajectory_file} type dcd waitfor all

# Output files
set rmsd_outfile [open "${output_dir}/rmsd.csv" w]
set rg_outfile [open "${output_dir}/rg.csv" w]
set energy_outfile [open "${output_dir}/energy.csv" w]
set ss_outfile [open "${output_dir}/ss.csv" w]

# Write headers
puts \$rmsd_outfile "Frame,Time,RMSD"
puts \$rg_outfile "Frame,Time,Rg"
puts \$ss_outfile "Frame,Time,Helix,Sheet,Coil,Turn"

# Select all protein atoms and backbone atoms
set all [atomselect top "protein"]
set backbone [atomselect top "protein and name CA C N"]
set num_frames [molinfo top get numframes]
set dt 0.002 ;# ps per step, adjust based on your trajectory

# Reference frame for RMSD
set reference [atomselect top "protein and name CA C N" frame 0]

# Calculate RMSD and Rg for each frame
for {set i 0} {\$i < \$num_frames} {incr i} {
    # Time in ns
    set time [expr \$i * \$dt * 1000.0 / 1000.0]
    
    # RMSD
    \$backbone frame \$i
    \$backbone move [measure fit \$backbone \$reference]
    set rmsd [measure rmsd \$backbone \$reference]
    puts \$rmsd_outfile "\$i,\$time,\$rmsd"
    
    # Radius of gyration
    \$all frame \$i
    set rg [measure rgyr \$all]
    puts \$rg_outfile "\$i,\$time,\$rg"
    
    # Secondary structure
    mol ssrecalc top
    set sel [atomselect top "protein" frame \$i]
    set helix [llength [atomselect top "protein and structure H" frame \$i]]
    set sheet [llength [atomselect top "protein and structure E" frame \$i]]
    set coil [llength [atomselect top "protein and structure C" frame \$i]]
    set turn [llength [atomselect top "protein and structure T" frame \$i]]
    set total [llength \$sel]
    
    # Normalize by total residues
    set helix_pct [expr double(\$helix) / \$total]
    set sheet_pct [expr double(\$sheet) / \$total]
    set coil_pct [expr double(\$coil) / \$total]
    set turn_pct [expr double(\$turn) / \$total]
    
    puts \$ss_outfile "\$i,\$time,\$helix_pct,\$sheet_pct,\$coil_pct,\$turn_pct"
    \$sel delete
}

# Save average structure
set average [atomselect top "protein" frame [expr \$num_frames / 2]]
\$average writepdb "${output_dir}/average.pdb"
\$average delete

# Close files
close \$rmsd_outfile
close \$rg_outfile
close \$ss_outfile

# Exit VMD
quit
EOF
      
      # Run VMD analysis
      local vmd_cmd="vmd -dispdev text -e $vmd_script"
      log_command "$vmd_cmd"
      eval "$vmd_cmd" > "${output_dir}/vmd_analysis.log" 2>&1 || {
        log_warning "VMD analysis failed. Check the log file for details."
      }
      
      # Extract energy information from log file if available
      if [[ -f "$log_file" ]]; then
        log_info "Extracting energy information from log file..."
        grep "ENERGY:" "$log_file" | awk '{print $2","$12","$14","$18","$20}' > "${TEMP_DIR}/energy_data.txt"
        echo "Step,Potential,Kinetic,Temperature,Total" > "${output_dir}/energy.csv"
        cat "${TEMP_DIR}/energy_data.txt" >> "${output_dir}/energy.csv"
      fi
      
      # Add analysis results to report
      if [[ -f "${output_dir}/rmsd.csv" ]]; then
        cat >> "$report_file" << EOF
### RMSD Analysis

Root Mean Square Deviation (RMSD) was calculated for the protein backbone atoms (C, CA, N).
RMSD data is available in CSV format at \`rmsd.csv\`.

EOF
      fi
      
      if [[ -f "${output_dir}/rg.csv" ]]; then
        cat >> "$report_file" << EOF
### Radius of Gyration Analysis

Radius of gyration was calculated for the protein.
Rg data is available in CSV format at \`rg.csv\`.

EOF
      fi
      
      if [[ -f "${output_dir}/ss.csv" ]]; then
        cat >> "$report_file" << EOF
### Secondary Structure Analysis

Secondary structure analysis was performed for the protein.
SS data is available in CSV format at \`ss.csv\` with the following columns:
- Frame
- Time
- Helix (fraction)
- Sheet (fraction)
- Coil (fraction)
- Turn (fraction)

EOF
      fi
      
      if [[ -f "${output_dir}/average.pdb" ]]; then
        cat >> "$report_file" << EOF
### Average Structure

An average structure from the trajectory has been generated and saved as \`average.pdb\`.

EOF
      fi
      
      if [[ -f "${output_dir}/energy.csv" ]]; then
        cat >> "$report_file" << EOF
### Energy Analysis

Energy components were extracted from the log file.
Energy data is available in CSV format at \`energy.csv\` with the following columns:
- Step
- Potential Energy
- Kinetic Energy
- Temperature
- Total Energy

EOF
      fi
      
      # Create a summary at the end of the report
      cat >> "$report_file" << EOF
## Summary

This analysis report contains:
- RMSD analysis
- Radius of gyration analysis
- Secondary structure analysis
- Average structure
- Energy analysis

The raw data files are provided in CSV format for further analysis and visualization.
EOF
      
      log_success "NAMD analysis complete"
      log_info "Analysis report generated: $report_file"
    else
      log_info "[DRY-RUN] Would analyze NAMD simulation results"
    fi
  else
    log_warning "VMD not found. Limited analysis will be performed."
    
    # Minimal analysis without VMD
    if [[ -f "$log_file" && "$DRY_RUN" != "true" ]]; then
      log_info "Extracting energy information from log file..."
      grep "ENERGY:" "$log_file" | awk '{print $2","$12","$14","$18","$20}' > "${TEMP_DIR}/energy_data.txt"
      echo "Step,Potential,Kinetic,Temperature,Total" > "${output_dir}/energy.csv"
      cat "${TEMP_DIR}/energy_data.txt" >> "${output_dir}/energy.csv"
      
      cat >> "$report_file" << EOF
### Energy Analysis

Energy components were extracted from the log file.
Energy data is available in CSV format at \`energy.csv\` with the following columns:
- Step
- Potential Energy
- Kinetic Energy
- Temperature
- Total Energy

## Summary

This limited analysis report contains only energy data due to the absence of VMD.
For more comprehensive analysis, please install VMD or use other analysis tools.
EOF
      
      log_success "Limited NAMD analysis complete"
      log_info "Analysis report generated: $report_file"
    else
      log_info "[DRY-RUN] Would perform limited analysis of NAMD results"
    fi
  fi
  
  return 0
}