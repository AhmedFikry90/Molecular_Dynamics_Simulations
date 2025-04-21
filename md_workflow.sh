#!/bin/bash
#
# Main controller script for Molecular Dynamics simulations workflow
# 
# Author: AhmedFikry90

# Set strict mode
set -euo pipefail

# Load configuration and utilities
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"
source "${SCRIPT_DIR}/utils/logger.sh"
source "${SCRIPT_DIR}/utils/validation.sh"
source "${SCRIPT_DIR}/utils/io_handler.sh"

# Display banner
display_banner() {
  log_info "================================================================="
  log_info "                Molecular Dynamics Workflow v1.0                 "
  log_info "================================================================="
  log_info "Started at $(date)"
  log_info "Working directory: $(pwd)"
  log_info "================================================================="
}

# Display usage information
display_usage() {
  cat << EOF
Usage: $(basename "$0") [OPTIONS] COMMAND

A comprehensive workflow for Molecular Dynamics simulations.

OPTIONS:
  -h, --help                Show this help message
  -c, --config FILE         Use custom configuration file
  -v, --verbose             Enable verbose output
  -f, --force               Force overwrite of existing files
  -d, --dry-run             Show commands without executing

COMMANDS:
  init PROJECT_NAME         Initialize a new project
  prepare                   Prepare molecule structure
  parameterize              Assign force field parameters
  setup                     Set up simulation system
  minimize                  Run energy minimization
  equilibrate               Run system equilibration
  production                Run production simulation
  analyze                   Analyze simulation results
  clean                     Clean temporary files
  all                       Run complete workflow

Examples:
  $(basename "$0") init my_protein
  $(basename "$0") --config custom_config.sh all
EOF
}

# Check for required dependencies
check_dependencies() {
  log_info "Checking dependencies..."
  
  # Basic utilities
  for cmd in awk sed grep cut; do
    if ! command -v "$cmd" &> /dev/null; then
      log_error "Required command not found: $cmd"
      exit 1
    fi
  done
  
  # Check for MD software based on configuration
  case "$MD_ENGINE" in
    gromacs)
      if ! command -v gmx &> /dev/null; then
        log_warning "GROMACS (gmx) not found in PATH"
        log_info "Please ensure GROMACS is installed or loaded via a module system"
      else
        log_success "GROMACS found: $(command -v gmx)"
      fi
      ;;
    amber)
      if ! command -v pmemd &> /dev/null; then
        log_warning "AMBER (pmemd) not found in PATH"
        log_info "Please ensure AMBER is installed or loaded via a module system"
      else
        log_success "AMBER found: $(command -v pmemd)"
      fi
      ;;
    namd)
      if ! command -v namd2 &> /dev/null; then
        log_warning "NAMD (namd2) not found in PATH"
        log_info "Please ensure NAMD is installed or loaded via a module system"
      else
        log_success "NAMD found: $(command -v namd2)"
      fi
      ;;
    *)
      log_warning "Unknown MD engine: $MD_ENGINE"
      log_info "Please specify a supported MD engine in config.sh"
      ;;
  esac
}

# Initialize a new project
init_project() {
  local project_name="$1"
  
  if [[ -z "$project_name" ]]; then
    log_error "Project name is required"
    display_usage
    exit 1
  fi
  
  if [[ -d "$project_name" ]]; then
    if [[ "$FORCE_OVERWRITE" == "true" ]]; then
      log_warning "Project directory already exists, overwriting as requested"
    else
      log_error "Project directory already exists: $project_name"
      log_info "Use -f/--force to overwrite existing project"
      exit 1
    fi
  fi
  
  log_info "Initializing project: $project_name"
  
  mkdir -p "$project_name"/{input,output,logs,temp}
  
  # Create project configuration file
  cat > "$project_name/project_config.sh" << EOF
#!/bin/bash
# Project configuration for $project_name
# Created: $(date)

# Project identification
PROJECT_NAME="$project_name"
PROJECT_DESCRIPTION=""

# Molecule information
STRUCTURE_FILE=""
TOPOLOGY_FILE=""

# Simulation parameters
SIMULATION_TIME=10        # ns
TIMESTEP=0.002            # ps
TEMPERATURE=300           # K
PRESSURE=1.0              # bar

# Output frequency
ENERGY_OUTPUT_FREQ=1000   # steps
TRAJECTORY_OUTPUT_FREQ=1000  # steps

# Resources
NUM_THREADS=$NUM_THREADS
NUM_GPUS=$NUM_GPUS
MPI_PROCESSES=$MPI_PROCESSES
EOF
  
  # Create README file
  cat > "$project_name/README.md" << EOF
# $project_name

Molecular Dynamics simulation project

## Project Structure

- \`input/\`: Input files and structures
- \`output/\`: Simulation results
- \`logs/\`: Log files
- \`temp/\`: Temporary files

## Configuration

Edit \`project_config.sh\` to configure simulation parameters.

## Running Simulations

To run a complete simulation workflow:

\`\`\`bash
../md_workflow.sh -c project_config.sh all
\`\`\`

## Notes


EOF
  
  log_success "Project initialized successfully: $project_name"
  log_info "Edit $project_name/project_config.sh to configure your simulation"
}

# Parse command line options
parse_args() {
  POSITIONAL_ARGS=()
  
  while [[ $# -gt 0 ]]; do
    case $1 in
      -h|--help)
        display_usage
        exit 0
        ;;
      -c|--config)
        if [[ -z "$2" || "$2" == -* ]]; then
          log_error "Option $1 requires an argument"
          exit 1
        fi
        CONFIG_FILE="$2"
        if [[ ! -f "$CONFIG_FILE" ]]; then
          log_error "Config file not found: $CONFIG_FILE"
          exit 1
        fi
        source "$CONFIG_FILE"
        shift 2
        ;;
      -v|--verbose)
        VERBOSE="true"
        shift
        ;;
      -f|--force)
        FORCE_OVERWRITE="true"
        shift
        ;;
      -d|--dry-run)
        DRY_RUN="true"
        shift
        ;;
      -*|--*)
        log_error "Unknown option: $1"
        display_usage
        exit 1
        ;;
      *)
        POSITIONAL_ARGS+=("$1")
        shift
        ;;
    esac
  done
  
  set -- "${POSITIONAL_ARGS[@]}"
  
  if [[ $# -lt 1 ]]; then
    log_error "No command specified"
    display_usage
    exit 1
  fi
  
  COMMAND="$1"
  shift
  COMMAND_ARGS=("$@")
}

# Execute commands
run_command() {
  case "$COMMAND" in
    init)
      if [[ ${#COMMAND_ARGS[@]} -lt 1 ]]; then
        log_error "Project name is required for init command"
        display_usage
        exit 1
      fi
      init_project "${COMMAND_ARGS[0]}"
      ;;
    prepare)
      source "${SCRIPT_DIR}/modules/prepare.sh"
      prepare_structure
      ;;
    parameterize)
      source "${SCRIPT_DIR}/modules/parameterize.sh"
      parameterize_structure
      ;;
    setup)
      source "${SCRIPT_DIR}/modules/setup.sh"
      setup_system
      ;;
    minimize)
      source "${SCRIPT_DIR}/modules/minimize.sh"
      minimize_energy
      ;;
    equilibrate)
      source "${SCRIPT_DIR}/modules/equilibrate.sh"
      equilibrate_system
      ;;
    production)
      source "${SCRIPT_DIR}/modules/production.sh"
      run_production
      ;;
    analyze)
      source "${SCRIPT_DIR}/modules/analyze.sh"
      analyze_results
      ;;
    clean)
      clean_temp_files
      ;;
    all)
      log_info "Running complete workflow..."
      source "${SCRIPT_DIR}/modules/prepare.sh"
      source "${SCRIPT_DIR}/modules/parameterize.sh"
      source "${SCRIPT_DIR}/modules/setup.sh"
      source "${SCRIPT_DIR}/modules/minimize.sh"
      source "${SCRIPT_DIR}/modules/equilibrate.sh"
      source "${SCRIPT_DIR}/modules/production.sh"
      source "${SCRIPT_DIR}/modules/analyze.sh"
      
      prepare_structure
      parameterize_structure
      setup_system
      minimize_energy
      equilibrate_system
      run_production
      analyze_results
      log_success "Complete workflow finished"
      ;;
    *)
      log_error "Unknown command: $COMMAND"
      display_usage
      exit 1
      ;;
  esac
}

# Clean temporary files
clean_temp_files() {
  log_info "Cleaning temporary files..."
  
  if [[ "$DRY_RUN" == "true" ]]; then
    log_info "[DRY-RUN] Would remove files in temp/"
  else
    find "temp/" -type f -not -name "*.keep" -delete
    log_success "Temporary files cleaned"
  fi
}

# Main function
main() {
  display_banner
  check_dependencies
  
  # Execute the requested command
  run_command
  
  log_info "Workflow completed successfully"
}

# Parse command line arguments and run main function
parse_args "$@"
main