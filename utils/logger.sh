#!/bin/bash
#
# Logging utilities for MD simulation workflow
#
# Author: AhmedFikry90

# Define color codes for terminal output
if [[ -t 1 ]]; then  # Check if stdout is a terminal
  COLOR_RESET="\033[0m"
  COLOR_RED="\033[0;31m"
  COLOR_GREEN="\033[0;32m"
  COLOR_YELLOW="\033[0;33m"
  COLOR_BLUE="\033[0;34m"
  COLOR_MAGENTA="\033[0;35m"
  COLOR_CYAN="\033[0;36m"
  COLOR_GRAY="\033[0;37m"
  COLOR_BOLD="\033[1m"
else
  # No colors if not a terminal
  COLOR_RESET=""
  COLOR_RED=""
  COLOR_GREEN=""
  COLOR_YELLOW=""
  COLOR_BLUE=""
  COLOR_MAGENTA=""
  COLOR_CYAN=""
  COLOR_GRAY=""
  COLOR_BOLD=""
fi

# Logging levels
LOG_LEVEL_NONE=0
LOG_LEVEL_ERROR=1
LOG_LEVEL_WARNING=2
LOG_LEVEL_INFO=3
LOG_LEVEL_DEBUG=4

# Initialize the log file
init_log_file() {
  if [[ -n "$LOG_FILE" ]]; then
    # Create log directory if it doesn't exist
    local log_dir=$(dirname "$LOG_FILE")
    if [[ ! -d "$log_dir" ]]; then
      mkdir -p "$log_dir"
    fi
    
    # Create or truncate log file
    echo "# MD Workflow Log - Started at $(date)" > "$LOG_FILE"
    echo "# =====================================" >> "$LOG_FILE"
  fi
}

# Log a message to stdout and/or log file
# Usage: _log_message LEVEL COLOR MESSAGE
_log_message() {
  local level="$1"
  local color="$2"
  local prefix="$3"
  local message="$4"
  
  # Skip if log level is too low
  if [[ "$level" -gt "$LOG_LEVEL" ]]; then
    return 0
  fi
  
  # Format message with timestamp if enabled
  local formatted_message=""
  if [[ "$LOG_TIMESTAMP" == "true" ]]; then
    formatted_message="[$(date '+%Y-%m-%d %H:%M:%S')] $prefix: $message"
  else
    formatted_message="$prefix: $message"
  fi
  
  # Print to stdout with color
  echo -e "${color}${formatted_message}${COLOR_RESET}"
  
  # Write to log file without color codes if specified
  if [[ -n "$LOG_FILE" ]]; then
    echo "$formatted_message" >> "$LOG_FILE"
  fi
}

# Log an error message
log_error() {
  _log_message $LOG_LEVEL_ERROR "$COLOR_RED" "ERROR" "$1"
}

# Log a warning message
log_warning() {
  _log_message $LOG_LEVEL_WARNING "$COLOR_YELLOW" "WARNING" "$1"
}

# Log an info message
log_info() {
  _log_message $LOG_LEVEL_INFO "$COLOR_BLUE" "INFO" "$1"
}

# Log a debug message
log_debug() {
  _log_message $LOG_LEVEL_DEBUG "$COLOR_GRAY" "DEBUG" "$1"
}

# Log a success message
log_success() {
  _log_message $LOG_LEVEL_INFO "$COLOR_GREEN" "SUCCESS" "$1"
}

# Log a command to be executed
log_command() {
  _log_message $LOG_LEVEL_DEBUG "$COLOR_CYAN" "EXEC" "$1"
}

# Execute a command with logging
# Usage: exec_cmd "command to execute"
exec_cmd() {
  local cmd="$1"
  
  log_command "$cmd"
  
  if [[ "$DRY_RUN" == "true" ]]; then
    log_info "[DRY-RUN] Would execute: $cmd"
    return 0
  fi
  
  if [[ "$VERBOSE" == "true" ]]; then
    eval "$cmd"
  else
    eval "$cmd" > /dev/null 2>&1
  fi
  
  local exit_code=$?
  if [[ $exit_code -ne 0 ]]; then
    log_error "Command failed with exit code $exit_code: $cmd"
    return $exit_code
  fi
  
  return 0
}

# Display a progress bar
# Usage: show_progress CURRENT TOTAL [DESCRIPTION]
show_progress() {
  local current=$1
  local total=$2
  local description="${3:-Processing}"
  local percent=$((current * 100 / total))
  local bar_length=50
  local filled_length=$((bar_length * current / total))
  
  # Skip if not in verbose mode or not a terminal
  if [[ "$VERBOSE" != "true" || ! -t 1 ]]; then
    return 0
  fi
  
  # Create the progress bar
  local bar=""
  for ((i=0; i<bar_length; i++)); do
    if [[ $i -lt $filled_length ]]; then
      bar+="▓"
    else
      bar+="░"
    fi
  done
  
  # Print the progress bar
  printf "\r%-20s [%s] %3d%%" "$description" "$bar" "$percent"
  
  # Print newline if complete
  if [[ $current -eq $total ]]; then
    echo
  fi
}

# Initialize log file on script load
init_log_file