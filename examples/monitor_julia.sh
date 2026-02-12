#!/bin/bash

# Julia Process Monitor (Bash version)
# Monitors a Julia script execution, tracking memory, CPU, GPU usage, and precompilation status

set -e

# Configuration
SCRIPT_TO_RUN=""
INTERVAL=5
PROJECT_PATH=""
JULIA_ARGS=""
GPU_LOG_INTERVAL=0.2  # How often to log GPU usage (in seconds)
GPU_LOG_FILE=""
GPU_LOGGER_PID=""

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color

# Function to print usage
usage() {
    cat << EOF
Usage: $0 <script.jl> [options]

Monitor a Julia script execution with detailed diagnostics.

Arguments:
  script.jl           Julia script to run and monitor

Options:
  -i, --interval N    Monitoring interval in seconds (default: 5)
  -p, --project PATH  Julia project path
  -h, --help          Show this help message

Features:
  - Memory and CPU monitoring at regular intervals
  - Continuous GPU usage logging (0.2s interval) to CSV file
  - Precompilation detection
  - CUDA library loading detection
  - GPU usage summary at completion

Examples:
  $0 run_tjlf_timing.jl
  $0 run_tjlf_timing.jl -i 10
  $0 run_tjlf_timing.jl --project=.. -i 5

EOF
    exit 1
}

# Parse arguments
if [ $# -lt 1 ]; then
    usage
fi

SCRIPT_TO_RUN="$1"
shift

while [ $# -gt 0 ]; do
    case "$1" in
        -i|--interval)
            INTERVAL="$2"
            shift 2
            ;;
        -p|--project)
            PROJECT_PATH="$2"
            shift 2
            ;;
        --project=*)
            PROJECT_PATH="${1#*=}"
            shift
            ;;
        -h|--help)
            usage
            ;;
        *)
            echo "Unknown option: $1"
            usage
            ;;
    esac
done

# Validate script exists
if [ ! -f "$SCRIPT_TO_RUN" ]; then
    echo -e "${RED}Error: Script not found: $SCRIPT_TO_RUN${NC}"
    exit 1
fi

# Build Julia command
if [ -n "$PROJECT_PATH" ]; then
    JULIA_ARGS="--project=$PROJECT_PATH"
fi

# Format bytes to human readable
format_bytes() {
    local bytes=$1
    if [ "$bytes" -ge 1073741824 ]; then
        printf "%.2f GB" $(echo "scale=2; $bytes / 1073741824" | bc)
    elif [ "$bytes" -ge 1048576 ]; then
        printf "%.2f MB" $(echo "scale=2; $bytes / 1048576" | bc)
    elif [ "$bytes" -ge 1024 ]; then
        printf "%.2f KB" $(echo "scale=2; $bytes / 1024" | bc)
    else
        printf "%d B" "$bytes"
    fi
}

# Format memory (from KB)
format_memory() {
    local kb=$1
    local bytes=$((kb * 1024))
    format_bytes "$bytes"
}

# Check if process is running
is_running() {
    local pid=$1
    kill -0 "$pid" 2>/dev/null
    return $?
}

# Get process state description
get_state_desc() {
    local state=$1
    case "$state" in
        R*) echo "Running" ;;
        S*) echo "Sleeping (interruptible)" ;;
        D*) echo "Disk sleep (uninterruptible I/O)" ;;
        Z*) echo "Zombie" ;;
        T*) echo "Stopped" ;;
        t*) echo "Tracing stop" ;;
        W*) echo "Paging" ;;
        X*) echo "Dead" ;;
        *) echo "Unknown" ;;
    esac
}

# Continuous GPU logger (runs in background)
continuous_gpu_logger() {
    local pid=$1
    local logfile=$2
    local interval=$3
    
    # Write header
    echo "timestamp,elapsed_ms,gpu0_util,gpu0_mem_used,gpu0_mem_total,process_on_gpu,process_gpu_mem" > "$logfile"
    
    local start_time=$(date +%s%3N)  # milliseconds
    
    while kill -0 "$pid" 2>/dev/null; do
        local current_time=$(date +%s%3N)
        local elapsed=$((current_time - start_time))
        local timestamp=$(date '+%Y-%m-%d %H:%M:%S.%3N')
        
        if command -v nvidia-smi &> /dev/null; then
            # Get GPU 0 utilization and memory
            local gpu_util=$(nvidia-smi --query-gpu=utilization.gpu --format=csv,noheader,nounits -i 0 2>/dev/null | head -1)
            local gpu_mem_used=$(nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits -i 0 2>/dev/null | head -1)
            local gpu_mem_total=$(nvidia-smi --query-gpu=memory.total --format=csv,noheader,nounits -i 0 2>/dev/null | head -1)
            
            # Check if process is on GPU
            local proc_check=$(nvidia-smi --query-compute-apps=pid,used_memory --format=csv,noheader 2>/dev/null | grep "^$pid" || true)
            local on_gpu="no"
            local proc_mem="0"
            
            if [ -n "$proc_check" ]; then
                on_gpu="yes"
                proc_mem=$(echo "$proc_check" | awk -F', ' '{print $2}' | tr -d ' MiB')
            fi
            
            echo "$timestamp,$elapsed,$gpu_util,$gpu_mem_used,$gpu_mem_total,$on_gpu,$proc_mem" >> "$logfile"
        fi
        
        sleep "$interval"
    done
}

# Monitor function
monitor_process() {
    local pid=$1
    local iteration=$2
    local elapsed=$3
    
    # Get ps info
    local ps_line=$(ps -p "$pid" -o rss,vsz,pcpu,state,etime,nlwp 2>/dev/null | tail -n 1)
    
    if [ -z "$ps_line" ]; then
        echo -e "${RED}Process no longer exists${NC}"
        return 1
    fi
    
    read -r rss vsz cpu_pct state etime threads <<< "$ps_line"
    
    # Get GPU info (if nvidia-smi is available)
    local gpu_available=false
    local gpu_utilization=""
    local gpu_memory_used=""
    local gpu_memory_total=""
    local gpu_processes=""
    
    if command -v nvidia-smi &> /dev/null; then
        gpu_available=true
        # Get GPU utilization and memory for all GPUs
        gpu_utilization=$(nvidia-smi --query-gpu=utilization.gpu --format=csv,noheader,nounits 2>/dev/null | paste -sd "," -)
        gpu_memory_used=$(nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits 2>/dev/null | paste -sd "," -)
        gpu_memory_total=$(nvidia-smi --query-gpu=memory.total --format=csv,noheader,nounits 2>/dev/null | paste -sd "," -)
        
        # Check if our Julia process is using the GPU
        gpu_processes=$(nvidia-smi --query-compute-apps=pid,process_name,used_memory --format=csv,noheader 2>/dev/null | grep "^$pid" || true)
    fi
    
    # Print header
    echo ""
    echo "================================================================================"
    printf "Update #%d - Elapsed: %ds (%.2f min) - %s\n" \
        "$iteration" "$elapsed" "$(echo "scale=2; $elapsed / 60" | bc)" "$(date '+%Y-%m-%d %H:%M:%S')"
    echo "--------------------------------------------------------------------------------"
    
    # Memory info
    echo -e "${CYAN}Memory:${NC}"
    echo "  RSS (Resident):  $(format_memory $rss)"
    echo "  VSZ (Virtual):   $(format_memory $vsz)"
    
    # CPU info
    echo -e "\n${CYAN}CPU & Threads:${NC}"
    echo "  CPU Usage:       ${cpu_pct}%"
    echo "  Threads:         $threads"
    
    # Process state
    local state_desc=$(get_state_desc "$state")
    echo -e "\n${CYAN}Process State:${NC}"
    echo "  State:           $state ($state_desc)"
    echo "  Uptime:          $etime"
    
    # GPU information
    if [ "$gpu_available" = true ]; then
        echo -e "\n${CYAN}GPU Status:${NC}"
        
        if [ -n "$gpu_processes" ]; then
            echo -e "  ${GREEN}✓ GPU IS BEING USED${NC}"
            echo "  Process GPU Memory: $(echo "$gpu_processes" | awk -F', ' '{print $3}') MiB"
            
            # Parse GPU utilization and memory (comma-separated for multiple GPUs)
            IFS=',' read -ra gpu_utils <<< "$gpu_utilization"
            IFS=',' read -ra gpu_mems_used <<< "$gpu_memory_used"
            IFS=',' read -ra gpu_mems_total <<< "$gpu_memory_total"
            
            local gpu_idx=0
            for util in "${gpu_utils[@]}"; do
                echo "  GPU $gpu_idx:"
                echo "    Utilization:   ${util}%"
                echo "    Memory:        ${gpu_mems_used[$gpu_idx]} / ${gpu_mems_total[$gpu_idx]} MiB"
                gpu_idx=$((gpu_idx + 1))
            done
        else
            echo -e "  ${RED}✗ GPU NOT BEING USED${NC} (no GPU processes for this PID)"
            
            # Still show GPU status
            IFS=',' read -ra gpu_utils <<< "$gpu_utilization"
            IFS=',' read -ra gpu_mems_used <<< "$gpu_memory_used"
            IFS=',' read -ra gpu_mems_total <<< "$gpu_memory_total"
            
            local gpu_idx=0
            for util in "${gpu_utils[@]}"; do
                echo "  GPU $gpu_idx: ${util}% util, ${gpu_mems_used[$gpu_idx]} / ${gpu_mems_total[$gpu_idx]} MiB"
                gpu_idx=$((gpu_idx + 1))
            done
        fi
    fi
    
    # I/O stats from /proc
    if [ -f "/proc/$pid/io" ]; then
        local read_bytes=$(grep "^read_bytes:" /proc/$pid/io 2>/dev/null | awk '{print $2}')
        local write_bytes=$(grep "^write_bytes:" /proc/$pid/io 2>/dev/null | awk '{print $2}')
        
        if [ -n "$read_bytes" ] && [ -n "$write_bytes" ]; then
            echo -e "\n${CYAN}I/O Activity:${NC}"
            echo "  Total Read:      $(format_bytes $read_bytes)"
            echo "  Total Write:     $(format_bytes $write_bytes)"
            
            # Calculate rate if we have previous values
            if [ -n "$LAST_READ" ] && [ -n "$LAST_WRITE" ]; then
                local read_delta=$((read_bytes - LAST_READ))
                local write_delta=$((write_bytes - LAST_WRITE))
                local read_rate=$((read_delta / INTERVAL))
                local write_rate=$((write_delta / INTERVAL))
                echo "  Read rate:       $(format_bytes $read_rate)/s"
                echo "  Write rate:      $(format_bytes $write_rate)/s"
            fi
            
            # Store for next iteration
            LAST_READ=$read_bytes
            LAST_WRITE=$write_bytes
        fi
    fi
    
    # Check for precompilation (look for .ji files)
    echo -e "\n${CYAN}Precompilation Status:${NC}"
    
    local ji_files=""
    local ji_count=0
    local cuda_libs_loaded=false
    
    if command -v lsof &> /dev/null; then
        ji_files=$(lsof -p "$pid" 2>/dev/null | grep "\.ji" | awk '{print $NF}' | sort -u)
        ji_count=$(echo "$ji_files" | grep -c "\.ji" || true)
        
        # Check for CUDA libraries
        local cuda_check=$(lsof -p "$pid" 2>/dev/null | grep -i "cuda\|cublas\|cufft\|cusolver" || true)
        if [ -n "$cuda_check" ]; then
            cuda_libs_loaded=true
        fi
        
        if [ "$ji_count" -gt 0 ]; then
            echo -e "  ${YELLOW}🔄 PRECOMPILATION DETECTED${NC}"
            echo "  Cache files (.ji) open: $ji_count"
            
            local count=0
            while IFS= read -r file; do
                if [ -n "$file" ]; then
                    count=$((count + 1))
                    if [ $count -le 5 ]; then
                        echo "    [$count] $(basename "$file")"
                    fi
                fi
            done <<< "$ji_files"
            
            if [ $ji_count -gt 5 ]; then
                echo "    ... and $((ji_count - 5)) more"
            fi
        else
            echo -e "  ${GREEN}✓ Not in precompilation${NC} (no .ji cache files open)"
        fi
        
        # Show CUDA library status
        if [ "$cuda_libs_loaded" = true ]; then
            echo -e "\n  ${GREEN}✓ CUDA libraries loaded${NC}"
            local cuda_lib_count=$(lsof -p "$pid" 2>/dev/null | grep -i "cuda\|cublas\|cufft\|cusolver" | wc -l)
            echo "    CUDA-related libraries: $cuda_lib_count"
        elif [ "$gpu_available" = true ]; then
            echo -e "\n  ${YELLOW}○ CUDA libraries not loaded${NC}"
        fi
        
        # Show .jl files being accessed
        local jl_files=$(lsof -p "$pid" 2>/dev/null | grep "\.jl" | grep -v "monitor_julia" | awk '{print $NF}' | sort -u | head -n 3)
        local jl_count=$(lsof -p "$pid" 2>/dev/null | grep "\.jl" | grep -v "monitor_julia" | wc -l)
        
        if [ "$jl_count" -gt 0 ]; then
            echo -e "\n  .jl files accessed: $jl_count"
            while IFS= read -r file; do
                if [ -n "$file" ]; then
                    echo "    - $(basename "$file")"
                fi
            done <<< "$jl_files"
        fi
    else
        echo "  lsof not available - cannot detect precompilation"
    fi
    
    # Warnings
    echo ""
    local warnings=0
    
    if [[ "$state" == D* ]]; then
        echo -e "${RED}⚠️  WARNING: Process in uninterruptible sleep (likely waiting on I/O)${NC}"
        warnings=$((warnings + 1))
    fi
    
    if [ "$iteration" -gt 2 ]; then
        local cpu_int=${cpu_pct%.*}
        if [ "$cpu_int" -lt 1 ]; then
            echo -e "${YELLOW}⚠️  WARNING: Low CPU usage - process may be stuck or waiting${NC}"
            warnings=$((warnings + 1))
        fi
    fi
    
    # Check memory (16 GB = 16777216 KB)
    if [ "$rss" -gt 16777216 ]; then
        echo -e "${YELLOW}⚠️  WARNING: High memory usage (>16GB)${NC}"
        warnings=$((warnings + 1))
    fi
    
    # GPU warning: CUDA libs loaded but no GPU process detected
    if [ "$gpu_available" = true ] && [ "$cuda_libs_loaded" = true ] && [ -z "$gpu_processes" ]; then
        echo -e "${YELLOW}⚠️  WARNING: CUDA libraries loaded but GPU not actively used by process${NC}"
        echo "    This may indicate GPU code path not executing or idle state"
        warnings=$((warnings + 1))
    fi
    
    if [ $warnings -eq 0 ]; then
        echo -e "${GREEN}No warnings${NC}"
    fi
    
    echo "================================================================================"
    
    return 0
}

# Cleanup function
cleanup() {
    echo ""
    echo "================================================================================"
    echo -e "${YELLOW}Monitoring interrupted by user (Ctrl+C)${NC}"
    echo "Note: The Julia process may still be running."
    if [ -n "$JULIA_PID" ]; then
        if is_running "$JULIA_PID"; then
            echo "Julia process (PID: $JULIA_PID) is still running."
            echo "Use 'kill $JULIA_PID' to terminate it if needed."
        fi
    fi
    
    # Kill GPU logger
    if [ -n "$GPU_LOGGER_PID" ]; then
        kill "$GPU_LOGGER_PID" 2>/dev/null || true
        wait "$GPU_LOGGER_PID" 2>/dev/null || true
    fi
    
    # Show GPU log summary if it exists
    if [ -n "$GPU_LOG_FILE" ] && [ -f "$GPU_LOG_FILE" ]; then
        show_gpu_summary
    fi
    
    echo "================================================================================"
    exit 130
}

# Show GPU log summary
show_gpu_summary() {
    echo ""
    echo "================================================================================"
    echo -e "${CYAN}Continuous GPU Logging Summary${NC}"
    echo "Log file: $GPU_LOG_FILE"
    
    if [ ! -f "$GPU_LOG_FILE" ]; then
        echo "No GPU log file found."
        return
    fi
    
    local total_samples=$(tail -n +2 "$GPU_LOG_FILE" | wc -l)
    local gpu_active_samples=$(tail -n +2 "$GPU_LOG_FILE" | awk -F',' '$6=="yes"' | wc -l)
    
    if [ "$total_samples" -gt 0 ]; then
        local pct_active=$(echo "scale=2; $gpu_active_samples * 100 / $total_samples" | bc)
        
        echo "Total samples:        $total_samples"
        echo "GPU active samples:   $gpu_active_samples ($pct_active%)"
        
        if [ "$gpu_active_samples" -gt 0 ]; then
            echo -e "\n${GREEN}✓ GPU WAS USED during execution${NC}"
            
            # Get max GPU memory used by process
            local max_proc_mem=$(tail -n +2 "$GPU_LOG_FILE" | awk -F',' '{print $7}' | sort -n | tail -1)
            echo "Max process GPU mem:  ${max_proc_mem} MiB"
            
            # Get average GPU utilization when process was active
            local avg_util=$(tail -n +2 "$GPU_LOG_FILE" | awk -F',' '$6=="yes" {sum+=$3; count++} END {if(count>0) printf "%.1f", sum/count; else print "0"}')
            echo "Avg GPU util (active): ${avg_util}%"
            
            # Show time periods when GPU was active
            echo -e "\n${CYAN}GPU Active Periods:${NC}"
            tail -n +2 "$GPU_LOG_FILE" | awk -F',' '
                $6=="yes" {
                    if (!in_period) {
                        start_time=$1
                        start_elapsed=$2
                        in_period=1
                    }
                    end_time=$1
                    end_elapsed=$2
                }
                $6=="no" && in_period {
                    printf "  %s to %s (%.2f s)\n", start_time, end_time, (end_elapsed-start_elapsed)/1000
                    in_period=0
                }
                END {
                    if (in_period) {
                        printf "  %s to %s (%.2f s)\n", start_time, end_time, (end_elapsed-start_elapsed)/1000
                    }
                }
            '
        else
            echo -e "\n${RED}✗ GPU WAS NEVER USED${NC}"
            echo "The process did not appear in GPU compute apps at any sampled time."
        fi
    else
        echo "No data samples collected."
    fi
    
    echo "================================================================================"
}

# Set up signal handler
trap cleanup SIGINT SIGTERM

# Print header
echo "================================================================================"
echo "Julia Process Monitor (Bash)"
echo "Script: $SCRIPT_TO_RUN"
echo "Monitoring interval: $INTERVAL seconds"
echo "Start time: $(date '+%Y-%m-%d %H:%M:%S')"
echo "================================================================================"
echo ""

# Set up GPU logging file
if command -v nvidia-smi &> /dev/null; then
    GPU_LOG_FILE="gpu_log_$(date '+%Y%m%d_%H%M%S').csv"
    echo "Continuous GPU logging enabled (interval: ${GPU_LOG_INTERVAL}s)"
    echo "GPU log file: $GPU_LOG_FILE"
else
    echo "nvidia-smi not found - GPU logging disabled"
fi
echo ""

# Start Julia process in background
echo "Starting: julia $JULIA_ARGS $SCRIPT_TO_RUN"
echo ""

julia $JULIA_ARGS "$SCRIPT_TO_RUN" &
JULIA_PID=$!

echo "Process PID: $JULIA_PID"
echo ""

# Initial delay
sleep 2

# Check if process started successfully
if ! is_running "$JULIA_PID"; then
    echo -e "${RED}Error: Julia process failed to start or exited immediately${NC}"
    wait "$JULIA_PID"
    EXIT_CODE=$?
    echo "Exit code: $EXIT_CODE"
    exit $EXIT_CODE
fi

# Start continuous GPU logger if available
if [ -n "$GPU_LOG_FILE" ]; then
    continuous_gpu_logger "$JULIA_PID" "$GPU_LOG_FILE" "$GPU_LOG_INTERVAL" &
    GPU_LOGGER_PID=$!
    echo "Continuous GPU logger started (PID: $GPU_LOGGER_PID)"
    echo ""
fi

# Monitor loop
ITERATION=0
START_TIME=$(date +%s)
LAST_READ=""
LAST_WRITE=""

while is_running "$JULIA_PID"; do
    ITERATION=$((ITERATION + 1))
    CURRENT_TIME=$(date +%s)
    ELAPSED=$((CURRENT_TIME - START_TIME))
    
    monitor_process "$JULIA_PID" "$ITERATION" "$ELAPSED"
    
    sleep "$INTERVAL"
done

# Process finished
wait "$JULIA_PID"
EXIT_CODE=$?

# Stop GPU logger
if [ -n "$GPU_LOGGER_PID" ]; then
    kill "$GPU_LOGGER_PID" 2>/dev/null || true
    wait "$GPU_LOGGER_PID" 2>/dev/null || true
fi

CURRENT_TIME=$(date +%s)
ELAPSED=$((CURRENT_TIME - START_TIME))

echo ""
echo "================================================================================"
echo -e "${GREEN}Process completed!${NC}"
printf "Total time: %ds (%.2f min)\n" "$ELAPSED" "$(echo "scale=2; $ELAPSED / 60" | bc)"
echo "Exit code: $EXIT_CODE"
echo "End time: $(date '+%Y-%m-%d %H:%M:%S')"
echo "================================================================================"

# Show GPU summary if available
if [ -n "$GPU_LOG_FILE" ] && [ -f "$GPU_LOG_FILE" ]; then
    show_gpu_summary
fi

exit $EXIT_CODE
