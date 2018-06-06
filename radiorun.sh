#!/bin/bash

# Configure run parameters.
DATA_DIR="/data/network/recent_data/jtsai"   # Directory path containing radio data file.
DATA_FILENAME="057974_001488582"             # Radio data filename.
DATA_PATH="${DATA_DIR}/${DATA_FILENAME}"     # Full path to the radio data file.

INTEGTIME=0.025            # Spectral integration time in seconds.
NUM_SAMPLES=4000           # Number of spectral samples to extract from the radio data file.
NUM_SAMPLESPERSEC=10       # Number of spectral samples per second of data to extract.  This overrides
                           # the extraction of a specified total number of samples if the value is
                           # greater than zero.
DM_START=100              # Starting DM for trials.
DM_END=200                # Ending DM for trials.

NUM_PROCS=$(nproc --all)   # Number of processes to use in the MPI environment.
MEM_LIMIT=4096            # Memory limit, in MBs, for loading files during the de-dispersion search
                           # stage.

WORK_DIR="/mnt/toaster/cyancey"
RESULTS_DIR="."


# Select whether we are using the release install of radiotrans or still using the developer version to
# debug issues.
INSTALL_DIR="${HOME}/local/radiotrans"
if [[ ${#} -gt 0 ]]; then
   while [ -n "${1}" ]
   do
      if [[ "${1}" =~ (-D|--DEBUG) ]]; then
         INSTALL_DIR="${HOME}/dev/radiotrans"
         break
      fi
      shift
   done
fi

# Create the working directory, if it doesn't exist.
if [ ! -d "${WORK_DIR}" ]; then
   mkdir "${WORK_DIR}"
   if [ ! -d "${WORK_DIR}" ]; then
      echo "Could not create working directory ${WORK_DIR}"
      exit 1
   fi
fi

# Build the command-line to run radiotrans.sh
CMD="${INSTALL_DIR}/radiotrans.sh"
CMD_OPTS=(--install-dir "${INSTALL_DIR}" --integrate-time ${INTEGTIME} --samples ${NUM_SAMPLES} \
      --samples-per-sec ${NUM_SAMPLESPERSEC} --nprocs ${NUM_PROCS} --memory-limit ${MEM_LIMIT} \
      --work-dir "${WORK_DIR}" --results-dir "${RESULTS_DIR}" --dm-start ${DM_START} \
      --dm-end ${DM_END})

# Run the radiotrans.sh script.
${CMD} ${CMD_OPTS[*]} "${DATA_PATH}"

exit 0
