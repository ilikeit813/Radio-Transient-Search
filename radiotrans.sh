#!/bin/bash
#
# radiotrans.sh
#
# Created by:     Cregg C. Yancey
# Creation date:  March 3 2018
#
# Modified by:
# Modified date:
#
# PURPOSE: Implement the full workflow documented in the readme file for generating a list of radio
#          radio transient signals from a given radio time-series data file.  Working files are dumped
#          into the current working directory, so make sure you are running in a directory that has
#          enough space.
#
# INPUT:
#     DATA_FILE: full path to the radio time-series data file.
#
# OUTPUT:
#     ppc_SNR_pol_%.1i_td_%.2i_no_%.05i.txt: radio transient file will have this pattern.
#     
# NOTE: The ${INSTALL_PATH} template variable is replaced by my git-export.sh script with the path to which
#       the package is installed.  If you don't have my git-export.sh script, then you can just manually 
#       replace all instances of ${INSTALL_PATH} as necessary.


# Make sure extended regular expressions are supported.
shopt -s extglob


# Helper functions to clean up files.
function delete_files()
{
   local FILES=("${@}")

   for curr in ${FILES[*]}
   do
      rm -f "${curr}"
   done
}
# end delete_files()
#
function transfer_results()
{
   local WORK_DIR="."
   local RESULTS_DIR="."

   if [ ${#} -gt 0]; then
      while [ -n "${1}" ]
      do
         case "${1}" in
            -w | --work-dir)
               WORK_DIR="${2}"
               shift; shift
               ;;
            -r | --results-dir)
               RESULTS_DIR="${2}"
               shift; shift
               ;;
            *)
               shift
               ;;
         esac
      done
   fi

   if [ "${WORK_DIR} != "${RESULTS_DIR} ]; then
      echo "Transferring results from ${WORK_DIR} to ${RESULTS_DIR}..."
      mv "${WORK_DIR}/lowspectrogram.png" "${RESULTS_DIR}/lowspectrogram.png"
      mv "${WORK_DIR}/highspectrogram.png" "${RESULTS_DIR}/highspectrogram.png"
      mv "${WORK_DIR}/lowbandpass" "${RESULTS_DIR}/lowbandpass"
      mv "${WORK_DIR}/highbandpass.png" "${RESULTS_DIR}/highbandpass.png"
      mv "${WORK_DIR}/*.txt" "${RESULTS_DIR}/*.txt"
   fi
}
# end transfer_results()



# Set the default install path.  NOTE: the string 'OPT-INSTALL_DIR' is replaced by the git-export.sh
# script with the path to directory in which radiotrans has been installed.
INSTALL_PATH="OPT-INSTALL_DIR"
if [[ ${INSTALL_PATH} == "OPT-INSTALL_DIR" ]]; then
   INSTALL_PATH="."
fi


USAGE='
radiotrans.sh

   radiotrans.sh [-d | --install-dir <path>] [-i | --integrate-time <secs>] <data_file>

   Implement the full workflow documented in the readme file for generating a list of radio
   radio transient signals from a given radio time-series data file.  Working files are dumped
   into the current working directory, so make sure you are running in a directory that has
   enough space.

   ARGUMENTS:
      <data_file>: path to the radio data file.

   OPTIONS:
      -h | --help:                  Display this help message.

      -d | --install-dir <path>:    Specify <path> as the path to the directory containing the
                                    radiotrans component modules.

      -i | --integrate-time <secs>: Specify <secs> as the spectral integration time, in seconds for each
                                    spectral sample.

      -s | --samples <num>:         Specify the number of spectral samples to extract from the radio
                                    data file.

      -p | --samples-per-sec <num>: Specify the number of spectral samples per second of data to extract
                                    from the radio data file.  This overrides the specification for the
                                    total number of samples.

      -n | --nprocs <num>:          Number of MPI processes to use.  This defaults to the number of
                                    processor cores on the machine.

      -m | --memory-limit <num>:    Memory limit, in MBs, per MPI process during the de-dispersive
                                    search.

      --dm-start <num>:             Starting dispersion-measure trial value.

      --dm-end <num>:               Ending dispersion-measure trial value.
'

AFFIRMATIVE='^(y|yes|yup|yea|yeah|ya)$'


# ==== MAIN WORKFLOW FOR RADIOTRANS.SH ===
#
#
RESUME_CMD_FILEPATH="./radiotrans_cmd.resume"
RESUME_VAR_FILEPATH="./radiotrans_var.resume"

DATA_PATH=        # Path to the radio time-series data file.
SPECTINTEGTIME=5  # Spectral integration time.
NUM_SAMPLES=10000 # Number of spectral samples to extract.
NUM_SAMPLESPERSEC=0  # Number of spectral samples per second to extract in the waterfall.  If this is
                     # greater than zero, then use it to override the total number of samples.
                     
NUM_PROCS=$(nproc --all)   # Number of concurrent processes to use under MPI
MEM_LIMIT=16      # Memory limit, in MBs, per MPI process when performing the de-dispersive search.
DM_START=300      # Starting dispersion-measure trial value.
DM_END=2000       # Ending dispersion-measure trial value.
WORK_DIR="."      # Working directory.
RESULTS_DIR="."   # Results directory.


# Parse command-line arguments.
if [[ ${#} -gt 0 ]]; then
   while [ -n "${1}" -a -z "${DATA_PATH}" ]
   do
      case "${1}" in
         -h | --help) # Display the help message and then quit.
            echo "${USAGE}"
            exit 0
            ;;
         -d | --install-dir) # Specify the install directory path to the radiotrans package.
            INSTALL_PATH="${2}"
            shift; shift
            ;;
         -i | --integrate-time) # Specify the spectral integration time.
            SPECTINTEGTIME=${2}
            # CCY - TODO: add checking that this is an actual number.
            shift; shift
            ;;
         -s | --samples) # Specify the number of spectral samples to extract for the coarse spectrogram.
            NUM_SAMPLES=${2}
            shift; shift
            ;;
         -p | --samples-per-sec) # Override the total sample count with the specified number of samples
                                 # per second.
            if [ ${2} -gt 0 ]; then
               NUM_SAMPLESPERSEC=${2}
            fi
            shift; shift
            ;;
         -n | --nprocs) # Specify the number of MPI processes to use
            NUM_PROCS=${2}
            # CCY - TODO: add checking that this is an actual number.
            if [ ${NUM_PROCS} -lt 1 ]; then
               NUM_PROC=1
            fi
            shift; shift
            ;;
         -m | --memory-limit) # Specify the memory limit when performing the de-dispersive search.
            MEM_LIMIT=${2}
            # CCY - TODO: add checking that this is an actual number.
            if [ ${MEM_LIMIT} -lt 1 ]; then
               MEM_LIMIT=1
            fi
            shift; shift
            ;;
         --dm-start) # Specify the starting dispersion-measure trial value.
            DM_START=${2}
            shift; shift
            ;;
         --dm-end) # Specify the ending dispersion-measure trial value.
            DM_END=${2}
            shift; shift
            ;;
         --work-dir) # Specify the working directory.
            WORK_DIR="${2}"
            shift; shift
            ;;
         --results-dir) # Specify the results directory.
            RESULTS_DIR="${2}"
            shift; shift
            ;;
         -*) # Unknown option
            echo "ERROR: radiotrans.sh -> Unknown option"
            echo "     ${1}"
            exit 1
            ;;
         *) # Get the data file path
            DATA_PATH="${1}"
            shift
            ;;
      esac
   done
else
   echo "WARNING: radiotrans -> Nothing specified to do"
   echo "${USAGE}"
   exit 0
fi

# Check that a valid data file path has been specified.
if [ -n "${DATA_PATH}" ]; then
   if [ ! -f "${DATA_PATH}" ]; then
      echo "ERROR: radiotrans.sh -> Data file path not found"
      echo "     ${DATA_PATH}"
      exit 1
   fi
else
   echo "ERROR: radiotrans.sh -> Must provide a path to the data file."
   exit 1
fi

# Check that specified install path exists and that all necessary components are contained.
if [ -d "${INSTALL_PATH}" ]; then
   package_modules=(dp.py drx.py waterfall.py chkwaterfall.py eyexam.py interpolate.py \
                     waterfallcombine.py bandpasscheck.py ft.py freqtint.py dv.py disper.py \
                     plot.py errors.py apputils.py)
   for module in ${package_modules[*]}; do
      MODULE_PATH="${INSTALL_PATH}/${module}"
      if [ ! -f "${MODULE_PATH}" ]; then
         echo "ERROR: radiotrans.sh -> Missing package module"
         echo "     ${MODULE_PATH}"
         exit 1
      fi
   done
else
   echo "ERROR: radiotrans.sh -> Install path does not exist"
   echo "     ${INSTALL_PATH}"
   exit 1
fi

# Check that the working directory exists.
if [ ! -d "${WORK_DIR}" ]; then
   echo "ERROR: radiotrans.sh -> working directory does not exist.  User may need to create it."
   echo "     ${WORK_DIR}"
   exit 1
fi
# Check that the results directory exists.
if [ ! -d "${RESULTS_DIR}" ]; then
   echo "ERROR: radiotrans.sh -> results directory does not exist.  User may need to create it."
   echo "     ${RESULTS_DIR}"
   exit 1
fi


# Source the resume functionality.
source ${INSTALL_PATH}/resume.sh


# IMPLEMENT RESUMABLE RADIOTRANS WORKFLOW 
#
echo "Starting radiotrans workflow:"
# Workflow resume labels.  These are to label each executable stage of the workflow for use with
# resumecmd.
#
LBL_WATERFALL1="WaterfallCoarse"
LBL_WATERFALL2="WaterfallDetailed"
LBL_CHKWATERFALLA="CheckWaterfallA"
LBL_CHKWATERFALLB="CheckWaterfallB"
LBL_EYEEXAM="Eye-exam"
LBL_INTERPOLATE="Interpolate"
LBL_COMBINE="WaterfallCombine"
LBL_BANDPASS="Bandpass"
LBL_SPECTROGRAM="Spectrogram"
LBL_FRAMETRIM="FrameTrim"
LBL_FREQTINT="FrequencyTimeIntegration"
LBL_RFICUT="SmoothRFIBandpass"
LBL_DEDISPERSLOW="De-dispersion_Low"
LBL_DEDISPERSHIGH="De-dispersion_High"
LBL_RESULTS="TransferResults"
LBL_CLEAN1="Cleanup1"
LBL_CLEAN2="Cleanup2"
LBL_CLEAN2B="Cleanup2B"
LBL_CLEAN3="Cleanup3"

# The calls to the individual python scripts are patterned after waterfall.sh.
#


# Generate the coarse spectrogram
echo "    Generating coarse spectral samples for coarse spectrogram from ${DATA_PATH}..."
resumecmd -l ${LBL_WATERFALL1} mpirun -np ${NUM_PROCS} python ${INSTALL_PATH}/waterfall.py \
   --integrate-time ${SPECTINTEGTIME} --samples ${NUM_SAMPLES} \
   --samples-per-sec ${NUM_SAMPLESPERSEC} --work-dir ${WORK_DIR} "${DATA_PATH}"
report_resumecmd


# Combine the individual waterfall spectral sample files into a singular coarse spectrogram file.
echo "    Combining coarse spectral sample files into singular coarse spectrogram file..."
COMBWATERFALL="combined_waterfall"
resumecmd -l ${LBL_COMBINE} -k ${RESUME_LASTCMD_SUCCESS} \
   mpirun -np 1 python ${INSTALL_PATH}/waterfallcombine.py \
   --outfile "${WORK_DIR}/${COMBWATERFALL}" --work-dir "${WORK_DIR}" "${DATA_PATH}"
report_resumecmd


echo "Cleaning up coarse spectrogram files (this may take a few minutes)..."
DATAFILENAME=$(basename ${DATA_PATH})
resumecmd -l ${LBL_CLEAN1} -k ${RESUME_LASTCMD_SUCCESS} \
   delete_files "${WORK_DIR}/waterfall${DATAFILENAME%.*}*.npy"
report_resumecmd


# View the spectrograms for the low and high frequency tunings and adjust the bandpass filtration for
# each tuning according to the user's specification.
LBL_LOWFCL="LOWFCL"
LBL_LOWFCH="LOWFCH"
LBL_HIGHFCL="HIGHFCL"
LBL_HIGHFCH="HIGHFCH"
resumevar -l ${LBL_LOWFCL} LOW_FCL 0         # Lower FFT index for the low tuning frequency.
resumevar -l ${LBL_LOWFCH} LOW_FCH 4095      # Upper FFT index for the low tuning frequency.
resumevar -l ${LBL_HIGHFCL} HIGH_FCL 0       # Lower FFT index for the high tuning frequency.
resumevar -l ${LBL_HIGHFCH} HIGH_FCH 4095    # Upper FFT index for the high tuning frequency.

if [[ ${RESUME_LASTCMD_SUCCESS} -ne 0 ]]; then
   FLAG_CONTINUE=1
   # Before we try to view the bandpass spectrograms, we need to pause the workflow here because the user
   # may or may not be online when we get to this point, and the user may or may not have an X11 server and
   # X11 forwarding setup.  Pausing here will give the user a chance to get online and perform any
   # necessary setup to view the bandpass plots before moving to the next phase.
   echo
   echo "    The next phase involves adjusting the bandpass filtering for RFI.  Plots of the bandpass"
   echo "    intensity and the spectrogram are generated, and the user can elect to have the plots be"
   echo "    displayed from within the workflow (requires X11 and X11 forwarding, if working remotely);"
   echo "    however, if this workflow is being run in the background with the user not in attendance,"
   echo "    it is recommended that the plots not be displayed through the workflow and that the user"
   echo "    use a utility like 'display' or 'imagemagick' to display the plot files (lowbandpass.png,"
   echo "    highbandpass.png, lowspectrogram.png, highspectrogram.png"
   echo
   echo "    Display plots from within the workflow (requires X11 and X11 forwarding)? [y/n] "
   read USER_ANSWER
   USE_PLOTS=$(echo "${USER_ANSWER}" | tr '[:upper:]' '[:lower:]')
   RESUME_FORCE_OPT=
   while [[ ${FLAG_CONTINUE} -eq 1 ]]
   do
      FLAG_CONTINUE=0
      echo
      echo "     Current FFT bandpass settings:"
      echo "         Low tuning lower index = ${LOW_FCL}"
      echo "         Low tuning upper index = ${LOW_FCH}"
      echo "         High tuning lower index = ${HIGH_FCL}"
      echo "         High tuning upper index = ${HIGH_FCH}"
      echo
      echo "    Generate plots to visualize bandpass? [y/n]"
      read USER_ANSWER
      GEN_PLOTS=$(echo "${USER_ANSWER}" | tr '[:upper:]' '[:lower:]')

      # Generate feedback plots if selected by user.
      if [[ "${GEN_PLOTS}" =~ ${AFFIRMATIVE} ]]; then
         echo "    Generating bandpass plots..."
         resumecmd -l ${LBL_BANDPASS} -k ${RESUME_LASTCMD_SUCCESS} -s ${RESUME_REPEAT} \
            mpirun -np 1 python ${INSTALL_PATH}/bandpasscheck.py --low-tuning-lower ${LOW_FCL} \
            --low-tuning-upper ${LOW_FCH} --high-tuning-lower ${HIGH_FCL} \
            --high-tuning-upper ${HIGH_FCH} --work-dir ${WORK_DIR} "${WORK_DIR}/${COMBWATERFALL}.npy"
         report_resumecmd

         echo "    Generating spectrogram plots..."
         resumecmd -l ${LBL_SPECTROGRAM} -k ${RESUME_LASTCMD_SUCCESS} -s ${RESUME_REPEAT} \
            mpirun -np 1 python ${INSTALL_PATH}/watchwaterfall.py --low-tuning-lower ${LOW_FCL} \
            --low-tuning-upper ${LOW_FCH} --high-tuning-lower ${HIGH_FCL} \
            --high-tuning-upper ${HIGH_FCH} --work-dir ${WORK_DIR} "${WORK_DIR}/${COMBWATERFALL}.npy"
         report_resumecmd

         # Display the plots if selected by the user.
         if [[ "${USE_PLOTS}" =~ ${AFFIRMATIVE} ]]; then 
            display "${WORK_DIR}/lowbandpass.png" &
            display "${WORK_DIR}/highbandpass.png" &

            display "${WORK_DIR}/lowspectrogram.png" &
            display "${WORK_DIR}/highspectrogram.png" &
         fi
      fi

      echo "     Do you wish to change the bandpass region? [y/n] "
      read USER_ANSWER
      USER_ANSWER=$(echo "${USER_ANSWER}" | tr '[:upper:]' '[:lower:]')
      if [[ "${USER_ANSWER}" =~ ${AFFIRMATIVE} ]]; then 
         echo "    Change low tuning bandpass? [y/n]"
         read USER_ANSWER
         USER_ANSWER=$(echo "${USER_ANSWER}" | tr '[:upper:]' '[:lower:]')
         if [[ "${USER_ANSWER}" =~ ${AFFIRMATIVE} ]]; then 
            echo "    Enter the lower FFT index for the low tuning (integer between 0 and 4095): "
            read USER_ANSWER
            resumevar -l ${LBL_LOWFCL} -f LOW_FCL ${USER_ANSWER}
            echo "    Enter the upper FFT index for the low tuning (integer between 0 and 4095): "
            read USER_ANSWER
            resumevar -l ${LBL_LOWFCH} -f LOW_FCH ${USER_ANSWER}
            FLAG_CONTINUE=1
            RESUME_FORCE_OPT="-f"
         fi

         # Change bandpass for high tuning, if desired.
         echo "    Change high tuning bandpass? [y/n]"
         read USER_ANSWER
         USER_ANSWER=$(echo "${USER_ANSWER}" | tr '[:upper:]' '[:lower:]')
         if [[ "${USER_ANSWER}" =~ ${AFFIRMATIVE} ]]; then 
            echo "    Enter the lower FFT index for the high tuning (integer between 0 and 4095): "
            read USER_ANSWER
            resumevar -l ${LBL_HIGHFCL} -f HIGH_FCL ${USER_ANSWER}
            echo "    Enter the upper FFT index for the high tuning (integer between 0 and 4095): "
            read USER_ANSWER
            resumevar -l ${LBL_HIGHFCH} -f HIGH_FCH ${USER_ANSWER}
            FLAG_CONTINUE=1
            RESUME_FORCE_OPT="-f"
         fi
      fi
   done
fi



# Extract information about the time-series for use in doing the de-dispersion and transient extraction
# with dv.py
echo "     Generating time-series information for de-dispersion..."
resumecmd -l ${LBL_FREQTINT} ${RESUME_FORCE_OPT} -k ${RESUME_LASTCMD_SUCCESS} \
   mpirun -np 1 python ${INSTALL_PATH}/freqtint.py "${DATA_PATH}" \
   --low-tuning-lower ${LOW_FCL} --low-tuning-upper ${LOW_FCH} \
   --high-tuning-lower ${HIGH_FCL} --high-tuning-upper ${HIGH_FCH} --work-dir ${WORK_DIR}
report_resumecmd


# Generate the detailed spectrogram.
echo "    Generating detailed spectral samples for de-dispersion from ${DATA_PATH}..."
resumecmd -l ${LBL_WATERFALL2} ${RESUME_FORCE_OPT} -k ${RESUME_LASTCMD_SUCCESS} \
   mpirun -np ${NUM_PROCS} python ${INSTALL_PATH}/waterfall.py \
   --integrate-time ${SPECTINTEGTIME} --samples ${NUM_SAMPLES} \
   --samples-per-sec ${NUM_SAMPLESPERSEC} --detailed --work-dir ${WORK_DIR} "${DATA_PATH}"
report_resumecmd


# Perform data smoothing, RFI cleaning, and bandpass filtering of detailed spectrogram.
echo "    Performing data smoothing, RFI cleaning, and bandpass filtering of detailed spectrogram..."
resumecmd -l ${LBL_RFICUT} ${RESUME_FORCE_OPT} -k ${RESUME_LASTCMD_SUCCESS} \
   mpirun -np ${NUM_PROCS} python ${INSTALL_PATH}/rficut.py "${DATA_PATH}" \
   --low-tuning-lower ${LOW_FCL} --low-tuning-upper ${LOW_FCH} \
   --high-tuning-lower ${HIGH_FCL} --high-tuning-upper ${HIGH_FCH} \
   --bandpass-window 10 --baseline-window 50 --work-dir ${WORK_DIR}
report_resumecmd

#PROFILE="${RESULTS_DIR}/rficut_profile"
#resumecmd -l ${LBL_RFICUT} ${RESUME_FORCE_OPT} -k ${RESUME_LASTCMD_SUCCESS} \
#   mpirun -np 1 python -m cProfile -o "${PROFILE}" ${INSTALL_PATH}/rficut.py "${DATA_PATH}" \
#   --low-tuning-lower ${LOW_FCL} --low-tuning-upper ${LOW_FCH} \
#   --high-tuning-lower ${HIGH_FCL} --high-tuning-upper ${HIGH_FCH} \
#   --bandpass-window 10 --baseline-window 50 --work-dir ${WORK_DIR}
#report_resumecmd


# Perform the de-dispersion for each of the tunings.
echo "     Performing de-dispersion on low tuning data..."
#NUM_PROCS=1
resumecmd -l ${LBL_DEDISPERSLOW} ${RESUME_FORCE_OPT} -k ${RESUME_LASTCMD_SUCCESS} \
   mpirun -np ${NUM_PROCS} python ${INSTALL_PATH}/dv.py "${DATA_PATH}" \
   --lower ${LOW_FCL} --upper ${LOW_FCH} --tuning 0 --frequency-file "${WORK_DIR}/lowtunefreq.npy" \
   --integration-time ${SPECTINTEGTIME} --memory-limit ${MEM_LIMIT} \
   --samples-per-sec ${NUM_SAMPLESPERSEC} --dm-start ${DM_START} --dm-end ${DM_END} \
   --work-dir ${WORK_DIR}
report_resumecmd

exit 1 # Debugging short-circuit

echo "     Performing de-dispersion on high tuning data..."
resumecmd -l ${LBL_DEDISPERSHIGH} ${RESUME_FORCE_OPT} -k ${RESUME_LASTCMD_SUCCESS} \
   mpirun -np ${NUM_PROCS} python ${INSTALL_PATH}/dv.py "${DATA_PATH}" \
   --lower ${HIGH_FCL} --upper ${HIGH_FCH} --tuning 1 --frequency-file "${WORK_DIR}/hightunefreq.npy" \
   --integration-time ${SPECTINTEGTIME} --memory-limit ${MEM_LIMIT} \
   --samples-per-sec ${NUM_SAMPLESPERSEC} --dm-start ${DM_START} --dm-end ${DM_END} \
   --work-dir ${WORK_DIR}
report_resumecmd


# Cleanup the detailed spectrogram files.  We don't want to leave these hanging around cause they take
# up a lot of disk space.  Also, cleanup an other intermediate files that we no longer need.
echo "Cleaning up detailed spectrogram files (this may take quite a while)..."
DATAFILENAME=$(basename ${DATA_PATH})
resumecmd -l ${LBL_CLEAN2} -k ${RESUME_LASTCMD_SUCCESS} ${RESUME_FORCE_OPT} \
   delete_files "${WORK_DIR}/master0_${DATAFILENAME%.*}*.npy"
resumecmd -l ${LBL_CLEAN2B} -k ${RESUME_LASTCMD_SUCCESS} ${RESUME_FORCE_OPT} \
   delete_files "${WORK_DIR}/master1_${DATAFILENAME%.*}*.npy"
report_resumecmd
echo "Cleaning up other files..."
DATAFILENAME=$(basename ${DATA_PATH})
resumecmd -l ${LBL_CLEAN3} -k ${RESUME_LASTCMD_SUCCESS} ${RESUME_FORCE_OPT} \
   rm -f "${WORK_DIR}/lowtunefreq.npy" "${WORK_DIR}/hightunefreq.npy" \
         "${WORK_DIR}/tint.npy" "${WORK_DIR}/${COMBWATERFALL}.npy"
report_resumecmd


# Transfer results to the results directory.
resumecmd -l ${LBL_CLEAN3} -k ${RESUME_LASTCMD_SUCCESS} ${RESUME_FORCE_OPT} \
   transfer_results --work-dir "${WORK_DIR}" --results-dir "${RESULTS_DIR}"
report_resumecmd

echo "Radio transient workflow done!"
exit 0

