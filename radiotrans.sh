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
                                    data file to generate the coarse spectrogram.

'

AFFIRMATIVE='^(y|yes|yup|yea|yeah|ya)$'


# Source the resume functionality.
source ${INSTALL_PATH}/resume.sh


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
NUM_PROCS=1       # Number of concurrent processes to use under MPI


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


# IMPLEMENT RESUMABLE RADIOTRANS WORKFLOW 
#
echo "Starting radiotrans workflow:"
# Workflow resume labels.  These are to label each executable stage of the workflow for use with
# resumecmd.
#
LBL_WATERFALL="Waterfall"
LBL_CHKWATERFALLA="CheckWaterfallA"
LBL_CHKWATERFALLB="CheckWaterfallB"
LBL_EYEEXAM="Eye-exam"
LBL_INTERPOLATE="Interpolate"
LBL_COMBINE="WaterfallCombine"
LBL_BANDPASS="Bandpass"
LBL_FRAMETRIM="FrameTrim"
LBL_FREQTINT="FrequencyTimeIntegration"
LBL_DEDISPERSION="De-dispersion"

# The calls to the individual python scripts are patterned after waterfall.sh.
#


# Generate spectral samples for the coarse spectrogram.
echo "    Generating spectral samples for coarse spectrogram from ${DATA_PATH}..."
resumecmd -l ${LBL_WATERFALL} mpirun -np ${NUM_PROCS} python ${INSTALL_PATH}/waterfall.py \
   --integrate-time ${SPECTINTEGTIME} --samples ${NUM_SAMPLES} \
   --samples-per-sec ${NUM_SAMPLESPERSEC} "${DATA_PATH}"


# CCY - The original readme for the workflow claims this checks through the waterfall spectral samples and
# determines if any are missing.  In the computing environment in which this workflow will be used, that
# is not the case.  Instead, what does is insert a selection of additonal spectral samples that were not
# extracted by waterfall.py.  These additional samples, of course, exist between the ones generated by
# waterfall.py.  However, this same result can now be achieved from waterfall.py by simply specifying a
# higher number of spectral samples to generate.  Therefore, this stage does not seem really necessary.
# Consequently, for now, until I can better determine its function, this stage will be omitted from the
# workflow.
#
#echo "    Generating additional intermediate waterfall spectral samples..."
#resumecmd -l ${LBL_CHKWATERFALLA} -k ${RESUME_LASTCMD_SUCCESS} \
#   mpirun -np ${NUM_PROCS} python ${INSTALL_PATH}/chkwaterfall.py \
#   --integrate-time ${SPECTINTEGTIME} --samples ${NUM_SAMPLES} "${DATA_PATH}"
#exit 0


# CCY - Removed this entire stage from the workflow.  Deep examination of the code lead me to the
# conclusion that this effort is not necessary.  However, I'm leaving the code here in case it does
# become necessary to re-introduce this stage, hopefully in an improved form.
#
#FLAG_CONTINUE=1
#RESUME_FORCE_OPT=
#while [[ ${FLAG_CONTINUE} -eq 1 ]] && [[ ${RESUME_LASTCMD_SUCCESS} -ne 0 ]]
#do
#   resumecmd -l ${LBL_CHKWATERFALLB} ${RESUME_FORCE_OPT} -k ${RESUME_LASTCMD_SUCCESS} \
#      mpirun -np ${NUM_PROCS} python ${INSTALL_PATH}/chkwaterfall.py "${DATA_PATH}"
#   resumecmd -l ${LBL_EYEEXAM} ${RESUME_FORCE_OPT} -k ${RESUME_LASTCMD_SUCCESS} \
#      mpirun -np ${NUM_PROCS} python ${INSTALL_PATH}/eyexam.py "${DATA_PATH}"
#   resumecmd ${RESUME_FORCE_OPT} -k ${RESUME_LASTCMD_SUCCESS} \
#      echo "Are there any missing frame files?[y/n]"
#   resumecmd ${RESUME_FORCE_OPT} -k ${RESUME_LASTCMD_SUCCESS} read USER_ANSWER
#   # Newer versions of bash have a more elegant way of converting strings to all lower or upper case,
#   # however, I'm doing it this way for compatibility.
#   resumecmd ${RESUME_FORCE_OPT} -k ${RESUME_LASTCMD_SUCCESS} \
#      USER_ANSWER=$(echo "${USER_ANSWER}" | tr '[:upper:]' '[:lower:]')
#   if [[ "${USER_ANSWER}" =~ ${AFFIRMATIVE} ]]; then
#      echo "    Interpolating missing frame files..."
#      resumecmd -l ${LBL_INTERPOLATE} ${RESUME_FORCE_OPT} -k ${RESUME_LASTCMD_SUCCESS} \
#         mpirun -np ${NUM_PROCS} python ${INSTALL_PATH}/interpolate.py "${DATA_PATH}"
#      # Force resumecmd to execute commands in the loop because we need to do all this again.
#      RESUME_FORCE_OPT="-f"
#   else
#      FLAG_CONTINUE=0
#   fi
#done


# Combine the individual spectral sample files into a singular coarse spectrogram file.
echo "    Combining spectral sample files into singular coarse spectrogram file..."
COMBWATERFALL="combined_waterfall"
resumecmd -l ${LBL_COMBINE} -k ${RESUME_LASTCMD_SUCCESS} \
   mpirun -np 1 python ${INSTALL_PATH}/waterfallcombine.py \
   --outfile "${COMBWATERFALL}" "${DATA_PATH}"


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

      # Generate feedback plots if selected by user.
      echo "    Generating bandpass plots..."
      resumecmd -l ${LBL_BANDPASS} -k ${RESUME_LASTCMD_SUCCESS} -s ${RESUME_REPEAT} \
         mpirun -np 1 python ${INSTALL_PATH}/bandpasscheck.py --low-tuning-lower ${LOW_FCL} \
         --low-tuning-upper ${LOW_FCH} --high-tuning-lower ${HIGH_FCL} \
         --high-tuning-upper ${HIGH_FCH} "${COMBWATERFALL}.npy"

      echo "    Generating spectrogram plots..."
      resumecmd -l ${LBL_BANDPASS} -k ${RESUME_LASTCMD_SUCCESS} -s ${RESUME_REPEAT} \
         mpirun -np 1 python ${INSTALL_PATH}/watchwaterfall.py --low-tuning-lower ${LOW_FCL} \
         --low-tuning-upper ${LOW_FCH} --high-tuning-lower ${HIGH_FCL} \
         --high-tuning-upper ${HIGH_FCH} "${COMBWATERFALL}.npy"

      # Display the plots if selected by the user.
      if [[ "${USE_PLOTS}" =~ ${AFFIRMATIVE} ]]; then 
         display lowbandpass.png &
         display highbandpass.png &

         display lowspectrogram.png &
         display highspectrogram.png &
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

# CCY - The stage with ft.py, from examining the code, is exactly the same as waterfall.py except that
# it trims the resultant waterfall to the bandpass regions that we have specified.  However, some errors
# were discovered in how it goes about this.  Essentially, when the code creates the internal
# n-dimensional array to hold the spectrogram for the tunnings, it sizes the number of FFT bins by the
# number used for the low-tuning bandpass.  It uses this for BOTH the low and high tuning FFT
# spectrograms.  But, by definition of the operation, we are allow to select very different spans in the
# FFT indices for the low tuning vs. the high tuning.  Thus, there is no guarantee that the number of
# FFT bins used for the low tuning are greater than or equal to the number of bins for the high tuning.
# Thus, the results are not necessarily valid because information could be lost.
#
# Regardless of this problem, this is actually an unnecessary step because we have already validated the
# spectrogram using the improved functionality in bandpasscheck.py and waterwaterfall.py.
# 
# Generate new frame files using the bandpass regions specified above.
#echo "     Generating new frame files using specified bandpass regions..."
#resumecmd -l ${LBL_FRAMETRIM} ${RESUME_FORCE_OPT} -k ${RESUME_LASTCMD_SUCCESS} \
#   mpirun -np ${NUM_PROCS} python ${INSTALL_PATH}/ft.py "${DATA_PATH}" \
#   --low-tuning-lower ${LOW_FCL} \
#   --low-tuning-upper ${LOW_FCH} --high-tuning-lower ${HIGH_FCL} --high-tuning-upper ${HIGH_FCH}



# Extract information about the time-series for use in doing the de-dispersion and transient extraction
# with dv.py
echo "     Generating time-series information for de-dispersion..."
resumecmd -l ${LBL_FREQTINT} ${RESUME_FORCE_OPT} -k ${RESUME_LASTCMD_SUCCESS} \
   mpirun -np 1 python ${INSTALL_PATH}/freqtint.py "${DATA_PATH}" \
   --low-tuning-lower ${LOW_FCL} \
   --low-tuning-upper ${LOW_FCH} --high-tuning-lower ${HIGH_FCL} --high-tuning-upper ${HIGH_FCH}

# CCY - Stop workflow here for debugging.
exit 0


# Perform the de-dispersion for each of the tunings.
echo "     Performing de-dispersion on low tuning data..."
resumecmd -l ${LBL_DEDISPERSIONA} ${RESUME_FORCE_OPT} -k ${RESUME_LASTCMD_SUCCESS} \
   mpirun -np ${NUM_PROCS} python ${INSTALL_PATH}/dv.py "${DATA_PATH}" \
   --lower ${LOW_FCL} --upper ${LOW_FCH} --tuning "low"
echo "     Performing de-dispersion on high tuning data..."
resumecmd -l ${LBL_DEDISPERSIONB} ${RESUME_FORCE_OPT} -k ${RESUME_LASTCMD_SUCCESS} \
   mpirun -np ${NUM_PROCS} python ${INSTALL_PATH}/dv.py "${DATA_PATH}" \
   --lower ${HIGH_FCL} --upper ${HIGH_FCH} --tuning "high"


# Handle cleanup, if everything went okay.  We will no longer need all those frame files and other
# temporary work files hanging around.
echo "Performing file cleanup..."
DATAFILENAME=$(basename ${DATA_PATH})
resumecmd -l "CLEANUP1" -k ${RESUME_LASTCMD_SUCCESS} rm -f "./waterfall${DATAFILENAME%.*}*.npy"

echo "Radio transient workflow done!"
exit 0

