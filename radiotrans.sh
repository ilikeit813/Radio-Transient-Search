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
# NOTE: The OPT-INSTALL_DIR template variable is replaced by my git-export.sh script with the path to which
#       the package is installed.  If you don't have my git-export.sh script, then you can just manually 
#       replace all instances of OPT-INSTALL_DIR as necessary.

USAGE='
radiotrans.sh

   radiotrans.sh <data_file>

   Implement the full workflow documented in the readme file for generating a list of radio
   radio transient signals from a given radio time-series data file.  Working files are dumped
   into the current working directory, so make sure you are running in a directory that has
   enough space.

   ARGUMENTS:
      <data_file>: full path to the data file.

   OPTIONS:
      -h | --help: Display this help message.

      '
AFFIRMATIVE='^(y|yes|yup|yea|yeah|ya)$'

DATA_PATH=        # Path to the radio time-series data file.


# Parse command-line arguments.
while [ -n "${1}" -a -z "${DATA_PATH}" ]
do
   case "${1}" in
      -h | --help) # Display the help message and then quit.
         echo "${USAGE}"
         exit 0
         ;;
      -*) # Unknown option
         echo "UNKNOWN OPTION: ${1}"
         exit 1
         ;;
      *) # Get the data file path
         DATA_PATH="${1}"
         shift
         ;;
   esac
done

# Check that a valid data file path has been specified.
if [ -n "${DATA_PATH}" ]; then
   if [ ! -f "${DATA_PATH}" ]; then
      echo "INVALID DATA FILE PATH: ${DATA_PATH}"
      exit 1
   fi
else
   echo "Must provide a path to the data file."
   exit 1
fi




# == Implement the radio transient workflow ==
# The calls to the individual python scripts are patterned after waterfall.sh.
#
NUM_PROCS=0
# Create the FFT frame files from the specified time-series.
echo "Starting radiotrans workflow:"
echo "    Generating frame files from ${DATA_PATH}..."
mpirun -np ${NUM_PROCS} python OPT-INSTALL_DIR/waterfall.py "${DATA_PATH}"


# Check that all necessary frame files exist.  If there are some missing, then we'll need to interpolate
# between existent frames to create the missing frames.  NOTE: we run chkwaterfall.py twice in the initial
# run to ensure that as many frame files as possible
# have been created.
echo "    Checking for missing frame files..."
mpirun -np ${NUM_PROCS} python OPT-INSTALL_DIR/chkwaterfall.py "${DATA_PATH}"
FLAG_CONTINUE=1
while [ ${FLAG_CONTINUE} -eq 1 ]
do
   mpirun -np ${NUM_PROCS} python OPT-INSTALL_DIR/chkwaterfall.py "${DATA_PATH}"
   mpirun -np ${NUM_PROCS} python OPT-INSTALL_DIR/eyexam.py "${DATA_PATH}"
   echo "Are there any missing frame files?[y/n]"
   read USER_ANSWER
   # Newer versions of bash have a more elegant way of converting strings to all lower or upper case,
   # however, I'm doing it this way for compatibility.
   USER_ANSWER=$(echo "${USER_ANSWER}" | tr '[:upper:]' '[:lower:]')
   if [[ "${USER_ANSWER}" =~ ${AFFIRMATIVE} ]]; then
      echo "    Interpolating missing frame files..."
      mpirun -np ${NUM_PROCS} python OPT-INSTALL_DIR/interpolate.py "${DATA_PATH}"
   else
      FLAG_CONTINUE=0
   fi
done

# Combine the individual frame files into a single FFT file.
mpirun -np ${NUM_PROCS} python OPT-INSTALL_DIR/waterfallcombine.py "${DATA_PATH}"

# View the spectrograms for the low and high frequency tunings and adjust the bandpass filtration for
# each tuning according to the user's specification.
LOW_FCL=0         # Lower FFT index for the low tuning frequency.
LOW_FCH=4095      # Upper FFT index for the low tuning frequency.
HIGH_FCL=0        # Lower FFT index for the high tuning frequency.
HIGH_FCH=4095     # Upper FFT index for the high tuning frequency.
FLAG_CONTINUE=1
# Before we try to view the bandpass spectrograms, we need to pause the workflow here because the user
# may or may not be online when we get to this point, and the user may or may not have an X11 server and
# X11 forwarding setup.  Pausing here will give the user a chance to get online and perform any
# necessary setup to view the bandpass plots before moving to the next phase.
echo "    Next phase is plotting of bandpass.  This requires X11 and X11 forwarding, if working"
echo "    remotely.  The workflow is currently paused to allow setting up X11 and X11 forwarding."
echo "    Enter any key when ready to proceed:"
read USER_ANSWER
echo "    The next phase involves adjusting the bandpass filtering for RFI.  This can be done"
echo "    interactively with plots (requires X11 and X11 forwarding, if working remotely, be"
echo "    setup), or if the user already knowns the bandpass range to use, it can be entered"
echo "    directly (skipping any display of plots)."
echo "    Adjust bandpass filtering with plots?[y/n] "
read USER_ANSWER
USER_ANSWER=$(echo "${USER_ANSWER}" | tr '[:upper:]' '[:lower:]')

if [[ "${USER_ANSWER}" =~ ${AFFIRMATIVE} ]]; then 
   echo "    Plotting tuning bandpass regions for transient extraction..."
   while [ ${FLAG_CONTINUE} -eq 1 ]
   do
      FLAG_CONTINUE=0
      mpirun -np ${NUM_PROCS} python OPT-INSTALL_DIR/bandpasscheck.py --low-tuning-lower ${LOW_FCL} \
         --low-tuning-upper ${LOW_FCH} --high-tuning-lower ${HIGH_FCL} --high-tuning-upper ${HIGH_FCH}
      echo "    Do the bandpass regions have acceptably low RFI?[y/n]"
      read USER_ANSWER
      USER_ANSWER=$(echo "${USER_ANSWER}" | tr '[:upper:]' '[:lower:]')
      if ! [[ "${USER_ANSWER}" =~ ${AFFIRMATIVE} ]]; then 
         echo "    Change low tuning bandpass?[y/n]"
         read USER_ANSWER
         USER_ANSWER=$(echo "${USER_ANSWER}" | tr '[:upper:]' '[:lower:]')
         if [[ "${USER_ANSWER}" =~ ${AFFIRMATIVE} ]]; then 
            echo "    Enter the lower FFT index for the low tuning: "
            read LOW_FCL
            echo "    Enter the upper FFT index for the low tuning: "
            read LOW_FCH
            FLAG_CONTINUE=1
         fi

         # Change bandpass for high tuning, if desired.
         echo "    Change high tuning bandpass?[y/n]"
         read USER_ANSWER
         USER_ANSWER=$(echo "${USER_ANSWER}" | tr '[:upper:]' '[:lower:]')
         if [[ "${USER_ANSWER}" =~ ${AFFIRMATIVE} ]]; then 
            echo "    Enter the lower FFT index for the high tuning: "
            read HIGH_FCL
            echo "    Enter the upper FFT index for the high tuning: "
            read HIGH_FCH
            FLAG_CONTINUE=1
         fi
      fi
   done
else
   FLAG_CONTINUE=1
   while [ ${FLAG_CONTINUE} -eq 1 ]
   do
      echo "     Current bandpass FFT settings:"
      echo "         Low tuning lower index = ${LOW_FCL}"
      echo "         Low tuning upper index = ${LOW_FCH}"
      echo "         High tuning lower index = ${HIGH_FCL}"
      echo "         High tuning upper index = ${HIGH_FCH}"
      echo "     Do you wish to change any of these?[y/n] "
      read USER_ANSWER
      USER_ANSWER=$(echo "${USER_ANSWER}" | tr '[:upper:]' '[:lower:]')
      if [[ "${USER_ANSWER}" =~ ${AFFIRMATIVE} ]]; then 
         # Change bandpass for low tuning, if desired.
         echo "    Change low tuning bandpass?[y/n]"
         read USER_ANSWER
         USER_ANSWER=$(echo "${USER_ANSWER}" | tr '[:upper:]' '[:lower:]')
         if [[ "${USER_ANSWER}" =~ ${AFFIRMATIVE} ]]; then 
            echo "    Enter the lower FFT index for the low tuning: "
            read LOW_FCL
            echo "    Enter the upper FFT index for the low tuning: "
            read LOW_FCH
         fi

         # Change bandpass for high tuning, if desired.
         echo "    Change high tuning bandpass?[y/n]"
         read USER_ANSWER
         USER_ANSWER=$(echo "${USER_ANSWER}" | tr '[:upper:]' '[:lower:]')
         if [[ "${USER_ANSWER}" =~ ${AFFIRMATIVE} ]]; then 
            echo "    Enter the lower FFT index for the high tuning: "
            read HIGH_FCL
            echo "    Enter the upper FFT index for the high tuning: "
            read HIGH_FCH
         fi
      else
         FLAG_CONTINUE=0
      fi
   done
fi

# Generate new frame files using the bandpass regions specified above.
mpirun -np ${NUM_PROCS} python OPT-INSTALL_DIR/ft.py "${DATA_PATH}" --low-tuning-lower ${LOW_FCL} \
      --low-tuning-upper ${LOW_FCH} --high-tuning-lower ${HIGH_FCL} --high-tuning-upper ${HIGH_FCH}

# Extract information about the time-series for use in doing the de-dispersion and transient extraction
# with dv.py
mpirun -np ${NUM_PROCS} python OPT-INSTALL_DIR/freqtint.py "${DATA_PATH}" --low-tuning-lower ${LOW_FCL} \
      --low-tuning-upper ${LOW_FCH} --high-tuning-lower ${HIGH_FCL} --high-tuning-upper ${HIGH_FCH}

# Perform the de-dispersion for each of the tunings.
mpirun -np ${NUM_PROCS} python OPT-INSTALL_DIR/dv.py "${DATA_PATH}" --lower ${LOW_FCL} \
      --upper ${LOW_FCH} --tuning "low"
mpirun -np ${NUM_PROCS} python OPT-INSTALL_DIR/dv.py "${DATA_PATH}" --lower ${HIGH_FCL} \
      --upper ${HIGH_FCH} --tuning "high"
