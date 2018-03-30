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


# Make sure extended regular expressions are supported.
shopt -s extglob


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




# == Script resuming capability ==
#
RESUME_CMD_FILEPATH="./default_cmd.resume"   # Default path to the resume-status file.
RESUME_VAR_FILEPATH="./default_var.resume"   # Default path to the resume-variable file.

# Command resume status.
#
RESUME_INVALID=-2
RESUME_NOEXIST=-1
RESUME_DONE=0
RESUME_FAIL=1
RESUME_REPEAT=2

# Setup trapping of abort signals from commands.
FLAG_ABORTSIG=0
function __abort_trap () {
   FLAG_ABORTSIG=1
}
trap __abort_trap 6

# Hash command to use for generating resume labels.
HASH_CMD="shasum"
HASH_OPTS=

# Output the resume status associated with the command label ${RESUME_LABEL}.  This should only be used
# by resumecmd.
function __getResumeStatus ()
{
   local RESUME_LABEL="${1}"
   local RESUME_CMD_FILEPATH="${2}"

   local LABEL=
   local STATUS=
   local DUMMY=

   if [ -n "${RESUME_CMD_FILEPATH}" ]; then
      if [ -f "${RESUME_CMD_FILEPATH}" ]; then
         # Validate the resume label.
         if [[ -z "${RESUME_LABEL}" ]] || [[ "${RESUME_LABEL}" =~ *([[:blank:]]|[[:space:]]) ]]; then
            echo "__getResumeStatus: Invalid resume label" >&2
            return 1
         fi

         # Look for the resume line with the specified label and get the associated status, if found.
         local DUMMY=
         while read -r LABEL DUMMY STATUS
         do
            if [[ "${LABEL}  -" == "${RESUME_LABEL}" ]]; then
               break;
            fi
         done < "${RESUME_CMD_FILEPATH}"

         if [[ "${LABEL}  -" == "${RESUME_LABEL}" ]]; then
            case ${STATUS} in
               ${RESUME_DONE}|${RESUME_FAIL}|${RESUME_REPEAT}) # Output the valid resume status.
                  echo ${STATUS}
                  ;;
               *) # Anything else is an invalid resume status.
                  echo ${RESUME_INVALID}
                  ;;
            esac
         else
            echo ${RESUME_NOEXIST}
         fi

      else
         echo "__getResumeStatus: Cannot find resume file ${RESUME_CMD_FILEPATH}" >&2
         return 1
      fi
   else
      echo "__getResumeStatus: No resume file path specified" >&2
      return 1
   fi
   return 0
}

# function resumecmd ()
#
# PURPOSE: Determines whether a given command with arguments should be executed based on its associated
# resume status that is saved in a resume-status file.  The resume-status file is specified by setting
# the global variable RESUME_CMD_FILEPATH, which is defaulted to "./default.resume"
# script. The usage is the following:
#
#     resumecmd [ -l <label> | -s <status> | -f | -h | --help ] <CMD> <CMDARGS>
#
# The -l and -s options need to be specified BEFORE the CMD and its subseauent arguments.  Using -h |
# --help first just causes resumecnd to display a help message and then return with status 0.
#
# ENVIRONMENT:
#     RESUME_CMD_FILEPATH = path to the file recording the resume status of each command. The default is set
#     to "./default.resume" upon the first execution of resumecmd.  To use a different file, just set
#     RESUME_CMD_FILEPATH to the path of the file before executing resumecmd.
#
# INPUTS:
#     CMD = the command to execute from the main script.
#     CMDARGS = arguments to the command to execute.
#
# OPTIONALS:
#     -l LABEL = a specific label attached to the hash of the command to uniquely identify it when
#                 looking up its resume status.
#     -s STATUS = a forced resume status on the command. This is usually made to be RESUME_REPEAT to
#                 force the command to be repeated on subsequent runs of the script.
#     -f = force execution of the command, independent of the current resume status.  The resume status
#           is updated after execution.
#
# OUTPUTS:
#    file with resume status of commands executed with this function. 
#
function resumecmd ()
{
   local RESUME_USAGE='
  NAME:
      resumecmd -- resume command

  PURPOSE: 
      Determines whether a given command with arguments should be executed based on its associated
      resume status that is saved in a resume-status file.  The resume-status file is specified by setting
      the global variable RESUME_CMD_FILEPATH, which is defaulted to "./default.resume"
      script. The usage is the following:
 
          resumecmd [ -l <label> | -s <status> | -f | -h | --help ] <CMD> <CMDARGS>
 
  The -l and -s options need to be specified BEFORE the CMD and its subseauent arguments.  Using -h |
  --help first just causes resumecnd to display a help message and then return with status 0.
 
  ENVIRONMENT:
      RESUME_CMD_FILEPATH = path to the file recording the resume status of each command. The default is set
      to "./default.resume" upon the first execution of resumecmd.  To use a different file, just set
      RESUME_CMD_FILEPATH to the path of the file before executing resumecmd.

      RESUME_LASTCMD_SUCCESS = flag denoting whether the command executed by resumecmd was successful.
 
  INPUTS:
      CMD = the command to execute from the main script.
      CMDARGS = arguments to the command to execute.
 
  OPTIONALS:
      -l LABEL = a specific label attached to the hash of the command to uniquely identify it when
                  looking up its resume status.
      -s STATUS = a forced resume status on the command. This is usually made to be RESUME_REPEAT to
                  force the command to be repeated on subsequent runs of the script.
      -f = force execution of the command, independent of the current resume status.  The resume status
            is updated after execution.
 
  OUTPUTS:
     file with resume status of commands executed with this function. 
'

   local RESUME_CMD=             # Command with arguments to execute.
   local RESUME_LABEL=           # Label to identify the command in the resume-status file.
   local RESUME_STATUS=          # Resume status of the command.
   local RESUME_FORCE_STATUS=    # Resume status to force the command to have after execution.
   local RESUME_FORCE_EXEC=0     # Flag denoting whether to force execute the command, independent of
                                 # its current resume status.

   local RESUME_KEY=1            # Triggering key denoting whether to bypass the execution of resumecmd.
                                 # If the value is 1, then resumecmd executes as normal.  If the value
                                 # is 0, then resumecmd returns immediately with return value 0.  This
                                 # is used for dealing with commands that are dependent on the
                                 # successful execution of prior commands.


   # Set the default resume-status filepath, if one has not been specified.
   if [ -z "${RESUME_CMD_FILEPATH}" ]; then
      RESUME_CMD_FILEPATH="./default_cmd.resume"
   fi

   # If resumecmd is given arguments, parse them; otherwise, just return from resumecmd.
   if [[ ${#} -gt 0 ]]; then
      while [[ -z ${RESUME_CMD} ]] && [[ -n "${1}" ]]
      do
         case "${1}" in
            -l) # Obtain the unique label to attach to the resume status.
               RESUME_LABEL="${2}"
               shift; shift
               ;;
            -s) # Force the resume status to the given status.
               RESUME_FORCE_STATUS="${2}"
               shift; shift
               ;;
            -f) # Force repeat the command, independent of any existing resume status.
               RESUME_FORCE_EXEC=1
               shift
               ;;
            -k) # Set the triggering key value.  A value of 0 will bypass the execution of resumecmd and
                # cause it to immediately return with return value 0.  Any other value will cause
                # resumecmd to execute as normal.
                RESUME_KEY=${2}
                shift; shift
                ;;
            -h | --help) # Display resumecmd help and then return.
               echo "${RESUME_USAGE}"
               return 0
               ;;
            *) # Obtain the command to execute and its arguments, and construct the label to identify it
               # in the resume-status file.
               RESUME_CMD=("${@}")
               RESUME_LABEL="${RESUME_LABEL}-$(echo "${RESUME_CMD[*]}" | ${HASH_CMD} ${HASH_OPTS[*]})"
               ;;
         esac
      done
   else
      return 0
   fi

   # Set {RESUME_LASTCMD_SUCCESS}, if it has not been set.
   if [ -z "${RESUME_LASTCMD_SUCCESS}" ]; then
      RESUME_LASTCMD_SUCCESS=0
   fi

   # Check the resume key for whether to execute the rest of resumecmd.
   if [[ ${RESUME_KEY} -eq 0 ]]; then
      return 0
   fi

   # If the resume file does not exist, create it.
   if [ ! -f "${RESUME_CMD_FILEPATH}" ]; then
      touch "${RESUME_CMD_FILEPATH}"
      # If the resume-status file could not created, then notify the user and execute the command
      # without saving the resume state.
      if [ ! -f "${RESUME_CMD_FILEPATH}" ]; then
         echo "resumecmd: Could not find or create resume file ${RESUME_CMD_FILEPATH}" >&2
         echo "resumecmd: Command executing without saved resume status => ${RESUME_CMD[*]}" >&2
         ${RESUME_CMD[*]}
         # Capture the success/fail state of the command.
         local CMD_RETURN=${?}
         RESUME_LASTCMD_SUCCESS=1
         if [[ ${CMD_RETURN} -ne 0 ]] || [[ ${FLAG_ABORTSIG} -eq 1 ]]; then
            RESUME_LASTCMD_SUCCESS=0
         fi
         FLAG_ABORTSIG=0
         return 1
      fi
   fi

   # Check the resume file for the current resume status of the command.
   RESUME_STATUS=$(__getResumeStatus "${RESUME_LABEL}" "${RESUME_CMD_FILEPATH}")

   # Execute the command and record the resulting resume status if its current resume status is
   # RESUME_REPAT or not RESUME_DONE, or if execution has been forced.
   if [[ ${RESUME_STATUS} -eq ${RESUME_REPEAT} ]] || [[ ${RESUME_STATUS} -ne ${RESUME_DONE} ]] || \
      [[ ${RESUME_FORCE_EXEC} -eq 1 ]]; then
      # Execute the command and obtain its return status.
      ${RESUME_CMD[*]}
      local CMD_RETURN=${?}
      local NEW_STATUS=

      # Capture the success/fail state of the command.
      RESUME_LASTCMD_SUCCESS=1
      if [[ ${CMD_RETURN} -ne 0 ]] || [[ ${FLAG_ABORTSIG} -eq 1 ]]; then
         RESUME_LASTCMD_SUCCESS=0
      fi
     
      # Check if resume status is forced and record the forced status as the new status. Otherwise, 
      # determine the appropriate status to record.
      if [ -n "${RESUME_FORCE_STATUS}" ]; then
         NEW_STATUS=${RESUME_FORCE_STATUS}
      else
         # If this is not a repeating command, the determine whether it succeeded or failed.
         NEW_STATUS=${RESUME_REPEAT}
         if [[ ${RESUME_STATUS} -ne ${RESUME_REPEAT} ]]; then
            if [[ ${CMD_RETURN} -ne 0 ]] || [[ ${FLAG_ABORTSIG} -eq 1 ]]; then
               NEW_STATUS=${RESUME_FAIL}
            else
               NEW_STATUS=${RESUME_DONE}
            fi
         fi
      fi

      # Record the resume status: add it if it doesn't exist; otherwise, update it.
      if [[ ${RESUME_STATUS} -eq ${RESUME_NOEXIST} ]]; then
         echo "${RESUME_LABEL}" ${NEW_STATUS} >> "${RESUME_CMD_FILEPATH}"
      else
         sed -i -e s@"${RESUME_LABEL}.*"@"${RESUME_LABEL} ${NEW_STATUS}"@ "${RESUME_CMD_FILEPATH}"
      fi
      FLAG_ABORTSIG=0
   fi
}

# Output the resume variable value associated with the variable label ${RESUME_LABEL}.  This should only 
# be used by resumevar
function __getResumeValue ()
{
   local RESUME_LABEL="${1}"
   local RESUME_FILEPATH="${2}"

   local LABEL=
   local VARVAUE=
   local DUMMY=

   if [ -n "${RESUME_FILEPATH}" ]; then
      if [ -f "${RESUME_FILEPATH}" ]; then
         # Validate the resume label.
         if [[ -z "${RESUME_LABEL}" ]] || [[ "${RESUME_LABEL}" =~ *([[:blank:]]|[[:space:]]) ]]; then
            echo "__getResumeValue: Invalid resume label" >&2
            return 2
         fi

         # Look for the resume line with the specified label and get the associated variable name and
         # its value, if found.
         while read -r LABEL DUMMY VARVALUE
         do
            if [[ "${LABEL}  -" == "${RESUME_LABEL}" ]]; then
               break;
            fi
         done < "${RESUME_FILEPATH}"

         if [[ "${LABEL}  -" == "${RESUME_LABEL}" ]]; then
            echo "${VARVALUE[*]}"
            return 0
         else
            echo
            return 1
         fi

      else
         echo "__getResumeValue: Cannot find resume file ${RESUME_FILEPATH}" >&2
         return 2
      fi
   else
      echo "__getResumeValue: No resume file path specified" >&2
      return 2
   fi
}

# function resumevar ()
#
# PURPOSE: Implements resuming the value of a specified variable from the resume-variable file.  If the
# value is not found in the resume-variable file, then the variable is assigned the specified value;
# that value is then saved to the resume-variable file.  Otherwise, the variable's value is
# reloaded from the resume-variable file regardless of the value specified.  The usage is the following:
#
#     resumevar [ -l <LABEL> | -f | -h | --help ] <VARNAME> <VARVALUE[*]>
#
# The -f option forces the variable to be assigned the specified value and update the resume-variable
# file with that value.  The -h and --help options display a help message for resumevar and then
# immediately returns.  The -l <LABEL> option adds a unique identifying label to distinguish different
# instances of the same name of the variable.
#
# INPUTS:
#     <VARNAME> = name of the variable.
#     <VARVALUE[*]> = value to assign to the variable.
#  
# OPTIONALS:
#     -l <LABEL> : add unique tag <LABEL> to distinguish one instance of a variable name from another.
#     -f : flag to force assignment and update of the variable's value.
#     -h | --help : display this help message and return.
#
function resumevar ()
{
   local VARNAME=
   local VARVALUE=
   local VARLABEL=
   local FORCE_VALUE=0
   local VAR_NOTFOUND=0

   # Set the default resume-variable filepath, if one has not been specified.
   if [ -z "${RESUME_VAR_FILEPATH}" ]; then
      RESUME_VAR_FILEPATH="./default_var.resume"
   fi

   if [[ ${#} -gt 0 ]]; then
      # Parse arguments.
      while [ -z "${VARNAME}" ] && [[ -n ${1} ]]
      do
         case ${1} in 
            -l) # Set the unique identifying label.
               VARLABEL="${2}"
               shift; shift
               ;;
            -f) # Set forcing the variable's value.
               FORCE_VALUE=1
               shift
               ;;
            -h | --help) # Display the help message and then return.
               ;;
            *) # Obtain the variable name and its value.  Also, create the label to identify the
               # variable in the resume-variable file.
               VARNAME="${1}"
               shift
               VARVALUE=("${@}")
               VARLABEL="${VARLABEL}-$(echo "${VARNAME}" | ${HASH_CMD} ${HASH_OPTS[*]})"
               ;;
         esac
      done

      # If the resume-variable filepath doesn't exist, then create it.
      if [ ! -f "${RESUME_VAR_FILEPATH}" ]; then
         touch "${RESUME_VAR_FILEPATH}"
         # If the resume-variable file can not be created, notify the user and assign the variable
         # without saving it.
         if [ ! -f "${RESUME_VAR_FILEPATH}" ]; then
            echo "resumevar: Could not create resume-variable file ${RESUME_VAR_FILEPATH}" >&2
            echo "resumevar: Assigning variable without resume => ${VARNAME}" >&2
            export ${VARNAME}="${VARVALUE[*]}"
            return 1
         fi
      fi

      # Find the variable and its value in the resume-variable file.
      RESUME_VALUE=$(__getResumeValue "${VARLABEL}" "${RESUME_VAR_FILEPATH}")
      VAR_NOTFOUND=${?}

      # If the variable's value is not found or is being forced, then assign the value that has been
      # specified and save it to the resume-variable file.  Otherwise, use the value loaded from the
      # resume-variable file.
      if [[ ${VAR_NOTFOUND} -ne 0 ]] || [[ ${FORCE_VALUE} -eq 1 ]]; then
         export ${VARNAME}="${VARVALUE[*]}"

         # Record the variable's value: add it if it doesn't exist; otherwise, update it.
         if [[ ${VAR_NOTFOUND} -ne 0 ]]; then
            echo "${VARLABEL} ${VARVALUE[*]}" >> "${RESUME_VAR_FILEPATH}"
         else
            sed -i -e s@"${VARLABEL}.*"@"${VARLABEL} ${VARVALUE[*]}"@ "${RESUME_VAR_FILEPATH}"
         fi
      else
         export ${VARNAME}="${RESUME_VALUE}"
      fi

   fi

   return 0
}


# === END SCRIPT RESUMING CAPABILITY ==





# ==== MAIN WORKFLOW FOR RADIOTRANS.SH ===
#
#
DATA_PATH=        # Path to the radio time-series data file.
RESUME_CMD_FILEPATH="./radiotrans_cmd.resume"
RESUME_VAR_FILEPATH="./radiotrans_var.resume"


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




# == IMPLEMENT RADIOTRANS WORKFLOW ==
#
# Workflow resume labels.
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
NUM_PROCS=0


# Create the FFT frame files from the specified time-series.
echo "Starting radiotrans workflow:"
echo "    Generating frame files from ${DATA_PATH}..."
resumecmd -l ${LBL_WATERFALL} mpirun -np ${NUM_PROCS} python OPT-INSTALL_DIR/waterfall.py "${DATA_PATH}"

exit 0

# Check that all necessary frame files exist.  If there are some missing, then we'll need to interpolate
# between existent frames to create the missing frames.  NOTE: we run chkwaterfall.py twice in the initial
# run to ensure that as many frame files as possible
# have been created.
echo "    Checking for missing frame files..."
resumecmd -l ${LBL_CHKWATERFALLA} -k ${RESUME_LASTCMD_SUCCESS} \
   mpirun -np ${NUM_PROCS} python OPT-INSTALL_DIR/chkwaterfall.py "${DATA_PATH}"
FLAG_CONTINUE=1
RESUME_FORCE_OPT=
while [[ ${FLAG_CONTINUE} -eq 1 ]] && [[ ${RESUME_LASTCMD_SUCCESS} -ne 0 ]]
do
   resumecmd -l ${LBL_CHKWATERFALLB} ${RESUME_FORCE_OPT} -k ${RESUME_LASTCMD_SUCCESS} \
      mpirun -np ${NUM_PROCS} python OPT-INSTALL_DIR/chkwaterfall.py "${DATA_PATH}"
   resumecmd -l ${LBL_EYEEXAM} ${RESUME_FORCE_OPT} -k ${RESUME_LASTCMD_SUCCESS} \
      mpirun -np ${NUM_PROCS} python OPT-INSTALL_DIR/eyexam.py "${DATA_PATH}"
   resumecmd ${RESUME_FORCE_OPT} -k ${RESUME_LASTCMD_SUCCESS} \
      echo "Are there any missing frame files?[y/n]"
   resumecmd ${RESUME_FORCE_OPT} -k ${RESUME_LASTCMD_SUCCESS} read USER_ANSWER
   # Newer versions of bash have a more elegant way of converting strings to all lower or upper case,
   # however, I'm doing it this way for compatibility.
   resumecmd ${RESUME_FORCE_OPT} -k ${RESUME_LASTCMD_SUCCESS} \
      USER_ANSWER=$(echo "${USER_ANSWER}" | tr '[:upper:]' '[:lower:]')
   if [[ "${USER_ANSWER}" =~ ${AFFIRMATIVE} ]]; then
      echo "    Interpolating missing frame files..."
      resumecmd -l ${LBL_INTERPOLATE} ${RESUME_FORCE_OPT} -k ${RESUME_LASTCMD_SUCCESS} \
         mpirun -np ${NUM_PROCS} python OPT-INSTALL_DIR/interpolate.py "${DATA_PATH}"
      # Force resumecmd to execute commands in the loop because we need to do all this again.
      RESUME_FORCE_OPT="-f"
   else
      FLAG_CONTINUE=0
   fi
done


# Combine the individual frame files into a single FFT file.
resumecmd -l ${LBL_COMBINE} -k ${RESUME_LASTCMD_SUCCESS} \
   mpirun -np ${NUM_PROCS} python OPT-INSTALL_DIR/waterfallcombine.py "${DATA_PATH}"

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
   RESUME_FORCE_OPT=
   if [[ "${USER_ANSWER}" =~ ${AFFIRMATIVE} ]]; then 
      echo "    Plotting tuning bandpass regions for transient extraction..."
      while [[ ${FLAG_CONTINUE} -eq 1 ]]
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
               read USER_ANSWER
               resumevar -l ${LBL_LOWFCL} -f LOW_FCL ${USER_ANSWER}
               echo "    Enter the upper FFT index for the low tuning: "
               read USER_ANSWER
               resumevar -l ${LBL_LOWFCH} -f LOW_FCH ${USER_ANSWER}
               FLAG_CONTINUE=1
               RESUME_FORCE_OPT="-f"
            fi

            # Change bandpass for high tuning, if desired.
            echo "    Change high tuning bandpass?[y/n]"
            read USER_ANSWER
            USER_ANSWER=$(echo "${USER_ANSWER}" | tr '[:upper:]' '[:lower:]')
            if [[ "${USER_ANSWER}" =~ ${AFFIRMATIVE} ]]; then 
               echo "    Enter the lower FFT index for the high tuning: "
               read USER_ANSWER
               resumevar -l ${LBL_HIGHFCL} -f HIGH_FCL ${USER_ANSWER}
               echo "    Enter the upper FFT index for the high tuning: "
               read USER_ANSWER
               resumevar -l ${LBL_HIGHFCH} -f HIGH_FCH ${USER_ANSWER}
               FLAG_CONTINUE=1
               RESUME_FORCE_OPT="-f"
            fi
         fi
      done
   else
      FLAG_CONTINUE=1
      while [[ ${FLAG_CONTINUE} -eq 1 ]]
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
               read USER_ANSWER
               resumevar -l ${LBL_LOWFCL} -f LOW_FCL ${USER_ANSWER}
               echo "    Enter the upper FFT index for the low tuning: "
               read USER_ANSWER
               resumevar -l ${LBL_LOWFCH} -f LOW_FCH ${USER_ANSWER}
               RESUME_FORCE_OPT="-f"
            fi

            # Change bandpass for high tuning, if desired.
            echo "    Change high tuning bandpass?[y/n]"
            read USER_ANSWER
            USER_ANSWER=$(echo "${USER_ANSWER}" | tr '[:upper:]' '[:lower:]')
            if [[ "${USER_ANSWER}" =~ ${AFFIRMATIVE} ]]; then 
               echo "    Enter the lower FFT index for the high tuning: "
               read USER_ANSWER
               resumevar -l ${LBL_HIGHFCL} -f HIGH_FCL ${USER_ANSWER}
               echo "    Enter the upper FFT index for the high tuning: "
               read USER_ANSWER
               resumevar -l ${LBL_HIGHFCH} -f HIGH_FCH ${USER_ANSWER}
               RESUME_FORCE_OPT="-f"
            fi
         else
            FLAG_CONTINUE=0
         fi
      done
   fi
fi

# Generate new frame files using the bandpass regions specified above.
echo "     Generating new frame files using specified bandpass regions..."
resumecmd -l ${LBL_FRAMETRIM} ${RESUME_FORCE_OPT} -k ${RESUME_LASTCMD_SUCCESS} \
   mpirun -np ${NUM_PROCS} python OPT-INSTALL_DIR/ft.py "${DATA_PATH}" \
   --low-tuning-lower ${LOW_FCL} \
   --low-tuning-upper ${LOW_FCH} --high-tuning-lower ${HIGH_FCL} --high-tuning-upper ${HIGH_FCH}

# Extract information about the time-series for use in doing the de-dispersion and transient extraction
# with dv.py
echo "     Generating time-series information for de-dispersion..."
resumecmd -l ${LBL_FREQTINT} ${RESUME_FORCE_OPT} -k ${RESUME_LASTCMD_SUCCESS} \
   mpirun -np ${NUM_PROCS} python OPT-INSTALL_DIR/freqtint.py "${DATA_PATH}" \
   --low-tuning-lower ${LOW_FCL} \
   --low-tuning-upper ${LOW_FCH} --high-tuning-lower ${HIGH_FCL} --high-tuning-upper ${HIGH_FCH}

# Perform the de-dispersion for each of the tunings.
echo "     Performing de-dispersion on low tuning data..."
resumecmd -l ${LBL_DEDISPERSIONA} ${RESUME_FORCE_OPT} -k ${RESUME_LASTCMD_SUCCESS} \
   mpirun -np ${NUM_PROCS} python OPT-INSTALL_DIR/dv.py "${DATA_PATH}" \
   --lower ${LOW_FCL} --upper ${LOW_FCH} --tuning "low"
echo "     Performing de-dispersion on high tuning data..."
resumecmd -l ${LBL_DEDISPERSIONB} ${RESUME_FORCE_OPT} -k ${RESUME_LASTCMD_SUCCESS} \
   mpirun -np ${NUM_PROCS} python OPT-INSTALL_DIR/dv.py "${DATA_PATH}" \
   --lower ${HIGH_FCL} --upper ${HIGH_FCH} --tuning "high"


# Handle cleanup, if everything went okay.  We will no longer need all those frame files and other
# temporary work files hanging around.
if [[ ${RESUME_LASTCMD_SUCCESS} -ne 0 ]]; then
   echo "Performing cleanup of temporary files..."
   rm -f ./waterfall.npy
   rm -f ./waterfall*.npy
   rm -f ./tInt.npy
   rm -f ./freq1.npy
   rm -f ./freq2.npy
   rm -f ./*.npy
fi
echo "Workflow done!"
exit 0

