# 
# resume.sh
#
# PURPOSE: Implements functions and variables that can be sourced into other scripts to allow command
# and variable resumption
#

# == Script resuming capability ==
#
if [ -z "${RESUME_CMD_FILEPATH}" ]; then
   RESUME_CMD_FILEPATH="./default_cmd.resume"   # Default path to the resume-status file.
fi
if [ -z "${RESUME_VAR_FILEPATH}" ]; then
   RESUME_VAR_FILEPATH="./default_var.resume"   # Default path to the resume-variable file.
fi

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
#     resumecmd [ -l <label> | -s <status> | -f | -k <key> | -h | --help ] <CMD> <CMDARGS>
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
#     -k KEY = trigger full execution of resumecmd based on the value of KEY.  If KEY = 0, then bypass
#              execution of resumecmd and return immediately.  Otherwise, fully execute resumecmd.  This
#              used to deal with commands that are dependent on prior commands executing successfully.
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
 
          resumecmd [ -l <label> | -s <status> | -f | -k <key> | -h | --help ] <CMD> <CMDARGS>
 
  The -l and -s options need to be specified BEFORE the CMD and its subseauent arguments.  Using -h |
  --help first just causes resumecnd to display a help message and then return with status 0.
 
  ENVIRONMENT:
      RESUME_CMD_FILEPATH = path to the file recording the resume status of each command. The default is set
                              to "./default.resume" upon the first execution of resumecmd.  To use a 
                              different file, just set
      RESUME_CMD_FILEPATH to the path of the file before executing resumecmd.

      RESUME_LASTCMD_SUCCESS = flag denoting whether the command executed by resumecmd was successful.
                                 This is used to deal with commands dependent on the successful
                                 execution of prior commands. 
 
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
      -k KEY = trigger full execution of resumecmd based on the value of KEY.  If KEY = 0, then bypass
               execution of resumecmd and return immediately.  Otherwise, fully execute resumecmd.  This
               used to deal with commands that are dependent on prior commands executing successfully.
 
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
            *) # Obtain the command to execute and its arguments; then construct the label to identify it
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
   RESUME_LASTCMD_SUCCESS=1
   if [[ ${RESUME_STATUS} -eq ${RESUME_REPEAT} ]] || [[ ${RESUME_STATUS} -ne ${RESUME_DONE} ]] || \
      [[ ${RESUME_FORCE_EXEC} -eq 1 ]]; then
      # Execute the command and obtain its return status.
      ${RESUME_CMD[*]}
      local CMD_RETURN=${?}
      local NEW_STATUS=

      # Capture the success/fail state of the command.
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



function report_resumecmd()
{
   if [ -n "${RESUME_LASTCMD_SUCCESS}" ]; then
      if [ ${RESUME_LASTCMD_SUCCESS} -eq 1 ]; then
         echo "resumecmd: task SUCCESS"
      else
         echo "resumecmd: task FAIL"
      fi
   else
      echo "resumecmd: task status UNKNOWN"
   fi
}

# === END SCRIPT RESUMING CAPABILITY ==

