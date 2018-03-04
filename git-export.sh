#!/bin/bash

# Script to export the current version of the named branch or tag in a git repository to the
# specified install directory as an installed application package.

USAGE='
git-export.sh [--install-dir | -d <root>] [--repos | -r <repos>] [--git-tag | -t <tag>] 
                       [--no-prefix | -n] <prefix> 

Exports a specified branch or tag <tag> from the git repository <repos> as an install under the
directory <root>/<prefix>/. If <tag> is not specified, then <tag> defaults to "master".  If <repos> is
not specified, then it is assumed that git-export.sh is being run from within the repository to export.
If <root> is not specified, then <root> defaults to "$HOME/local/".  If <prefix> is not specified,
then <prefix> defaults to the current directory name (this assumes that git-export.sh is being run from
the repository directory, and it is appropriately name).

   --install-root | -r <root> : set the root install directory to <root>.  Default = "$HOME/local/"

   --repos | -p <repos> : path to the repository to export.  Default = "./"

   --git-tag | -t <tag> : branch or tag to export from repository.  Default = "master"

   --no-prefix | -n : do not use a prefix when exporting.  This will cause the contents of the
                      repository to be dumped directly into "<root>/" instead of "<root>/<prefix>/"

   --help | -h : display this help.

NOTE: This export removes any special .git* files that are committed in the repository.  If it is
desired to keep such files, then the user should perform the export directly with "git archive".
'

INSTALL_ROOT=
GIT_TAG=
GIT_DIR=
PREFIX=
FLAG_NOPREFIX=0


# Parse command-line
while [ -n "$1" ]
do 
   case "$1" in
      --install-root | -r) # Set the install root directory.
         if [ -z "${INSTALL_ROOT}" ]; then
            INSTALL_ROOT="$2"
         fi
         shift; shift
         ;;
      --repos | -p) # Specify the repository directory.
         if [ -z "${GIT_DIR}" ]; then
            GIT_DIR="$2"
         fi
         shift; shift
         ;;
      --git-tag | -t) # Specify the git tag or branch to export.
         if [ -z "${GIT_TAG}" ]; then
            GIT_TAG="$2"
         fi
         shift; shift
         ;;
      --no-prefix | -n) # Do not use a package prefix to the export.
         FLAG_NOPREFIX=1
         shift
         ;;
      --help | -h) # Display usage.
         echo "${USAGE}"
         exit 0
         ;;
      -*) # Unknown option. 
         echo "UNKNOWN OPTION: $1"
         echo "${USAGE}"
         exit 0
         ;;
      *) # The first non-option argument encountered is considered the install package prefix.
         if [ -z "${PREFIX}" ]; then
            PREFIX="$1"
         fi
         shift
         ;;
   esac
done

# Determine defaults for values not set.
if [ -z "${INSTALL_ROOT}" ]; then
   INSTALL_ROOT="${HOME}/local/"
fi
if [ -z "${GIT_TAG}" ]; then
   GIT_TAG="master"
fi
if [ -z "${PREFIX}" ]; then
   CURR_DIR=`pwd`
   PREFIX="${CURR_DIR##*/}"
fi

# Remove the package prefix, if specified.
if [ ${FLAG_NOPREFIX} -eq 1 ]; then
   PREFIX="."
fi
#
#


# Go into the git repository directory, if one specified.
if [ -n "${GIT_DIR}" ]; then
   if [ -d "${GIT_DIR}" ]; then
      pushd "${GIT_DIR}"
   else
      echo "${GIT_DIR} not found or is not a directory"
      exit 1
   fi
fi

# Check that the current directory is a git repository.
EXIT_CODE=0
if [ -d ".git" ]; then
   # create a temporary archive of the git repository for export, then deploy that archive to the install
   # directory.
   echo "Creating git archive..."
   ARCHIVE_FILE="temp_export.tar"
   git archive ${GIT_TAG} --prefix="${PREFIX}/" --format=tar -o ${ARCHIVE_FILE}
   if [ ! -d "${INSTALL_ROOT}" ]; then
      mkdir "${INSTALL_ROOT}"
   fi

   # Extract the archive into the install directory and delete the temporary archive.
   echo "Extracting git archive to ${INSTALL_ROOT}/..."
   tar -xvf "${ARCHIVE_FILE}" -C "${INSTALL_ROOT}"
   rm -f "${ARCHIVE_FILE}"

   # In any *.sh shell scripts, change '/home/cyancey/bin' to install path.
   INSTALL_DIR="${INSTALL_ROOT}"
   if [ "${PREFIX}" != "." ]; then
      INSTALL_DIR="${INSTALL_ROOT}/${PREFIX}"
   fi
   echo "Replacing OPT-INSTALL_DIR template in all *.sh scripts with ${INSTALL_DIR}..."
   for  FILEPATH in `ls ${INSTALL_DIR}/*.sh`
   do
      echo "    Replacing in ${FILEPATH}."
      sed -i s@OPT-INSTALL_DIR@${INSTALL_DIR}@ "${FILEPATH}"
   done

   # Remove any special .git* files, like .gitignore, from the installed export.
   # Also, remove git-export.sh from the install directory, if it was included with the repository as a
   # utility dependency.  This is to avoid collisions on ${PATH} and other accidents.
   echo "Removing .git* files and git-export.sh inclusion from installed export..."
   for FILENAME in `ls ${INSTALL_DIR}/.git*`
   do
      echo "    Removing ${FILENAME}."
      rm -f "${FILENAME}"
   done
   GITEXPORT_SH=`find ${INSTALL_DIR} -name git-export.sh`
   if [ -f "${GITEXPORT_SH}" ]; then
      echo "    Removing ${GITEXPORT_SH}."
      echo "    If you need git-export.sh, copy it manually to an appropriate directory on your PATH."
      rm -f "${GITEXPORT_SH}"
   fi

else
   echo "Not in a git repository."
   EXIT_CODE=1
fi


if [ -n "${GIT_DIR}" ]; then
   popd
fi
exit ${EXIT_CODE}

