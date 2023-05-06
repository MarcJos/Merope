# At the CEA
PYBIND_REPO=https://www-git-cad.intra.cea.fr/DEC/collaboratif/mj263790/copy_pybind
GCC10_ENV=/soft/pleiades/testing/compilers/gcc/gcc-12.1/buster/x86_64/env.sh

pushd . > /dev/null
SCRIPT_PATH="${BASH_SOURCE[0]}"
if ([ -h "${SCRIPT_PATH}" ]); then
  while([ -h "${SCRIPT_PATH}" ]); do cd `dirname "$SCRIPT_PATH"`;
  SCRIPT_PATH=`readlink "${SCRIPT_PATH}"`; done
fi
cd `dirname ${SCRIPT_PATH}` > /dev/null
SCRIPT_PATH=`pwd`;
popd  > /dev/null
export PYTHONPATH=${SCRIPT_PATH}/BUILD-DIR/src:${SCRIPT_PATH}/BUILD-DIR/Interface_python:$PYTHONPATH

source $GCC10_ENV
