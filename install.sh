#!/bin/bash

function project_build ()
{
  echo "Type=${1}, Tests=${2}, Static=${3}, Dev=${4}, j=${5}"

  mkdir -p build
  cd build

  cmake .. -DCMAKE_BUILD_TYPE=${1} \
           -DWITH_TESTS=${2} \
           -DWITH_COVERAGE=${2} \
           -DSTATIC_BUILD=${3} \
           -DDEV_BUILD="${4}"
  make -j${5}
  alrighty=$(($? == 0))  # is the return code of make 0 ?
  if [[ ${1} == "Coverage" && ${alrighty} == 1 ]]; then
    make coverage-report
  elif [[ ${6} == 1 && ${alrighty} == 1 ]]; then
    ctest --verbose
  fi
}

function usage ()
{
  echo "fimpera build script."
  echo "Usage: "
  echo "  ./install.sh [-r str] [-t int] [-j int] [-d] [-s]"
  echo "Options: "
  echo "  -r <Release|Debug> -> build type {Release}."
  echo "  -t <0|1|2>         -> tests: 0 = disabled, 1 = compile, 2 = compile and run {2}."
  echo "  -j <INT>           -> nb threads {8}."
  echo "  -d                 -> dev build {disabled}."
  echo "  -s                 -> static build {disabled}."
  echo "  -h                 -> show help."
  exit 1
}

rel="Release"
deb="Debug"
cov="Coverage"

mode="Release"
dev="OFF"
static="OFF"
tests=2
tests_str="ON"
tests_run=1
jopt=8


while getopts "r:k:t:c:j:onwmsh" option; do
  case "$option" in
    r)
      mode=${OPTARG}
      [[ ${mode} != ${rel} && ${mode} != ${deb} && ${mode} != ${cov} ]] && usage
      ;;
    t)
      tests=${OPTARG}
      [[ ${tests} == 0 ]] || [[ ${tests} == 1 ]] || [[ ${tests} == 2 ]] || usage
      [[ ${tests} == 0 ]] && tests_str="OFF"
      [[ ${tests} == 1 ]] && tests_run=0
      ;;
    j)
      jopt=${OPTARG}
      ;;
    s)
      static="ON"
      ;;
    d)
      dev="ON"
      ;;
    h)
      usage
      ;;
    *)
      usage
      ;;
  esac
done

project_build ${mode} ${tests_str} ${static} ${dev} ${jopt} ${tests_run}