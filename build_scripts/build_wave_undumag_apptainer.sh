#!/bin/bash

do_flag=''
b_flag=''
files=''
verbose='false'

print_usage() {
  printf "Usage: ..."
}

while getopts 'd:' flag; do
  case "${flag}" in
    d) do_flag="${OPTARG}" ;;
    *) print_usage
       exit 1 ;;
  esac
done

if [ "$do_flag" == uw_appt ]; then
	echo "Downloading Undumag and Wave into apptainer"
else
	echo "Downloading Other Things"
fi
