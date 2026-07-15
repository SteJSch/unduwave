#!/bin/bash

do_flag=''
b_flag=''
files=''
verbose='false'
buildWAVE='false'
buildWAVEW='false'
buildWAVEL='false'
buildUndu='false'
buildUNDUWAVE='false'

print_usage() {
  printf 'Usage: %s -d buw|bu|bw|bpy|bpya\n' "$0"
  printf '  -d  buw - build undumag and wave, bu - build undumag, bw - build wave,\n'
  printf '      bpy - build the python project, bpya - build all external and python\n'
  exit 1
}

if [[ "${1:-}" == "--help" ]]; then
  print_usage
  exit 0
fi

while getopts 'd:' flag; do
  case "${flag}" in
    d) do_flag="${OPTARG}" ;;
    *) print_usage; exit 1 ;;  esac
done

if [ "$do_flag" == uw_appt ]; then
	echo "Downloading Undumag and Wave into apptainer"
fi
if [ "$do_flag" == buw ]; then
	buildWAVE='true'
	buildUndu='true'
fi
if [ "$do_flag" == bu ]; then
	buildUndu='true'
fi
if [ "$do_flag" == bw ]; then
	buildWAVE='true'
	buildWAVEL='true'
	buildWAVEW='true'
fi	
if [ "$do_flag" == bwl ]; then
	buildWAVEL='true'
fi	
if [ "$do_flag" == bww ]; then
	buildWAVEW='true'
fi	
if [ "$do_flag" == bpy ]; then
	buildUNDUWAVE='true'
fi
if [ "$do_flag" == bpya ]; then
	buildUNDUWAVE='true'
	buildWAVE='true'
	buildUndu='true'
fi

if [ "$buildUndu" == 'true' ]; then
	echo "Building Undumag For Windows"
	pwd=$PWD
	direc_undu=$pwd'/../unduwave/External-Software/UNDUMAG'
	cd $direc_undu
	cd lib
	rm *
	cd ..
	cp $pwd/compile_undumag_incl_win.sh $direc_undu/shell/
	cp $pwd/make_undumag_win2.py $direc_undu/python/
	python3 python/make_undumag_win2.py
	echo "Building Undumag For Linux"
	cd lib
	rm *
	cd ..
	python3 python/make_undumag.py	
fi
if [ "$buildWAVEL" == 'true' ]; then
	echo "Building WAVE For Linux"
	pwd=$PWD
	direc_wave=$pwd'/../unduwave/External-Software/WAVE'
	cd $direc_wave
	cd lib
	rm *
	cd ..
	python3 python/make_wave.py	
fi
if [ "$buildWAVEW" == 'true' ]; then
	echo "Building WAVE For Windows"
	pwd=$PWD
	direc_wave=$pwd'/../unduwave/External-Software/WAVE'
	cd $direc_wave
	cd lib
	rm *
	cd ..
	cp $pwd/make_wave_win2.py $direc_wave/python/
	python3 python/make_wave_win2.py
fi
if [ "$buildUNDUWAVE" == 'true' ]; then
	pwd=$PWD
	direc_unduwave=$pwd'/../'
	cd $direc_unduwave
	python3 -m build
	python3 -m pip install --editable .
fi






















