#!/bin/bash
#
# Filename: get_data.sh
# Author: Artem Baranovskii
# Date: 15th November 2023
# Description: Script to download SEACell-processed scRNA-seq data to use in the development
#
# Usage: ./get_data.sh [options] arguments
#
# Options:
#   -h, --help        Show brief help message
#   -d, --datasets    Specify the names of the datasets to download; Can be any of: Jansky, JanskyAdrenal, Kildisiute, Dong, Bedoya, and all"
#   -o, --output DIR  Specify a target directory to download the data
#
# Example:
#   ./get_data.sh -d Jansky Dong -o /path/to/output
#
set -e
while (( "$#" )); do
    case "$1" in
        -h|--help)
            echo "options:"
            echo "-h, --help              Show brief help message"
            echo "-d, --datasets          Specify the names of the datasets to download; Can be any of: Jansky, Jansky-Adrenal, Kildisiute, Dong, Bedoya, and all"
            echo "-o, --output DIR        Specify a directory to download the data to"
            exit 0
        ;;
        -d|--datasets)
            shift 
            if test $# -gt 0; then
                dats_names=()
				args=( "$@" )
				set -- "${args[@]}"
				while (( $# )); do
					if [ ${1:0:1} == "-" ]; then
						break
					fi
					#echo "Path: $1"
					dats_names+=($1)
					shift
				done
				unset args
            else
                echo "No dataset names provided"
				exit 1
            fi
        ;;
        -o|--out_path)
            shift
            if test $# -gt 0; then
                out_path=$1
            else
                echo "No output DIR specified"
                exit 1
            fi
            shift
		;;
        *)
            echo "bad option"
            exit 1
        ;;
    esac
done
###
if [ -z ${datasets+x} ]; then
	echo "No dataset names provided"
	exit 1
fi
if [ -z ${out_path+x} ]; then
	echo "No output DIR specified"
	exit 1
fi
parent="https://bimsbstatic.mdc-berlin.de/akalin/AAkalin_Anglemana/"
for x in ${datasets[@]}; do
    wget -P $out_path $parent/$x".tar.gz"
    tar -xzf $out_path/$x".tar.gz" -C $out_path
    rm $out_path/$x".tar.gz"
done