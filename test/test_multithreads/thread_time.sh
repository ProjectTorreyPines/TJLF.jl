#!/bin/bash

if [ "$#" -eq "1" ]; then
    if ! [[ "$1" =~ ^[0-9]+$ ]]; then
        echo "Give a number please"
        exit 1
    fi

    export JULIA_NUM_THREADS=$1

    
    julia time_threads.jl

elif [ "$#" -eq "2" ]; then
    if ! [[ "$1" =~ ^[0-9]+$ ]] || ! [[ "$2" =~ ^[0-9]+$ ]]; then
        echo "Give numbers please"
        exit 1
    fi

    if [ "$1" -le "$2" ]; then
	start="$1"
	end="$2"
    else
	start="$2"
	end="$1"
    fi

    for ((i=start;i<=end;i++)); do
	export JULIA_NUM_THREADS=$i
	julia time_threads.jl
    done
elif [ "$#" -gt "2" ]; then

    for i in "$@"; do
	if ! [[ "$i" =~ ^[0-9]+$ ]]; then
            echo "$i is not a number"
        else
	     export JULIA_NUM_THREADS=$i
         julia time_threads.jl
	fi
    done

else
    echo "Command line arguments"
    exit 1
fi
