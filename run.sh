#! /bin/bash

echo "This script will run all the track analysis in the background (nohup), and it is hard to be canceled."
read -p "Are you sure you want to continue? (y/N): " choice

if [[ "$choice" != "y" && "$choice" != "Y" ]]; then
    echo "Operation canceled."
    exit 1
fi

if [ -d "log" ]; then
    rm -rf log/*
else
    mkdir log
fi

for i in {1..12}; do
    data_file=$(printf "Data/DAT%06d" $i)
    log_file="log/DAT$(printf "%06d" $i).log"
    
    nohup time ./test "$data_file" > "$log_file" 2>&1 &
done

