#!/usr/bin/env bash

time for line in $(cat all_recount3.txt); do esearch -db sra -query $line  | efetch -format runinfo | cut -d ',' -f 1 | grep -v "Run" ; done > recount3.txt
