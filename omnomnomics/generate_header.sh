#!/bin/bash
cat TAGDIRlist.txt | xargs -l basename | cut -f "$1" -d "$2" | awk 'BEGIN{ORS="\t"}{print $0}' > clean.header.tmp
