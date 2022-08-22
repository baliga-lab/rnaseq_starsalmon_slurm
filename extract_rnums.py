#!/usr/bin/env python3

INFILE = "Past_PoritesAstreoides_metadata.csv"

with open(INFILE) as infile:
    infile.readline()
    for line in infile:
        comps = line.strip().split(",")
        rnum = comps[-1].replace('"', '')
        print(rnum)
