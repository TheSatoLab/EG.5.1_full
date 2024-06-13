#!/usr/bin/env python

import sys, re
argvs = sys.argv

f = open(argvs[1])

f.__next__()

print("Id\tmut")
Id_mut_d = {}
for line in f:
  line = line.strip().split("\t")
  Id = line[1]
  # Substitution
  mut_line = line[29].replace('(','').replace(')','')
  if mut_line != '':
    mut_l = mut_line.split(',')
    for mut in mut_l:
      mut = mut.replace('*','stop') # Stop mutation
      print("\t".join([Id,mut]))
  # Deletion
  mut_line = line[30].replace('(','').replace(')','')
  if mut_line != '':
    mut_l = mut_line.split(',')
    for mut in mut_l:
      mut = mut.replace('-','del')
      print("\t".join([Id,mut]))
  # Insertion
  mut_line = line[31].replace('(','').replace(')','')
  if mut_line != '':
    mut_l = mut_line.split(',')
    for mut in mut_l:
      mut = mut.split(':')[0]+':ins'+mut.split(':')[1]+mut.split(':')[2]
      print("\t".join([Id,mut]))
