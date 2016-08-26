__author__ = 'Debha'

import csv
import numpy as np
from operator import itemgetter


class Mutation:
    def __init__(self, sub, chr, start, end, ref, seq, type, gene, length):
        self.sub = sub
        self.chr = chr
        self.start = start
        self.end = end
        self.ref = ref
        self.seq = seq
        self.type = type
        self.gene = gene
        self.length = length

def read_data(file):
    f = open(file)
    csv_f = csv.reader(f)
    next(csv_f, None) #skip the header
    data = {}

    for row in csv_f:
        subject = row[6].strip("-")
        chr = row[0]
        start = int(row[1])
        end = int(row[2])
        ref = row[3]
        seq = row[4]
        length = 0

        if len(ref) == len(seq): mut_type = 'snv'
        else: mut_type = 'indel'

        mut = Mutation(subject, chr, start, end, ref, seq, mut_type, '', length)

        if subject not in data:
            data[subject] = [mut]
        else:
            data[subject].append(mut)

    f.close()

    return data

def read_map(file):
    f = open(file)
    csv_f = csv.reader(f)
    next(csv_f, None) #skip the header
    gene_map = {}

    for row in csv_f:
        name = row[4] #use 0 for ensembleID and 4 for gene name
        chr = row[1]
        start = int(row[2])
        end = int(row[3])

        if chr not in gene_map:
            gene_map[chr] = [[name, start, end]]
        else:
            gene_map[chr].append([name, start, end])

    for chr in gene_map:
        gene_map[chr] = sorted(gene_map[chr], key=itemgetter(1))

    return gene_map

def gene_muts(data, gene_map):
    gene_mut = {}
    for sub in data:
        for m in data[sub]:
            gene_name = ''
            for gene in gene_map[m.chr]:
                if m.start >= gene[1] and m.end <= gene[2]:
                    gene_name = gene[0]
                    m.gene = gene_name #add gene name to data
                    m.length = gene[2]-gene[1]
                    break

            if gene_name not in gene_mut and gene_name != '':
                if m.type == 'snv':
                    gene_mut[gene_name] = [1, 0]
                if m.type == 'indel':
                    gene_mut[gene_name] = [0, 1]
            elif gene_name in gene_mut:
                if m.type == 'snv':
                    gene_mut[gene_name][0] += 1
                if m.type == 'indel':
                    gene_mut[gene_name][1] += 1

    return gene_mut

def make_GO_map(file):
    f = open(file, mode="r")
    GO_map = {}

    for line in f:
        split_line = line.split("\t")
        for n in range(len(split_line)):
            GO_term = split_line[n]
            if "\n" in GO_term:
                GO_term = GO_term.rstrip("\n")
            if n > 1:
                try:
                    GO_map[split_line[0]].append(GO_term)
                except KeyError:
                    GO_map[split_line[0]] = [GO_term]
    f.close()

    return GO_map




