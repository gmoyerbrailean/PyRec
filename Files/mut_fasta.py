####### LICENSE #######
# This code is part of the Recombineering module, written by Gregory
# Moyerbrailean at Michigan State University, Department of Microbiology 
# and Molecular Genetics. 
# Copyright (C) 2010 Gregory Moyerbrailean
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

'''Output a fasta file with the mutated gene sequence for each mutation'''

import os,glob
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def main(oligo,fname,org):    
    # Step 1: Open temp.fasta and extract gene information
    handle = open('../../Files/extras/temp.fasta')
    wt_seq,gene_id,gene_desc = fasta_to_gene(handle)

    # Step 2: find oligo[:34] or so in temp.fasta
    wt_str = str(wt_seq)
    loc = wt_str.find(oligo.seq[:34].upper())

    if loc < 0:
        oligo = oligo.rev_comp()
        loc = wt_str.find(oligo.seq[:34].upper())
    
    # Step 3: replace part of the wt sequence with the oligo
    mt_str = wt_str[:loc] + oligo.seq + wt_str[loc+len(oligo.seq):]

    # Step 4: Generate fasta and output to a mutated fasta directory
    fasta_str = gene_to_fasta(mt_str,gene_id,gene_desc)

    os.chdir('../../mutated_genes/')
    fname = fname +  '.fasta'

    # Only write 1 file per enzyme per oligo
    files = glob.glob('*')
    if not fname in files:
        fd = open(fname,'w')
        fd.write(fasta_str)
    else:
        pass
    
    os.chdir('../oligos/%s'%org)

    return mt_str

def fasta_to_gene(handle):
    '''Convert a fasta file to a gene string'''
    
    gene = SeqIO.read(handle,'fasta')
    gene_seq = gene.seq
    gene_id = str(gene.id)
    gene_desc = str(gene.description)
    gene_desc = gene_desc.replace(gene_id+' ','') # Remove the ID
    return gene_seq,gene_id,gene_desc

def gene_to_fasta(gene_str,gene_id,gene_desc):
    '''Convert a gene string to a fasta formatted string'''

    gene_seq = Seq(gene_str)
    gene_seq_r = SeqRecord(gene_seq,id=gene_id,description=gene_desc)
    fasta_str = gene_seq_r.format('fasta')
    return fasta_str
