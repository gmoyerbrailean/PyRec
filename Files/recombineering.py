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

'''This module facilitates oligo generation for DNA recombineering

For each gene given, this module delegates the task and passes information
between the other program modules.'''

from Bio import SeqIO
import BLASTing as blst
import subprocess as sp
import dna_parsing,oligos,export

def iterSeqs(f_name,beg,end,org,cutSites,no_select,debug=False,SKIP_BLAST=True):
    if debug:
        print "In recombineering.iterSeqs"

    geneCnt = 0
    f_type = f_name.split('.')[-1] # Grab the file extension
    handle = open(f_name) # handle = open('../'+f_name)
    print
    print "Running file:",f_name
    genes_per_file = handle.read().count('>')
    if debug:
        print "File type:",f_type
        print "Number of genes found:",genes_per_file

    handle.seek(0) # Back to beginning of file
    seqParser = SeqIO.parse(handle,f_type)
    
    for gene in seqParser:
        geneID = gene.description.split()[0]

        # Create a fasta record to BLAST the genome with
        tempFile = open('Files/extras/temp.fasta','w')
        tempStr = ">%s\n%s" % (gene.description,gene.seq.upper())
        tempFile.write(tempStr)
        tempFile.close()

        # # If the genome is not closed, we must BLAST the genes
        # # against a similar strain
        # if org in openGenomes:
        #     SKIP_BLAST = False

        # Extract gene information from the BLASTing module
        if SKIP_BLAST:
            # If the organism has a closed genome, the local BLAST can
            # be skipped, thus avoiding the blast+ tools requirement
            n,d,p,m,e,w = blst.SearchGenome('Files/extras/temp.fasta',org,debug)
        else:
            n,d,p,m,e,w = blst.BlastGenome('Files/extras/temp.fasta',org,debug)
        blstInfo = [n,d,p,m,e,w]

        lag = getLagging(d,p)

        if debug:
            print
            print "*"*80
            print
            print "Gene ID:",geneID
            print "Gene length:",len(str(gene.seq))
            print "Num:",n
            print "Dir:",d
            print "Pos:",p
            print
            print "Parsing from",beg,"to",end

        Restrict(f_name,genes_per_file,gene,lag,blstInfo,beg,end,
                 cutSites,org,no_select,debug)
        
        # Count the overall number of genes parsed        
        geneCnt += 1
    print
    print "Oligo generation complete. Genes found: %d" % geneCnt
    print "If you are finished press Quit to exit, otherwise you may "
    print "continue to generate oligos."

def Restrict(f_name,gNum,gene,lag,blstInfo,
             beg,end,cutSites,org,no_selection,debug):

    geneStr = str(gene.seq).upper()
##    geneComp = str(gene.seq.reverse_complement().upper())
    geneComp = str(gene.seq.complement().upper())

    if lag == 'Ambiguous':

        # Create one file for coding oligos
        oligos = dna_parsing.main(geneStr,beg,end,cutSites,lag,
                                  org,no_selection,debug)
        if oligos and set(oligos) != set([False]):
            export.export(oligos,blstInfo,lag,gene,gNum,f_name,org)
        else:
            export.export(oligos,blstInfo,lag,gene,gNum,f_name,org,sp=True)

        # And one file for the template oligos
        ##oligos = dna_parsing.main(geneComp,beg,end,cutSites,lag,
        oligos = dna_parsing.main(geneStr,beg,end,cutSites,lag,
                                  org,no_selection,debug,b=True)
        if oligos and set(oligos) != set([False]):
            export.export(oligos,blstInfo,lag,gene,gNum,f_name,org,b=True)
        else:
            export.export(oligos,blstInfo,lag,gene,gNum,f_name,org,b=True,
                          sp=True)

    else:
        # Choose the sequence to use based on the lagging strand
##        if lag == 'Coding':
##            sequence = geneStr
##        else:
##            sequence = geneComp
        sequence = geneStr
            
        oligos = dna_parsing.main(sequence,beg,end,cutSites,lag,
                                  org,no_selection,debug)
        if oligos and set(oligos) != set([False]):
            export.export(oligos,blstInfo,lag,gene,gNum,f_name,org)
        else:
            export.export(oligos,blstInfo,lag,gene,gNum,f_name,org,sp=True)

def getLagging(d,p):
    '''Determine the lagging strand at current position'''
    if d == 'forward':
        if p == 'Early':
            lag = 'Template'
        elif p == 'Late':
            lag = 'Coding'
        else:
            lag = "Ambiguous"
    else:
        if p == 'Early':
            lag = 'Coding'
        elif p == 'Late':
            lag = 'Template'
        else:
            lag = "Ambiguous"
    return lag
