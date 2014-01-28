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

'''Output the oligos that have been generated to a csv file.'''

import csv,os,glob
from Bio import SeqIO
from classes import *
import mut_fasta

def export(oInfo,bInfo,strand,seqRec,gNum,fname,org,sp=False,b=False):
    '''Process results from the program and write them to a csv file'''

    # Gather information
    name,filename = altName(seqRec.id,fname,org,gNum,b)

    # Move to the organism's directory and write the csv file
    home = os.getcwd()
    oligoDir = 'oligos'
    try:
        os.chdir(oligoDir+'/%s'%org)
    except OSError:
        os.mkdir(oligoDir+'/%s'%org)
        os.chdir(oligoDir+'/%s'%org)

    writeFile(filename,name,oInfo,bInfo,strand,seqRec,org,sp,b)

    os.chdir(home) # Back to the PyRec root folder
    # os.chdir('Files') # Back to the PyRec directory

def writeFile(filename,name,oInfo,bInfo,strand,seqRec,org,sp,b):
    '''Generate a writer object and write data to a csv file'''

    # Open the file and create a writer object
    filename += '.csv'
    fd = open(filename,'w')
    wtr = csv.writer(fd)

    # Begin writing
    template(wtr,bInfo,seqRec,strand)

    if sp:
        special(wtr)
    else:
        num = 0
        for o in oInfo:
            if o:
                if b:
                    n = name + '_B%s'
                else:
                    n = name + '%s'

                # Use the mut_fasta module to create a
                # fasta file of the mutated gene, and a string
                # for reference in the oligo database
                e = enzLookup(o[0].enz[0])
                fasta_name = n % '_' + e
                mut_gene = mut_fasta.main(o[0],fasta_name,org)

                # Write the oligos to the csv file
                info = oligoLines(wtr,n,o,mut_gene)
                wtr.writerow(info)
                #wtr.writerows([[""],[""]])
                num += 1
    fd.close()

def template(wtr,bInfo,seqRec,strand):
    '''Information to provide for each export''' 

    # Extract the BLAST information
    num,dr,pos,match,e,bWarn = bInfo
    
    # Gene information
    geneID = seqRec.id
    geneDesc = seqRec.description
    geneSeq = seqRec.seq
    
    wtr.writerow(["Oligo and primers designed using the recombineering "\
                   "module for python."])
    wtr.writerow(["Written and designed by Gregory Moyerbrailean."])
    wtr.writerows([[""],[""]])
    wtr.writerows([["Gene Information:"],["ID:",geneID],\
                   ["Description:",geneDesc],["Length:",len(seqRec.seq)]])
    wtr.writerow([""])
    wtr.writerows([["BLAST Information:"],["Match:",match],["E-Value:",e],\
                   ["Location",num,bWarn],["Direction",dr],["Position",pos],
                   ["Lagging",strand]])
    wtr.writerows([[""],[""]])
    wtr.writerow(["Oligo Designs:"])

    # Write the header for the oligo information
    wtr.writerow(["Oligo","","Sequence","False Positive?","Mutated Gene",
                  "Length","GC%","T-anneal","Oligo fwd","Sequence","Length",
                  "GC%","T-anneal","Oligo rev","Sequence","Length","GC%",
                  "T-anneal","Product Size","","Oligo restriction","",
                  "Sequence","Length","GC%","T-anneal"])
    

def oligoLines(wtr,name,oligos,mut_gene):
    '''Format the information for each oligo generated'''

    oligo,fwdP,mP,revP,pcr_len = oligos
    e = enzLookup(oligo.enz[0])

    oN = name % '_' + e
    fwdN = name % '_fwd_' + e
    mpN = name % '_mama_' + e
    revN = name % '_rev_' + e

##    mut_made_at = ["","Mutation made at nt: %s" % oligo.loc]
##    header_line = ["Oligo","Sequence","",
##                   "Forward","Sequence","Length","T-anneal",
##                   "MAMA","Sequence","Length","T-anneal",
##                   "Reverse","Sequence","Length","T-anneal",
##                   "Product"]

    if type(fwdP)==Primer:
        fwdL = [fwdN,fwdP.seq.lower(),len(fwdP),fwdP.GC(),fwdP.Tm()]
    else:
        fwdL = [fwdN,fwdP,"","",""]

    if type(mP)==Primer:
        mpL = [mpN,mP.seq.lower(),len(mP),mP.GC(),mP.Tm()]
    else:
        mpL = [mpN,mP,"","",""]
        
    if type(revP)==Primer:
        revL = [revN,"",revP.seq.lower(),len(revP),revP.GC(),revP.Tm()]
    else:
        revL = [revN,revP,"","","",""]

    oligo_line = [oN,"",oligo.seq,oligo.falsepos,mut_gene,len(oligo),oligo.GC(),
                  oligo.Tm()]
    #oligo_line = [oN,"",oligo.seq,"","Gene","","",""]
    oligo_line.extend(fwdL)
    oligo_line.extend(mpL)
    oligo_line.extend([pcr_len,""])
    oligo_line.extend(revL)

    #return mut_made_at,header_line,oligo_line
    return oligo_line

def geneDict(org):
    '''Generate a dictionary of geneID's and S-codes'''
##    fLoc = '../../Files/genome/%s/geneIDs.txt'%org
    fLoc = 'Files/genome/%s/geneIDs.txt'%org
    fd = open(fLoc)
    r = fd.read()
    fd.close()
    gLst = r.split('\r')
    geneD = {}
    for i in gLst:
        try:
            ID,S = i.split('\t')
        except ValueError:
            print i.split('\t')
        else:
            geneD[ID] = S
    return geneD

def geneLookup(ID,fname,org,gNum):
    '''Look up the naming scheme for the organism'''
    try:
        geneD = geneDict(org)
        name = geneD[ID]
    except KeyError:
        if gNum == 1:
            # Use the fasta file name as the output file name if only one gene
            name = fname.split('/')[-1] 
            name = name.replace('.fasta','')
        else:
            # If < one genes, use the geneID for each gene
            name = ID
    return name

def altName(ID,fname,org,gNum,b):
    name = geneLookup(ID,fname,org,gNum)
    filename = name
    if b:
        filename += '_B'
    return name,filename

def openWtr(filename):
    filename += '.csv'
    fd = open(filename,'w')
    wtr = csv.writer(fd)
    return wtr

def enzLookup(enzSeq):
    '''Return an abbreviation for a given enzyme'''
    if enzSeq == 'GGATCC':
        return 'B'
    elif enzSeq == 'GAATTC':
        return 'E'
    elif enzSeq == 'AAGCTT':
        return 'H'
    else:
        return enzSeq

def special(wtr):
    '''Information to write when oligos cannot be generated'''
    wtr.writerows([["**** I'm sorry, no oligos could be designed with the"\
                   " current specifications. ****"],
                    ["**** Please consult the gene sequence or modify the "\
                     "specifications"]])
