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

'''Create pcr primers based on the oligos generated in dna_parsing.

When a 90mer recombineering oligo is generated, the next step is to design
three primer oligos: a forward primer, a reverse primer, and a MAMA-PCR primer.
The forward and reverse need to be approximately 500nt away from the
recombineering oligo, and so the primers are generated using a genomic scaffold
sequence. This can be a whole genome sequence for a closed genome, or a series
of contig files for an open genome. The MAMA primer is generated from the
recombineering oligo sequence. The primers are designed between 15 and 30 nt
in length, and to have similar melting temperatures for optimal PCR runs.'''

from Bio import SeqIO
from classes import *

def primerGen(oligo,org,geneStr,debug):
    if debug:
        print "in oligos.primerGen"
        print "Oligo",oligo.seq

    f = searchScaffolds(oligo,oligo.enz,org,debug)

    if f:
        num,subScaff,dr = f # Stop passing 'num'?        
        fwdP = upstream(subScaff,debug)

        if fwdP:

            # The recombineering oligo will be oriented 5' - 3'
            mamaP = MAMAprimer(oligo,fwdP,dr,debug)
            if not mamaP:
                mamaP = "Error: could not design reverse primer"
            else:
                # Depending on dr, re-orient to 5' - 3'
                if dr == '+':
                    #mamaP = mamaP.rev_comp()
                    mamaP = mamaP.comp()

            # Generate 'downstream'
            revP = downstream(fwdP,subScaff,debug)

            if revP:
                revP = revP.rev_comp()
                pcr_len = revP.location - fwdP.location + len(revP)
            else:
                revP = "Error: could not design reverse primer"
                pcr_len = 0

            return oligo,fwdP,mamaP,revP,pcr_len

        else:
            if debug:
                print "Error: Could not design forward primer"
            return False
    else:
        if debug:
            print "Error: oligo failed scaffold search"
        return False

def searchScaffolds(oligo,cutSite,org,debug):
    if debug:
        print "in oligos.searchScaffolds"

    # Based on the organism selected, determine the correct scaffold
    sFile = 'Files/genome/' + org + '/' + org + '.fasta'

    seqHandle = open(sFile)
    parser = SeqIO.parse(seqHandle,'fasta')

    # Manipulate the oligos instead of the contigs to increase efficiency
    a = oligo.seq[:34].upper()
    b = oligo.rev_comp().seq[-34:].upper()

    cntrl,dr = quickCheck(a,b,sFile,debug)
    
    if cntrl:
        for scaff in parser:
            sStr = str(scaff.seq)
            if dr == '+':
                num = sStr.find(oligo.seq[:34].upper())
            else:
                num = sStr.find(oligo.rev_comp().seq[:34].upper())

            if num > 0:
                if 'contig' in scaff.description or \
                   'Contig' in scaff.description:
                    pass

                # If the scaffold is a closed genome (in contrast to a contig),
                # the circular nature must be accounted for.
                else:
                    sStr = circularize(sStr)
                    num += 1000 # To compensate

                subS = createSubScaff(num,sStr,cutSite)
                if subS:
                    return (num,subS,dr)
                else:
                    if debug:
                        print "Error: cutsite not unique"
                    return False
                
        # If nothing was returned after searching the scaffold(s)
        if debug:
            print "Error: Oligo not found in scaffold"
        return False
    else:
        if debug:
            print oligo.seq
            print "Error: Oligo failed scaffold check"
        return False
                

def quickCheck(a,b,sFile,debug):
    '''Quickly search through scaffold strings

    The goal is to quickly ascertain wether the target sequence is present once,
    more than once, or not at all before committing to search through all the
    individual scaffolds. Also, allows the assurance that the correct scaffold is
    being used, and that the sequence that was searched for is uniquie'''
    if debug:
        print "in oligos.quickCheck"

    handle = open(sFile)
    allScaffs = handle.read()

    allScaffs = allScaffs.replace('\n','')
    
    handle.close()
    aN = allScaffs.find(a);aC = allScaffs.count(a)
    bN = allScaffs.find(b);bC = allScaffs.count(b)
    
    x=0,0
    if aN>0:
        x = 1,'+'
    elif bN>0:
        x = 1,'-'
        
##    # Make sure the sequence is unique
##    if aC>1 or bC>1:
##        x = 0,0
##    if aN > 0 and bN > 0:
##        x = 0,0
        
    return x

def circularize(scaffStr):
    '''Modify the scaffolds to account for the circular genomes

    Take the "first" and "last" 1000 nucleotides and add them to the "end"
    and "beginning" of the genome, respectively. This is only applicable
    with a closed circular genome, and not for contig scaffolds.'''

    front = scaffStr[:1000]
    rear = scaffStr[-1000:]

    return rear + scaffStr + front

def createSubScaff(num,sStr,cut):
    '''Search 600 bp up- and downstream from the mutation site.

    If the cutsite is present in the wt, move on to the next oligo. The only
    exception to this rule is if the cutsite is within the first or last
    100nt. In this case, the wt digest would still be distinguishable from the
    mutant digest.'''

    mid = num + 45 # Middle of oligo
    factor = 550

    while factor <= 650:
        u = mid - factor # 600, 500, 400
        d = mid + factor

        # This (below) may occur when using contigs
        if u < 0 or d > len(sStr):
            return False
        
        subScaff = sStr[u:d]
        if uniqueCutSite(subScaff,cut):
            return subScaff
        factor += 50 
    
    return False # cutSite not unique

def uniqueCutSite(sub,cut):
    '''Determine if the cutsite sequence is unique'''

    if cut[0] in sub[100:-100]:
        return False
    else:
        return True

def upstream(subScaff,debug):
    '''Searches about 500nt upstream of the mutated site for a primer'''
    subSeq = subScaff[0:100]
    primerL = primerSearch(0,100,subSeq)
    if primerL:
        fwdP = primerL.t_60()
        fwdP.location = subScaff.find(fwdP.seq)
        return fwdP
    else:
        return 0
    
def downstream(fwdP,subScaff,debug):
    '''Searches about 500bp downstream of the mutated site for a primer'''
    subSeq = subScaff[len(subScaff)-100:len(subScaff)]
    primerL = primerSearch(0,100,subSeq)
    if primerL:
        revP = primerL.t_match(fwdP.Tm()) # Need similar Tm to fwdP
        revP.location = subScaff.find(revP.seq)
        return revP
    else:
        return 0

def MAMAprimer(oligo,fwdP,dr,debug):
    '''Searches for a MAMA-PCR screening primer'''

    if dr == '+':
        subOligo = oligo.seq[oligo.mut:oligo.mut+30]
    else:
        
        # This is ONLY true for 5 consecutive mismatches
        # If any other mismatch system is used, this will need updating.
        subOligo = oligo.seq[oligo.mut-25:oligo.mut+5]

        # Take the reverse of the oligo--mamaSearch will shorten the
        # 5' (non-mismatch) end as necessary
        subOligo = subOligo[::-1]
        
    subOligo = subOligo.upper()
    primerL = mamaSearch(subOligo)

    if primerL:
        mamaP = primerL.t_match(fwdP.Tm()) # Need similar Tm to fwdP
        tOligo = oligo.seq.upper()
        mamaP.location = tOligo.find(mamaP.seq.upper())
        return mamaP
    else:
        return 0
        
def primerSearch(start,length,sub):
    '''Search a specified sequence for potential primer sequences'''
    primers = DNA_List()
    if length >= 15:
        for i in range(length):
            for n in range(15,31):
                if i+n < length:
                    primer = Primer(sub[i:i+n],'',start+i)
                    if 55 <= primer.Tm() < 65: # Ideal is 60, +-5 is acceptable
                        primers.append(primer)
                    # Quit early if a near-60 primer is found  
                    if 59.9 < primer.Tm() < 60.1:
                        return primers
    if primers:
        return primers
    else:
        return 0

def mamaSearch(sub):
    '''Modified version of primerSearch'''
    primers = DNA_List()
    if len(sub) >= 15:
        for n in range(15,31):
            primer = Primer(sub[:n][::-1],'') # [::-1] -> re-orient to 5' - 3'
            if 55 <= primer.Tm() < 65:
                primers.append(primer)
    return primers
