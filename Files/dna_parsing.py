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

'''Scans a gene for mismatch sites and designs recombineering oligos

When a gene is being scanned, the wildtype sequence is compared to a target
sequence engineered with a stop site and a restriction site.  If the correct
mismatches can be made to the wildtype sequence, a 90mer will be generated with
the mismatched nucleotides at the center.'''

from Bio import SeqIO
from classes import *
from false_pos import *
import oligos,export


def main(geneStr,s,e,cutSites,lag,org,no_selection,debug,b=False):
    if debug:
        print "in dna_parsing.main"

    # Generate the list of cutsites
    enzL,mutL = lists(cutSites)

    # Generate a list of useable oligos
    oL = []
    for mutSite in mutL:

        # geneStr is already oriented as the lagging strand
        oligoTup = scan(geneStr,mutSite,s,e,lag,debug,b)
        if oligoTup:
            oL.extend(oligoTup)

    # Sort oligos by location in gene
    if oL:
        oL.sort()
        results = []
        enzL2 = []
        for o in oL:
            res = oligos.primerGen(o[1],org,geneStr,debug)
            if res:
                
                if no_selection: # Return everything that is generated
                    results.append(res)
                else: # Restrict to one oligo per enzyme

                    if res[0].enz[0] in enzL2:
                        pass
                    else:
                        results.append(res)
                        enzL2.append(res[0].enz[0])

                    # Quick hack to end the search early
                    if len(enzL2) == len(enzL):
                        break
                    
        return results
    else:
        if debug:
            print "Error: dna parser failed to find mutable sites"
        return False
    
def scan(geneStr,mutSite,start,end,lag,debug,b): # change name?

    results = []

    if end > len(geneStr):
        end = len(geneStr)

    while start%3:
        start -= 1 # Ensure "begin" value is in-frame

    for i in range(start,end,3):
        
        geneSub = geneStr[i:i+9]
        
        if len(geneSub) != 9:
            continue # e.g., end of the gene
        modSite = nFormat(mutSite,geneSub)

        # Compare the wt strand to the desired mutated sequence
        if compare(geneSub,modSite):
            oStr,s,e = oFormat(geneSub,modSite[0])
            
            oligo = Oligo(geneStr[i-40:i].lower() + oStr + \
                          geneStr[i+9:i+50].lower(),'Coding',s+i+1,mutSite[2])

            if lag == 'Template':
                oligo = oligo.rev_comp()
            elif lag == 'Ambiguous' and b == True:
                oligo = oligo.rev_comp()

            if len(oligo) == 90:

                # Run through the false positive module
                oligo = screen(oligo,geneStr[i-40:i+50])

                if tcChk(oligo,geneStr[i-40:i+50]):
                    results.append((oligo.loc,oligo))
    if results:
        return results
    else:
        return False
    
def lists(eL):
    '''Generate lists of enzymes and mutation sites'''

    eL = eL.split(',') # From the text entry field
    if '' in eL:
        eL.remove('')
    eL = [n.strip() for n in eL] # Remove any spaces
    sL = ['TAA','TAG','TGA']

    # Create combinations of restriction and stop sites
    mL = expansion(eL,sL)
    
    return eL,mL

def expansion(enzymes,stops):
    '''Create a list that matches enzyme cut sites with stop codons

    For each combination, add to the list the combined string, 'f' or 'b',
    depending on whether the stop was added to the front or back, and a
    tuple of the original stop and enzyme site, to report in the excel file.'''

    # NOTE: The necessity for indicating 'f' vs 'b' has been depreciated as
    # of version 3.0.
    
    mL = []
    for e in enzymes:
        for s in stops:
            
            # The sites can appear adjacent
            mL.append((s+e,'f',(e,s))) # ((s+e),'f')
            mL.append((e+s,'b',(e,s)))

            # They can overlab by 1nt, in front
            if e[0] == s[-1]:
                mL.append((s+e[1:]+'N','f',(e,s)))

            # or behind
            if e[-1] == s[0]:
                mL.append(('N' + e+s[1:],'b',(e,s)))

            # Or they can overlap by 2nt, in front
            if e[:2] == s[-2:]:
                mL.append((s+e[2:]+'NN','f',(e,s)))

            # or behind
            if e[-2:] == s[:2]:
                mL.append(('NN' + e+s[2:],'b',(e,s))) # s[-1]
    return mL

def nFormat(site,sub):
    '''Modify placeholders to reflect the wildtype sequence'''
    if 'N' in site[0]:
        pos = site[1]
        seq = sFormat(site[0],sub)
        return (seq,pos)
    else:
        return (site[0],site[1])

def sFormat(site,wt):
    '''Replace 'N' by wt sequence'''
    t = ''
    for i in range(len(site)):
        if site[i] == 'N':
            t += wt[i]
        else:
            t += site[i]
    return t            

def oFormat(sub,site):
    '''Format the oligo: lowercase with uppercase mismatches'''
    t = ''
    x = 0
    for i in range(len(sub)):
        if sub[i] == site[i]:
            t += site[i].lower()
        else:
            t += site[i].upper()
            if x==0:
                s = i # Grab the first mutation index
                x = 1
            e = i # Grab the last mutation index
    return t,s,e

def compare(sub,site):        
    '''Determine if the proper mismatches can be made to the wt strand

    Someday this could be editable by user, but since the 5 mismatches
    seem to work most efficiently, we will continue using only those
    patterns for now.'''

    def mis(s,site):
        '''Create a string of mismatched indicies'''
        m = ''
        for x in range(len(s)):
            if s[x] != site[x]:
                m += str(x)
        return m

    pL = ['01234','12345','23456','34567','45678'       # 5 mismatches
          #'0123','1234','2345','3456','4567','5678'    # 4 mismatches
          ]    
    site = site[0]
    res = mis(sub,site)
    if res in pL:
        return True
    else:
        return False

def tcChk(oligo,geneSub):
    '''Test the oligo for a t-t or c-c mismatch'''
    x = 0
    for i in range(len(oligo)):
        mut,wt = oligo.seq[i],geneSub[i]
        if mut != wt:
            if mut == 'T' and wt == 'A':
                x = 1
            elif mut == 'C' and wt == 'G':
                x = 1
    return x
