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

'''Handles the local BLAST and the parsing of the results.

The BLAST uses the NCBI blast+ command-line tools to run a local BLAST against
the organism's genome. In the event that a closed genome is not available for
a species, the genome of a closely related strain can be used in its place.
When a hit has been found, the parser function will extract and return relevant
information regarding the corresponding gene.

Alternatively, the user may specify to disable the BLAST function. In this case,
the module will use the scaffold files to extract the necessary information.
The user therefore does not have to have the blast+ command-line tools.
However, the user will also not be able to run organisms such as L. reuteri
against a similar genome, as this method requires exact gene matches.'''

import subprocess
from Bio.Blast.Applications import NcbiblastnCommandline as ncl
from Bio.Blast import NCBIXML as nxml

def BlastGenome(queryFile,genome,debug,outputFile='Files/extras/temp_blast.xml'):
    if debug:
        print "In BLASTing.BlastGenome"

    # Modify the genome filename to reflect the path to the genome
    genome = genome.replace(' ','')
    genomePath = 'Files/genome/' + genome + '/' + genome

    ## Call blast+ from python
    cline = ncl(query=queryFile,db=genomePath,out=outputFile,outfmt=5)
    ret_code = subprocess.call(str(cline),shell=True)

    if ret_code:
        print 'BLASTing file "%s" returned error code %s' % (queryFile,ret_code)

    temp = open(queryFile).read()
    geneID = temp.split()[0]
    geneID = geneID.lstrip('>')
    result = nxml.read(open(outputFile))
    
    # If the blast returns no results, it will be treated as a gene
    # in the ambiguous region and oligos will be made from both strands
    if result.alignments:
        return parseRecord(result,genomePath,debug)
    else:
        return 0,0,'Ambiguous','No Match','N/A'

def parseRecord(xmlfile,genomePath,debug):
    if debug:
        print "In BLASTing.parseRecord"
    
    result = nxml.read(open('Files/extras/temp_blast.xml'))
    hit = result.alignments[0].hit_def
    e = result.descriptions[0].e
    if debug:
        print "Blast match: ",hit
        print "E-value: ",e
        
    hitL = hit.split()
    hitID = hitL[0]
    t = [n for n in hitL if '..' in n]
    hitInfo = t[0]
    num1,num2 = hitInfo.split('..')
    num2 = num2[:num2.find('(')]
    num1,num2 = int(num1),int(num2)
    strand = hitInfo[hitInfo.find('('):]

    
    # Determine the direction, relative location, and position of the gene
    direction = getDirection(hitInfo)
    termUpper,termLower = getRelativeLocation(genomePath)
    pos = getLocation(num1,termUpper,termLower)

    # TODO
    # Integrate warning for multiple hits
    
    return num1,direction,pos,hit,e,''

def SearchGenome(queryFile,genomeName,debug):
    from Bio import SeqIO

    genomePath = 'Files/genome/'+genomeName+'/'+genomeName
    genome = openGenome(genomePath)
    high,low = getRelativeLocation(genomePath)
    
    gene = SeqIO.read(open(queryFile),'fasta')
    geneStr = str(gene.seq)
    geneComp = str(gene.seq.reverse_complement())
    count = 0

    if geneStr in genome:
        direction = 'forward'
        n = genome.find(geneStr)
        pos = getLocation(n,high,low)
        count += genome.count(geneStr)
    elif geneComp in genome:
        direction = 'reverse'
        n = genome.find(geneComp)
        pos = getLocation(n,high,low)
        count += genome.count(geneComp)
    else:
        return 0,0,'Ambiguous','No Match','N/A',''

    # If the gene sequence is present more than once, issue a warning
    bWarn = 'Warning: Gene sequence detected multiple times in genome'
        
    return n,direction,pos,'No BLAST data','No BLAST data',bWarn

def getRelativeLocation(genomePath):
    l,t = getTermRegion(genomePath+'.txt')
    buff = 0.05 * l
    high = t + buff
    low = t - buff
    return high,low 

def getTermRegion(path):
    fd = open(path)
    info = fd.read()
    l,t = info.split('\n')
    l,t = int(l),int(t)
    return l,t    

def getDirection(line):
    if '(+)' in line:
        d = 'forward'
    elif '(-)' in line:
        d = 'reverse'
    return d

def getLocation(num,high,low):
    if num < low:
        p = 'Early'
    elif num > high:
        p = 'Late'
    else:
        p = 'Ambiguous'
    return p

def openGenome(gpath):
    fd = open(gpath+'.fasta')
    g = fd.read()
    g = g.replace('\n','')
    return g
