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

'''A container for the classes used by multiple programs to handle DNA'''

class Oligo(object):
    '''Initialize an oligo object'''
    def __init__(self,sequence,strand,loc,enzyme,fp=''):
        self.seq = sequence
        self.strand = strand
        self.loc = loc # Oligo location in the gene
        self.mut = self.mark(sequence) # Mutation location in the oligo
        self.enz = enzyme
        self.falsepos = fp
    def __str__(self):
        return self.seq
    def __len__(self):
        return len(self.seq)
    def __getslice__(self,i,j):
        return self.seq[i:j]
    # Decided to use Bio.Seq methods
##    def __reversed__(self):
##        return Oligo(self.seq[::-1])

    def mark(self,seq):
        '''Return the location of the first mismatch'''

        for i in range(len(seq)):
            if seq[i].isupper():
                return i

    def comp(self):
        '''Return the complement of the oligo

        Creates a Biopython Seq object, and uses the
        Seq object complement method'''
        from Bio.Seq import Seq
        from Bio.Alphabet import IUPAC
        s = Seq(self.seq,IUPAC.unambiguous_dna)
        s = s.complement()
        if self.strand in ['coding','Coding']:
            strand = 'Template'
        else:
            strand = 'Coding'
        return Oligo(str(s),strand,self.loc,self.enz)

    def rev_comp(self):
        '''Return the reverse complement of the oligo

        Creates a Biopython Seq object, and uses the
        Seq object reverse_complement method'''
        from Bio.Seq import Seq
        from Bio.Alphabet import IUPAC
        s = Seq(self.seq,IUPAC.unambiguous_dna)
        s = s.reverse_complement()
        if self.strand in ['coding','Coding']:
            strand = 'Template'
        else:
            strand = 'Coding'
        return Oligo(str(s),strand,self.loc,self.enz)

    def rev(self):
        '''Return the reverse sequence of the oligo'''
        return Oligo(self.seq[::-1],self.strand,self.loc,self.enz)

    def GC(self):
        '''Return the % G-C content of the oligo'''
        temp_seq = self.seq.upper()
        GC = ((temp_seq.count('G')+temp_seq.count('C')) \
              /float(len(temp_seq))) * 100
        return GC

    def Tm(self):
        '''Return the T-anneal for the primer object'''
        temp_seq = self.seq.upper()
        GC = ((temp_seq.count('G')+temp_seq.count('C')) \
              /float(len(self.seq))) * 100
        Tm = (0.41*GC)+(69.3-(650/float(len(temp_seq))))
        return Tm

class Primer(object):
    '''Initialize a primer object'''
    def __init__(self,sequence,strand,location=0):
        self.seq = sequence
        self.strand = strand
        self.location = location
    def __str__(self):
        return self.seq
    def __len__(self):
        return len(self.seq)
    def __reversed__(self):
        return Primer(self.seq[::-1])

    def comp(self):
        '''Return the complement of the primer

        Creates a Biopython Seq object, and uses the
        Seq object complement method'''
        from Bio.Seq import Seq
        from Bio.Alphabet import IUPAC
        s = Seq(self.seq,IUPAC.unambiguous_dna)
        s = s.complement()
        return Primer(str(s),self.strand,self.location)

    def rev_comp(self):
        '''Return the reverse complement of the primer

        Creates a Biopython Seq object, and uses the
        Seq object reverse_complement method'''
        from Bio.Seq import Seq
        from Bio.Alphabet import IUPAC
        s = Seq(self.seq,IUPAC.unambiguous_dna)
        s = s.reverse_complement()
        return Primer(str(s),self.strand,self.location)

    def rev(self):
        '''Reverse the orientation of the primer'''
        return Primer(self.seq[::-1],self.strand,self.location)

    def GC(self):
        '''Return the % G-C content of the oligo'''
        temp_seq = self.seq.upper()
        GC = ((temp_seq.count('G')+temp_seq.count('C')) \
              /float(len(temp_seq))) * 100
        return GC

    def Tm(self):
        '''Return the T-anneal for the primer object'''
        temp_seq = self.seq.upper()
        GC = ((temp_seq.count('G')+temp_seq.count('C')) \
              /float(len(self.seq))) * 100
        Tm = (0.41*GC)+(69.3-(650/float(len(temp_seq))))
        return Tm

class DNA_List(list):
    '''Initializes a list for DNA sequences

    Has all the properties of a regular list, with two additions:
    t_60 will return the sequence in the list closest to 60'C, and
    t_match will do the same with a specified temperature'''
    def __init__(self):
        list.__init__(self)
    def t_60(self):
        diffs = [abs(n.Tm()-60) for n in self]
        indx = diffs.index(min(diffs))
        return self[indx]
    def t_match(self,tm):
        diffs = [abs(n.Tm()-tm) for n in self]
        indx = diffs.index(min(diffs))
        return self[indx]

class CodonError(Exception):
    def __init__(self,value):
        self.value = value
    def __str__(self):
        return repr(self.value)

class Codon(object):
    '''Initialize a Codon class for nucleotides'''
    def __init__(self,splice):
        '''Initialize a Codon class for nucleotides'''
        if len(splice) == 3:
            self.nt1 = splice[0].upper()
            self.nt2 = splice[1].upper()
            self.nt3 = splice[2].upper()
            self.codon = self.nt1 + self.nt2 + self.nt3
        else:
            raise CodonError(1)
    def __str__(self):
        '''Prints the codon instance'''
        return '"%s"' % self.codon
    def __eq__(self,other):
        '''Checks equality of two codons'''
        if type(other) == str:
            if self.codon == other.upper():
                return True
            else:
                return False
        if type(other) == Codon:
            if self.codon == other.codon:
                return True
            else:
                return False

    def mut(self,first,last):
        '''C.mut(sub start, sub end) --> mutated codon

        Used when the codon is the non-stop codon being mutated, and the
        sequence does not matter. Will always attempt to make as many
        t-t or c-c mutations as it can. Works for both the coding and
        template strands.'''
        t = [self.nt1,self.nt2,self.nt3]
        for i in range(first,last+1):
            if t[i] == 'A':
                t[i] = 'T'
            elif t[i] == 'G':
                t[i] = 'C'
            else: # t[i] in 'CT':
                t[i] = 'G'
        m = ''.join(t)
        return m
