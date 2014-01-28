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

'''Locate wildtype sequences near the mutation that may cause false positives.

If a similar sequence to the mutation can be found in the nearby wildtype
sequence, there is a chance that the MAMA-PCR will give a false-positive result.
This module scans for those sequences and notes if one is found.'''

class FalsePositiveError(Exception):
    def __init__(self,value):
        self.value = value
    def __str__(self):
        return self.ValStr()
    def ValStr(self):
        if self.value == 1:
            return "Warning: Entire mismatch found in nearby wildtype sequence"
        elif self.value == 2:
            return "Warning: First or last 4nt of mismatch found "\
                   "in nearby wildtype sequence"
        elif self.value == 3:
            return "Warning: First or last 3nt of mismatch found in "\
                   "nearby wildtype sequence"
        elif self.value == 4:
            # "catch-all" error
            return "Warning: Similar sequence to mutation detected in "\
                   "nearby wildtype sequence"

def screen(oligo,wt_seq_str):
    mut_int = oligo.mut
    mut_str = oligo.seq[mut_int:mut_int+5]
    wt_sub_str = wt_seq_str[mut_int-10:mut_int+15]

    try:
        for i in range(len(wt_sub_str)-5):
            wt_temp = wt_sub_str[i:i+5]

        # Determine if the same mismatched nucleotides are present
        # in the nearby wildtype sequence
            if wt_temp == mut_str:
                raise FalsePositiveError(1)

        # Determine if any 3 or any 4 mismatched nucleotides are present
        # in the surrounding wildtype sequence
            elif wt_temp[:4] == mut_str[:4] or wt_temp[-4:] == mut_str[-4:]:
                raise FalsePositiveError(2)
            elif wt_temp[:3] == mut_str[:3] or wt_temp[-3:] == mut_str[-3:]:
                raise FalsePositiveError(3)

        # Determine if there is only a 1 or 2 nt difference
            diff = 0
            for nt in range(5):
                if wt_temp[nt] != mut_str[nt]:
                    diff += 1
            if diff == 1 or diff == 2:
                raise FalsePositiveError(4)

    except FalsePositiveError,e:
        oligo.falsepos = e
        return oligo
    else:
        return oligo
