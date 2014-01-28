'''Place this file in the oligos folder you wish to count, and run it.
The program will return the number of oligos that could vs could not be
generated'''

import glob

files = glob.glob('*.csv')
y,n = 0,0
for f in files:
    csv = open(f).read()
    if "I'm sorry" in csv:
        n += 1
    else:
        y += 1



y1 = "(%d)" % (y*100/(y+n))
n1 = "(%d)" % (n*100/(y+n))


print "Genes with Oligos:",y,y1
print "Genes without:",n,n1
print "Total:",y+n
