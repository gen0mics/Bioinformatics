#!/usr/bin/python -tt
# Needleman-Wunch Algorithm with specific penalties for 5', 3' and internal gaps

import sys, math, re


#### reads in fasta file
def read_fasta(fp):
     name, seq = None, []
     for line in fp:
         line = line.rstrip()
         if line.startswith(">"):
             if name:    yield (name, ''.join(seq))
             name, seq = line, []
         else:
             seq.append(line)
     if name: yield (name, ''.join(seq))

with open(sys.argv[1]) as fp:
    seqs = []
    for name, seq in read_fasta(fp):
        seqs.append(seq)

#### These are the sequences to be analyzed
seq1 = seqs[0]
seq2 = seqs[1]

def global_gap():
    #### gap policy
    def gap(x, y, y1):
        if (x == y):
            return gap3
        elif ( x == y1):
            return gap5
        else:
            return gapInt

    #### default gap value
    gap5 = 3.
    gap3 = 4.
    match = 0.
    mismatch = 3.
    gapInt = 5.

    #### print input
    print "Needleman-Wunch Algorithm"
    print "SEQUENCE 1:", seq1; print "SEQUENCE 2:", seq2
    print "5' gap          : ", gap5
    print "3' gap          : ", gap3
    print "gap penalty     : ", gapInt
    print "match score     : ", match
    print "mismatch score  : ", mismatch 

    #### Initiate and calculate value
    lseq1 = len(seq1); lseq2 = len(seq2)
    row = lseq2+1; col = lseq1+1
    
    val = []
    for i in range(row):
        val.append([i*gap5]) # gap value for the first column

    for j in range(1,col):
        val[0].append(j*gap5) # gap value in the first row

    for i in range(1,row):
        for j in range(1,col):
            three = []
            if (seq2[i-1] == seq1[j-1]):
                three.append(val[i-1][j-1] + match)
            else:
                three.append(val[i-1][j-1] + mismatch) # Match or Mismatch
            three.append(val[i-1][j] + (gap3 if j == lseq1 else gapInt))
            three.append(val[i][j-1] + (gap3 if i == lseq2 else gapInt))
	    val[i].append(min(three))
    #### print value (if you want to see the matrix of scores
    #for i in range(row):
    #    for j in range(col):
    #        print val[i][j], '\t',
    #    print ''
    score = val[lseq2][lseq1]


    #### trace back
    sequ1 = ''
    sequ2 = ''
    i = lseq2
    j = lseq1
    while (i > 0 or j > 0):
        if (i>0 and j>0 and val[i][j] == val[i-1][j-1] + (match if seq2[i-1]==seq1[j-1] else mismatch)):
            sequ1 += seq1[j-1]
            sequ2 += seq2[i-1]
            i -= 1; j -= 1
	elif (i>0 and val[i][j] == val[i-1][j] + gap(j, lseq1, 0)):
            sequ1 += '-'
            sequ2 += seq2[i-1]
            i -= 1
        elif (j>0 and val[i][j] == val[i][j-1] + gap(i, lseq2, 0)):
            sequ1 += seq1[j-1]
            sequ2 += '-'
            j -= 1

    sequ1r = ' '.join([sequ1[j] for j in range(-1, -(len(sequ1)+1), -1)])
    sequ2r = ' '.join([sequ2[j] for j in range(-1, -(len(sequ2)+1), -1)])

    print "Sequence 1: ", sequ1r
    print "Sequence 2: ", sequ2r
    print "Score     : ", score

if __name__ == "__main__":
    global_gap()
