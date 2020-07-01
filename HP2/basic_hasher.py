import sys
import argparse
import numpy as np
import time
import zipfile
import random
from itertools import product
from collections import defaultdict
import re

def parse_reads_file(reads_fn):
    """
    :param reads_fn: the file containing all of the reads
    :return: outputs a list of all paired-end reads

    HINT: This might not work well if the number of reads is too large to handle in memory
    """
    try:
        with open(reads_fn, 'r') as rFile:
            print("Parsing Reads")
            first_line = True
            count = 0
            all_reads = []
            for line in rFile:
                count += 1
                if count % 1000 == 0:
                    print(count, " reads done")
                if first_line:
                    first_line = False
                    continue
                ends = line.strip().split(',')
                all_reads.append(ends)
        return all_reads
    except IOError:
        print("Could not read file: ", reads_fn)
        return None


def parse_ref_file(ref_fn):
    """
    :param ref_fn: the file containing the reference genome
    :return: a string containing the reference genome
    """
    try:
        with open(ref_fn, 'r') as gFile:
            print("Parsing Ref")
            first_line = True
            ref_genome = ''
            for line in gFile:
                if first_line:
                    first_line = False
                    continue
                ref_genome += line.strip()
        return ref_genome
    except IOError:
        print("Could not read file: ", ref_fn)
        return None


def splitreads(inreads):
    reads = []
    for i in inreads:
        reads.append(i[0])
        reads.append(i[1])
    return reads

def hashgenome(genome):
    gdict = defaultdict(list)
    for i in range(0,len(genome)-10):
        gdict[genome[i:i+10]].append(i)
    return gdict


def divide(dna, d):
    """ Divide a string into small pieces so that at least one fragment is free of errors. d is the number of errors"""
    length = len(dna)
    l = length // (d + 1)
    k = length % (d + 1)  ## higher d will split patterns into shorter length
    result = []
    i = 0
    while i < length:  ## not very important just doing some splitting
        if k > 0:
            result.append((dna[i:i + l + 1], i))
            k -= 1
            i += l + 1
        else:
            result.append((dna[i:i + l], i))
            i += l
    return result


def HammingDistance(s1, s2):  ##! can use anything else, but Hamming is very fast and easy
    if len(s1) != len(s2): return float('inf')
    l = []
    for i in range(len(s1)):
        if s1[i] != s2[i]:
            l.append(i)
    return l

# match read to genome if possible
def match(pattern, gdict):
    return gdict[pattern] if pattern in gdict else -1

def APM(pattern, d, Text, gdict,indellist):
    """ Approximate Pattern Matching """
    result = set() # will only keep unique elements
    l = len(pattern)
    fragments = divide(pattern, d)
    matchcount = 0
    for f in fragments:  ## each fragment of each pattern
        s, i = f
        # print(s)
        rest = pattern[:i] + pattern[i + len(s):]  ## get the rest of pattern (note fragment came from pattern)... may be we can speed this up...
        pat = match(s,gdict) # must be perfect match
        if pat!=-1: # check the read surrounding the perfect match to see if there's at most than d errors
             # number of positions that meet error threshold
            for p in pat:  ##! use BWT to look where this small fragment of pattern match. note that B.lookup is just doing FM indexing, change it into your own function
                ##! note: @p should be an interger saying where the match is found, can match at many places
                # check if any position p has the rest of the read matching w <=d mismatches
                target = Text[p - i:p] + Text[p + len(s):p - i + l]  ## get the text (genome position) that matches rest.
                h = HammingDistance(rest, target)
                if type(h)!=float and len(h) <= d:  ## keep results if we have only "a few" errors.
                    result.add(p - i)
                    matchcount+=1
                else:
                    indellist.append(p-i)

    if matchcount>0: # at least one pattern matches, then return the result set
        # print(str(result),'matches')
        return result
    else: # no pattern matches, then the highlist is indels
        # print(str(indellist), 'are unmatched, adding to indels' )
        return set(indellist)



def solve(patterns, d, dna, gdict, SNPlist, indel_rdict):  ##! assume some input: long-string, patterns-to-match, some-error-rate
    result = defaultdict(list)

    for p in patterns:
        print(p)
        indellist = []
        a = APM(p, d, dna,gdict,indellist)
        if a!=set(indellist): # add positions to reads dictionary
            for i in a:
                result[p].append(i)
                # print(i,'added')
                # add snps
                h = HammingDistance(p,dna[i:i+50])
                if type(h)!= float:
                    for sn in h:
                        # print(sn)
                        SNPlist.add(sn+i)
                # SNPdict[p].append()
        elif a==set(indellist):
            for i in a: # add positions to indel dictionary
                indel_rdict[p].append(i)
                # print("indel dict added ", i)
        print("-"*50)
    # return ' '.join(sorted(map(str, result)))
    return result



def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    return [ atoi(c) for c in re.split(r'(\d+)', text) ]




# make SNP dictionary to find true snps
def SNP(snplist,reference,rdict):
    """
    :param snplist: list of snp positions
    :param reference: reference genome
    :param rdict:  directory of reads and their position + snps
    :return: snpdict: ref snp and position, and snps of other reads
    """
    SNPdict = defaultdict(list)
    for read in rdict:
        # print(read)
        # print(rdict[read])
        r = rdict[read][0]
        # print("create all indices: ", r, r + 50)
        # indices = [i for i in range(r[0],r[0]+50)]
        # print(indices)
        for d in range(r, r + 50):
            if d in snplist:
                ref = reference[d]
                donor = read[d - r]
                SNPdict[str(d) + "|" + str(ref)].append(donor)
        # else:
        #     print("misalignment > 0")
    return SNPdict

# find actual snps
def countSNP(snpdict, coverage):
    finalSNPs = []
    for snp in snpdict:
        if len(snpdict[snp])>coverage/6:
            refsnp = snp.split('|')[1]
            pos = snp.split('|')[0]
            A = ['A',0]
            G = ['G',0]
            T = ['T',0]
            C = ['C',0]
            maxsnp = ''
            maxcount=0
            for i in snpdict[snp]:
                if i=='A':
                    A[1]+=1
                if i=='G':
                    G[1]+=1
                if i=='T':
                    T[1]+=1
                if i=='C':
                    C[1]+=1
            for a in [A,G,T,C]:
                if a[1]>maxcount:
                    maxcount = a[1]
                    maxsnp = a[0]
            donsnp = maxsnp
            if refsnp!=donsnp:
                if int(pos)>=0:
                    finalSNPs.append([refsnp,donsnp,int(pos)])
    return finalSNPs

""" Indel Functions"""

def localalignment(ref, read):
    i = 0 # read pointer
    j = 0 # ref pointer
    m = -1 # mismatch penalty
    e = -2 # indel penalty
    o = -10 # starting indel penalty
    score = [[0 for i in range(len(read)+1)] for j in range(len(ref)+1)] # matrix of scores
    backtrack = [["0" for i in range(len(read)+1)] for j in range(len(ref)+1)] # backtracking matrix
    # initialize first row of scores
    score[0][0] = 0
    backtrack[0][0] = "source"
    # store the maximum position in matrix, for free ride to sink
    curr_max_row = 0
    curr_max_col = 0
    for j in range(1,len(score[0])):
        score[0][j] = score[0][j-1]+e+o if j==1 else score[0][j-1]+e
        backtrack[0][j] = "left"
    # find rest of scores using dp
    for i in range(1,len(score)):
        for j in range(len(score[0])):
            diag = score[i - 1][j - 1] + 1 if read[i - 1] == ref[j - 1] else score[i - 1][j - 1] + m
            right = score[i][j-1]+o+e if backtrack[i][j-1]!="left" else score[i-1][j]+e
            down = score[i-1][j]+o+e if backtrack[i-1][j]!="up" else score[i-1][j]+e
            if j==0:
                score[i][j] = down
            else:
                if i==len(read) and j==len(ref):
                    score[i][j] = max(right,down,diag,0,score[curr_max_row][curr_max_col])
                else:
                    score[i][j] = max(right,down,diag,0)
            # assign backtrack
            if score[i][j]==right:
                backtrack[i][j] = "left"
            elif score[i][j]==down:
                backtrack[i][j] = "up"
            elif score[i][j]==diag:
                backtrack[i][j] = "diag"
            elif score[i][j]==0:
                backtrack[i][j] = "source"
            else: #score[i][j]==score[curr_max_row][curr_max_col]
                backtrack[i][j] = "%s,%s" %(str(curr_max_row),str(curr_max_col))
            if i!=len(read) and j!= len(ref) and score[i][j]>score[curr_max_row][curr_max_col]:
                curr_max_row = i
                curr_max_col = j
    path = get_path(backtrack,len(read),len(ref), read, ref)
    return path

# backtracking path
def get_path(backtrack, n, m, read, ref):
    new_ref = ""
    new_read = ""
    i = n
    j = m
    while i!=0 and j!=0:
        # print(backtrack[i][j])
        if backtrack[i][j]=="left":
            j-=1
            new_ref = ref[j]+new_ref
            new_read = '-'+new_read
        elif backtrack[i][j]=="up":
            i-=1
            new_ref = '-'+new_ref
            new_read = read[i]+new_read
        elif backtrack[i][j]=="diag":
            i-=1
            j-=1
            new_read = read[i]+new_read
            new_ref = ref[j]+new_ref
        elif backtrack[i][j]=="source":
            new_read = str(i) + new_read
            new_ref = str(j) + new_ref
            i=0
            j=0
        else:
            s = backtrack[i][j].split(',')
            new_read = str(int(s[0])-1) + new_read
            new_ref = str(int(s[1])-1) + new_ref
            i = int(s[0])
            j = int(s[1])
    return [new_ref, new_read]

# find indels for all reads in indel dictionary
def find_indels(indel_rdict,reference,insertions,deletions):
    indels = defaultdict(list)
    count = 0
    for read in indel_rdict:
        for r in indel_rdict[read]:
            if r+50 <=len(reference):
                ref = reference[r:r+50]
                la = localalignment(ref,read)
                if "-" in la[0]:
                    num = re.search("^\d*",la[0]).group()
                    num2 = re.search("^\d*",la[1]).group()
                    if num=='':
                        num='0'
                    dash = re.findall("-+",la[0][len(num):])
                    for i in dash:
                        localindex = la[0][len(num):].index(i)
                        ind = r + int(num) + la[0][len(num):].index(i)
                        preseq = la[1][len(num2):]
                        seq = preseq[localindex:localindex+len(i)]
                        if "%s,%s" %(seq,str(ind)) not in insertions:
                            insertions.append("%s,%s" %(seq,str(ind)))
                if "-" in la[1]:
                    num = re.search("^\d*", la[0]).group()
                    num2 = re.search("^\d*", la[1]).group()
                    if num2=='':
                        num2='0'
                    dash = re.findall("-+", la[1][len(num2):])
                    for i in dash:
                        localindex = la[1][len(num2):].index(i)
                        ind = r + int(num2) + la[1][len(num2):].index(i)
                        preseq = la[0][len(num):]
                        seq = preseq[localindex:localindex+len(i)]
                        if "%s,%s" %(seq,str(ind)) not in deletions:
                            deletions.append("%s,%s" %(seq,str(ind)))
            # else:
            # print(r, "index out of range")
        # print(read, "analyzed")
        count+=1
        if count%1000 == 0:
            print(count, "reads analyzed")
    return indels


"""
    TODO: Use this space to implement any additional functions you might need

"""

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='basic_hasher.py takes in data for homework assignment 2 consisting '
                                     'of a genome and a set of reads and aligns the reads to the reference genome, '
                                     'then calls SNPS and indels based on this alignment.')
    parser.add_argument('-g', '--referenceGenome', required=True, dest='reference_file',
                        help='File containing a reference genome.')
    parser.add_argument('-r', '--reads', required=True, dest='reads_file',
                        help='File containg sequencing reads.')
    parser.add_argument('-o', '--outputFile', required=True, dest='output_file',
                        help='Output file name.')
    parser.add_argument('-t', '--outputHeader', required=True, dest='output_header',
                        help='String that needs to be outputted on the first line of the output file so that the\n'
                             'online submission system recognizes which leaderboard this file should be submitted to.\n'
                             'This HAS to be one of the following:\n'
                             '1) practice_W_3_chr_1 for 10K length genome practice data\n'
                             '2) practice_E_1_chr_1 for 1 million length genome practice data\n'
                             '3) hw2undergrad_E_2_chr_1 for project 2 undergrad for-credit data\n'
                             '4) hw2grad_M_1_chr_1 for project 2 grad for-credit data\n')
    args = parser.parse_args()
    reference_fn = args.reference_file
    reads_fn = args.reads_file

    # reads_fn = "reads_hw2undergrad_E_2_chr_1.txt"
    input_reads = parse_reads_file(reads_fn)

    if input_reads is None:
        sys.exit(1)

    # reference_fn = "ref_hw2undergrad_E_2_chr_1.txt"
    reference = parse_ref_file(reference_fn)

    if reference is None:
        sys.exit(1)

    # split reads and hash genome
    reads = splitreads(input_reads)  # each read len 50
    print("reads split")
    gdict = hashgenome(reference)
    print("genome hashed")
    # gfile = open("gdict.txt","w")
    # for g in gdict.items():
    #     gfile.write(g)
    # gfile.close()
    # make 5 fragments
    d = 4  # at most 4 errors in a read

    # dna = reference + "$"
    # patterns = reads


    # make empty SNP dictionary
    # make list of snps to use for later
    snplist = set()
    indel_rdict = defaultdict(list)
    print(time.asctime(time.localtime(time.time())))
    rdict = solve(reads,d,reference, gdict, snplist, indel_rdict)
    print(time.asctime(time.localtime(time.time())))
    print("indel_rdict: ", len(indel_rdict))
    print("rdict: ", len(rdict))
    # rdictfile = open("rdict.txt","w")
    # for r,pos in rdict.items():
    #     rdictfile.write("%s|%s" %(r,",".join(pos)))
    # rdictfile.close()
    # make snp dictionary
    snpdict = SNP(snplist,reference,rdict)
    # find actual snps
    snps = countSNP(snpdict,30)

    """ Finding Indels """
    ins = []
    dels = []
    print(time.asctime(time.localtime(time.time())))
    indels = find_indels(indel_rdict,reference,ins,dels)
    print(time.asctime(time.localtime(time.time())))

    insertions = []
    deletions = []
    for i in ins:
        insertions.append([i.split(",")[0],int(i.split(",")[1])])
    for d in dels:
        deletions.append([d.split(",")[0],int(d.split(",")[1])])

    """
        TODO: Call functions to do the actual read alignment here

    """
    # snps = [['A', 'G', 3425]]
    # insertions = [['ACGTA', 12434]]
    # deletions = [['CACGG', 12]]

    output_fn = args.output_file
    zip_fn = output_fn + '.zip'
    with open(output_fn, 'w') as output_file:
        output_file.write('>' + args.output_header + '\n>SNP\n')
        for x in snps:
            output_file.write(','.join([str(u) for u in x]) + '\n')
        output_file.write('>INS\n')
        for x in insertions:
            output_file.write(','.join([str(u) for u in x]) + '\n')
        output_file.write('>DEL\n')
        for x in deletions:
            output_file.write(','.join([str(u) for u in x]) + '\n')
    with zipfile.ZipFile(zip_fn, 'w') as myzip:
        myzip.write(output_fn)
