import sys
import argparse
import time
import zipfile
from collections import defaultdict


def parse_reads_file(reads_fn):
    """
    :param reads_fn: the file containing all of the reads
    :return: outputs a list of all paired-end reads
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


"""
    TODO: Use this space to implement any additional functions you might need

"""

# get key given a value
def get_key(val,gdict):
    for key, value in gdict.items():
        for i in value:
            if val == i:
                return key
    return "key doesn't exist"


# find number of matches btwn two strings
def mismatches(refstring,samplestring):
    mm = []
    index = 0
    for a,b in zip(refstring,samplestring):
        if a!=b:
            mm.append(index)
        index+=1
    return mm

# split paired end reads
def splitreads(inreads):
    reads = []
    for i in inreads:
        reads.append(i[0])
        reads.append(i[1])
    return reads


# make dictionary with genome: key = 10bp, sequence val =  positions list
def hashgenome(genome):
    gdict = defaultdict(list)
    for i in range(0,len(genome)-10):
        gdict[genome[i:i+10]].append(i)
    return gdict


# find read position on reference genome
def map(reads,gdict):
    readmap = defaultdict(list)
    e = 4
    for r in reads:
        ID = []
        nonnull = []
        frags = []
        null = []
        for i in range(0,50,10):
            fr = r[i:i+10]
            frags.append(fr)
            if fr in gdict:
                ID.append(gdict[fr])
                nonnull.append(int(i/10))
            else:
                ID.append("null")
                null.append(int(i/10))
        if len(nonnull)>1: # reads must have at least 2 reads aligned
            readmap[r].append(ID)
    return readmap


# discover the null sequence in the reads mapping
# fix any mistaken alignments
def fixrdict(rdict,gdict,SNPlist):
    fixeddict = defaultdict(list)
    for read in rdict:
        finalID = ['n', 'n', 'n', 'n', 'n']
        error_num = 0
        # find single-position entries in read and put into finalID
        c = 0 # counter
        singles = [] # read indices for singular position
        mults = [] # read indices for mulitple position
        nulls = [] # read indices for null
        ID = rdict[read][0]
        for p in ID:
            if p=='null':
                nulls.append(c)
            elif len(p)==1:
                finalID[c] = p[0]
                singles.append(c)
            else:
                mults.append(c)
            c+=1
        print("single position entries added: ", finalID)
        #
        # test validity of single values
        print('singles val', singles)
        if len(singles)>0: # testing only if there are singles
            rangeval = [i*10 for i in singles]
            print('range val', rangeval)
            poscheckval = [finalID[i]-(i*10) for i in singles] # should all be the same number in the list
            print('poscheckval', poscheckval)
            if sum(poscheckval)/len(singles):
                print('singles values match')
            else:
                print('singles values do not match')
        # resolve the multiples
        for m in mults:
            print("resolving mults")
            found = "False"
            for s in singles:
                diff = abs(s-m)*10
                if finalID[s]+diff in ID[m]:
                    finalID[m]=finalID[s]+diff
                    print("mult matched", finalID[m])
                    found="True"
                    break
                elif finalID[s]-diff in ID[m]:
                    finalID[m]=finalID[s]-diff
                    print("mult matched", finalID[m])
                    found = "True"
                    break
                else:
                    print(finalID[s],' does not connect to ', ID[m])
        # resolve the nulls
        possible_snps = []
        for n in nulls:
            print("resolving nulls")
            counter = 1
            for s in singles:
                diff = abs(s-n)*10
                if n<s: # get null's index and count mismatches
                    mis = mismatches(get_key(finalID[s]-diff,gdict),read[n*10:(n+1)*10])
                    finalID[n] = finalID[s]-diff
                    if counter==1:
                        error_num+=len(mis)
                        counter+=1
                        for a in mis:  # add SNP positions to snp list
                            possible_snps.append(finalID[s]-diff+a)
                        finalID[n] = finalID[s] - diff
                        break
                elif n>s:
                    mis = mismatches(get_key(finalID[s] + diff, gdict), read[n * 10:(n + 1) * 10])
                    if counter == 1:
                        error_num += len(mis)
                        counter+=1
                    for a in mis: # add SNP positions to snp list
                        possible_snps.append(finalID[s]+ diff+a)
                    finalID[n] = finalID[s] + diff
                    break
        if error_num<=4:
            # appending read
            if finalID != ['n', 'n', 'n', 'n', 'n']:
                fixeddict[read].append(finalID)
                for pos in possible_snps:
                    if pos not in SNPlist:
                        SNPlist.append(pos)
                print("FinalID: ", finalID)
                print("SNPs: ", possible_snps)
        else:
            print(finalID, " has ", error_num, " errors")
    return fixeddict


# map SNP dictionary
def SNP(SNPlist,rdict,reference):
    SNPdict = defaultdict(list)
    for read in rdict:
        print(read)
        print(rdict[read])
        mult = len(rdict[read])
        r = rdict[read][0]
        i = 0
        misaligned = 0
        # while i<len(r)-1:
        #     if r[i+1]-10!=r[i]:
        #         print("incorrect alignment")
        #         misaligned+=1
        #     i+=1
        # if misaligned==0:
        print("create all indices: ", r[0], r[0]+50)
        #indices = [i for i in range(r[0],r[0]+50)]
        #print(indices)
        for d in range(r[0],r[0]+50):
            if d in SNPlist:
                ref = reference[d]
                donor = read[d-r[0]]
                SNPdict[str(d)+"|"+str(ref)].append(donor)
        # else:
        #     print("misalignment > 0")
    return SNPdict



def countSNP(snpdict):
    finalSNPs = []
    for snp in snpdict:
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




if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='basic_aligner.py takes in data for homework assignment 1 consisting '
                                     'of a genome and a set of reads and aligns the reads to the reference genome, '
                                     'then calls SNPs based on this alignment')
    parser.add_argument('-g', '--referenceGenome', required=True, dest='reference_file',
                        help='File containing a reference genome.')
    parser.add_argument('-r', '--reads', required=True, dest='reads_file',
                        help='File containg sequencing reads.')
    parser.add_argument('-o', '--outputFile', required=True, dest='output_file',
                        help='Output file name.')
    parser.add_argument('-t', '--outputHeader', required=True, dest='output_header',
                        help='String that needs to be outputted on the first line of the output file so that the '
                             'online submission system recognizes which leaderboard this file should be submitted to.'
                             'This HAS to be practice_W_1_chr_1 for the practice data and hw1_W_2_chr_1 for the '
                             'for-credit assignment!')
    args = parser.parse_args()
    reference_fn = args.reference_file
    reads_fn = args.reads_file

    # reference_fn = "ref_hw1_W_2_chr_1.txt"
    # reads_fn = "reads_hw1_W_2_chr_1.txt"
    reads_fn = "reads_practice_W_3_chr_1.txt"
    reference_fn = "ref_practice_W_3_chr_1.txt"
    input_reads = parse_reads_file(reads_fn)
    if input_reads is None:
        sys.exit(1)
    reference = parse_ref_file(reference_fn)
    if reference is None:
        sys.exit(1)

    """
        TODO: Call functions to do the actual read alignment here
        
    """
    # split reads and make hashtable of all kmers in genome
    reads = splitreads(input_reads)
    print("reads split")
    gdict = hashgenome(reference)
    print("genome hashed")



    # get positions of read fragments
    # fragments are L=10, 5 fragments per read
    rdict = map(reads,gdict)
    "print rdict defined"

    # make SNPlist
    SNPlist = []
    # revise rdict to get rid of nulls and misaligned reads
    revised_rdict = fixrdict(rdict,gdict,SNPlist)
    # find SNP dict and all SNPS from the reads
    SNPdict = SNP(SNPlist,revised_rdict,reference)

    #snps = [['A', 'G', 3425]]

    # get SNPS that differ from reference genome
    snps = countSNP(SNPdict)

    output_fn = args.output_file
    zip_fn = output_fn + '.zip'
    with open(output_fn, 'w') as output_file:
        header = '>' + args.output_header + '\n>SNP\n'
        output_file.write(header)
        for x in snps:
            line = ','.join([str(u) for u in x]) + '\n'
            output_file.write(line)

        tails = ('>' + x for x in ('STR', 'CNV', 'ALU', 'INV', 'INS', 'DEL'))
        output_file.write('\n'.join(tails))

    with zipfile.ZipFile(zip_fn, 'w') as myzip:
        myzip.write(output_fn)
