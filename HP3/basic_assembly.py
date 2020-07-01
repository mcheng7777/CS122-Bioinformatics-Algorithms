from os.path import join
import sys
import time
from collections import defaultdict, Counter
import sys
import os
import zipfile
import argparse
sys.path.insert(0, os.path.abspath(".."))
sys.path.insert(0, os.path.abspath("../.."))


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


"""
    TODO: Use this space to implement any additional functions you might need
    # function to make k spectrum pair nodes

    
    
    # function to make adjacency list 

"""
# find indegree of a node
def inDegree(degdict, v):
	return degdict[v][0]

# find outdegree of a node
def outDegree(degdict, v):
	return degdict[v][1]

# evaluate if node is non-branching
def nonBranchNode(degdict, v):
	return (inDegree(degdict, v) == 1) and (outDegree(degdict, v) == 1)


# split apart paired reads
def splitreads(inreads):
    reads = []
    for i in inreads:
        reads.append(i[0])
        reads.append(i[1])
    return reads

# form a dictionary with keys=reads, value = list of kmer composition of read
def rdict_construction(rlist):
    k=30
    l=50
    rdict = defaultdict(list)
    for i in range(0,len(rlist)):
        j=0
        while j+k <=l:
            kmer = rlist[i][j:j + k]
            # account for input reads that do not have full length read pairs (second pair only)
            kmer = kmer+('-'*(k-len(kmer)))
            rdict[i].append(kmer)
            j+=1
    return rdict

# form a dictionary such that for each read, the value is the [ID of kmer's first occurance, and position along the read]
def kmer_pos_dict(rdict):
    """
    purpose is to only have kmer in debruijn graph once, and to find every occurance of the kmer in all reads
    param rdict - constructed read-kmer dictionary
    param kmer_pos - dictionary for all read IDs that possess the kmer
    return kmer_occ - dictionary of read id's , key = read ID, value = list of [first readID,position in first readID]
        - essentially it's the nodes of the debruijn graph

    """
    kmer_occ = defaultdict(list)
    enc_kmers = defaultdict(list) # dictionary of encountered kmers, key = kmer, val = [read ID,position]
    for ID in rdict:
        print ("searching read ", ID)
        r = rdict[ID]
        # each kmer of the read will be given it's first occurence position
        for i in range(0,len(r)): # need to traverse rdict[ID] by index
            kmer = r[i]
            print("kmer is ", kmer)
            kmer_occ[kmer].append([ID,i])
            # if kmer in enc_kmers: # first occurrence of the kmer already recorded
            #     print(kmer, " already encountered")
            #     kmer_occ[ID].append(enc_kmers[kmer])
            # else: # first occurrence of kmer not yet recorded
            #     print(kmer, " being added")
            #     enc_kmers[kmer] = [ID,i]
            #     kmer_occ[ID].append(enc_kmers[kmer])
            # kmer_pos[kmer].append(ID) # add ID to kmer's total positions dictionary
    return kmer_occ



# count total number of kmers
def kcounts(rdict):
    counts = defaultdict(list)
    for ID in rdict:
        for i in rdict[ID]:
            if i in counts:
                counts[i]+=1
            else:
                counts[i]=1
    return counts



# form degree dictionary key = k-1mer, val = [indegree,outdegree]
# function to make k spectrum pair nodes
def node_construction(edges,degdict):
    """
    :param: edges - list of kmer edges
    :param: degdict - degree dictionary, key = k-mer, val = [indegree, outdegree]
    :return: gdict - k-mer dictionary, key = "from" node,  val = lists of "to" nodes
    """
    gdict = defaultdict(list)  # dict, key="from" node,  val="to" node
    for i in edges:
        pref = i[0:-1]
        suff = i[1:]
        gdict[pref].append(suff)
        if pref in degdict:
            degdict[pref][1]+=1 # outdegree +=1
            # print("prefix added \n")
        else:
            degdict[pref] = [0,1]
            # print("prefix added \n")
        if suff in degdict:
            degdict[suff][0]+=1 # indegree +=1
            # print("suffix added \n")
        else:
            degdict[suff] = [1,0]
            # print("suffix added \n")
    return gdict

# find paths based off of input reads (string reconstruction)
def paths(gdict,kpos,starts,degdict):
    """
    param: gdict - adjacency list
    param: kpos - kmer positions in rdict
    param: starts - list of k-1mers to start path
    param: nodelength - length of k-1mer
    param: degdict - degree dictionary
    return: contigs - list of contigs
    """
    # after edge is visited, use values.pop(0) to take out the end node
    # use del gdict[key] to delete key
    k=30
    paths = []
    for n in starts:
        contig = []
        contig.append(n)
        done = "false"
        i = n
        while done=="false":
            if len(gdict[i])==1 or len(contig)==1: # nonbranching node
                next_node = gdict[i][0]
                contig.append(next_node)
                # update degrees
                degdict[i][1]-=1
                degdict[next_node][0]-=1
                # delete edge and update i
                gdict[i].pop(0)
                i = next_node
                print("next node: ", i)
            elif len(gdict[i])>1: # branching node
                c = 0
                # find last kmer in contig in the kpos dict. Check upcoming node to kpos dict one position over
                next_node = gdict[i][0]
                last = contig[-2][0]+i # cuz i is the last added k-1mer. we want right before i
                print("last contig kmer is ", last) # debug purpose
                while c<len(kpos[last]): # may not work lol
                    ID = kpos[last][c][0] # ID of read containing kmer
                    pos = kpos[last][c][1] # position of rdict kmer is located
                    if pos<50-k:
                        if rdict[ID][pos+1]==(i+next_node[-1]):
                            contig.append(next_node) # adds the next node in path
                            # update degrees
                            degdict[i][1] -= 1
                            degdict[next_node][0] -= 1
                            # delete edge and update i
                            gdict[i].pop(0)
                            if degdict[i][0] == 0 and degdict[i][1]>0: # add to starts if i has no incoming edges but still outgoing edges
                                for newstart in gdict[i]:
                                    starts.append(i)
                            i = next_node
                            c = len(kpos[last]) # end while loop
                        else: c+=1
                    else:
                        c+=1
            else: # end node, path ended
                print('what is contig {}'.format(contig))
                r = contig[0]
                r += ''.join(con[-1] for con in contig[1:])  ## add string together.
                paths.append(r)
                done = "true"
    return paths

# find starts of contigs
def contig_starts(degdict):
    cstarts = []
    counter = 0
    for v in degdict:
        if inDegree(degdict,v)==0:
            cstarts.append(v)
        if counter%500==0:
            print(counter, " start nodes added")
        counter+=1
    return cstarts

# create contigs from the kmers in the position dictionary
def find_contigs(adjlist,starts,degdict):
    """
    param: pos_dict - dictionary of kmers and their corresponding read ID
    return: list of contigs

    """
    result = []
    for start in starts:
        print('-' * 50)
        print('\n\nstart a "head" node {}'.format(start))
        for v in adjlist[start]:
            nextV = v
            path = [start, nextV]  ## take any path
            print('\npath to take {}'.format(path))
            print('is this non-branching path {}'.format(nonBranchNode(degdict, nextV)))
            while nonBranchNode(degdict, nextV):
                print('path is direct, so we continue to walk until we are able to branch')
                nextV = adjlist[nextV][0]
                path.append(nextV)
                print(path)

            ##
            print('what is path after "while" {}'.format(path))
            r = path[0]
            r += ''.join(p[-1] for p in path[1:])  ## add string together.
            result.append(r)
    return result




if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='basic_assembly.py takes in data for homework assignment 3 consisting '
                                                 'of a set of reads and aligns the reads to the reference genome.')
    parser.add_argument('-r', '--reads', required=True, dest='reads_file',
                        help='File containg sequencing reads.')
    parser.add_argument('-o', '--outputFile', required=True, dest='output_file',
                        help='Output file name.')
    parser.add_argument('-t', '--outputHeader', required=True, dest='output_header',
                        help='String that needs to be outputted on the first line of the output file so that the\n'
                             'online submission system recognizes which leaderboard this file should be submitted to.\n'
                             'This HAS to be one of the following:\n'
                             '1) spectrum_A_1_chr_1 for 10K spectrum reads practice data\n'
                             '2) practice_A_2_chr_1 for 10k normal reads practice data\n'
                             '3) hw3all_A_3_chr_1 for project 3 for-credit data\n')
    args = parser.parse_args()
    reads_fn = args.reads_file

    input_reads = parse_reads_file(reads_fn)
    if input_reads is None:
        sys.exit(1)

    """
            TODO: Call functions to do the actual assembly here

    """

    # split paired end reads
    reads = splitreads(input_reads)
    # create a reads, kmer composition dictionary
    rdict = rdict_construction(reads)

    # find positions of each kmer [ID,position]
    kpos_dict = kmer_pos_dict(rdict)

    # store kmer occurances in reads
    counts = kcounts(rdict)

    #max(counts.values())

    #sorted(counts.items(), key= lambda kv: (kv[1], kv[0]))

    # find kmers with over 10 coverage
    kmers = defaultdict(list)
    unused = defaultdict(list)
    for k in counts:
        if counts[k] > 5: # insert value
            kmers[k] = counts[k]
        else:
            unused[k] = counts[k]

    # sorted(kmers.items(), key=lambda kv: (kv[1], kv[0]))

    # find coverage numbers, get feel for coverage
    coverage_bins = defaultdict()
    for v in kmers.values():
        if int(v/10) in coverage_bins:
            coverage_bins[int(v/10)]+=1
        else:
            coverage_bins[int(v / 10)] = 1

    # sorted(coverage_bins.items(), key=lambda kv: (kv[0], kv[1]))

    # coverage around 20
    # kmer's appearance in genome = round (count/20)
    # make kmer edges
    edges = []
    for key,value in kmers.items():
        if value >30:
            for i in range(0,int(value/30)):
                edges.append(key)
        else:
            edges.append(key)

    # make adjlist
    degdict = defaultdict(list) # degrees of k-1mers
    debgraph = node_construction(edges,degdict)

    for key,value in degdict.items():
        if value[0]>1 and value[1]>1:
            print(key,value)

    # find starts of contigs
    starts = contig_starts(degdict)

    # get contigs
    # contigs = find_contigs(debgraph,starts,degdict)
    contigs = paths(debgraph,kpos_dict,starts,degdict)

    for i in contigs:
        print(len(i))


    output_fn = args.output_file
    zip_fn = output_fn + '.zip'
    with open(output_fn, 'w') as output_file:
        output_file.write('>' + args.output_header + '\n')
        output_file.write('>ASSEMBLY\n')
        output_file.write('\n'.join(contigs))
    with zipfile.ZipFile(zip_fn, 'w') as myzip:
        myzip.write(output_fn)
