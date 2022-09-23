"""
Transform a BAM file to a count table based on barcode information and genomic bins
"""

import time
import argparse
import sys
import os
import re
import math
from collections import OrderedDict
import pysam
import numpy as np
from scipy import sparse
import scipy.io as sio
from bx.intervals.intersection import Intersecter, Interval
import subprocess
    
def timing(function, *args):
    """                              
    Run a fonction and return the run time and the result of the function
    If the function requires arguments, those can be passed in
    """
    startTime = time.time()
    result = function(*args)
    print '%s function took %0.3f ms' % (function.func_name, (time.time() - startTime) * 1000)
    return result


def load_BED(in_file, featuresOverCoord=False, verbose=False):
    """
    Read a BED file and store the intervals in a tree
    Intervals are zero-based objects. The output object is a hash table with
    one search tree per chromosome

    BED file are half-open, meaning that a bin ]100, 200] covered the bases 101 to 200
    
    in_file = input file [character]
    verbose = verbose mode [logical]
    """
    x = {}
    if verbose:
        print "## Loading BED file '", in_file, "'..."
    featureNames=[]
    nline = 0
    with open(in_file) as bed_handle:
        for line in bed_handle:
            if nline > 0 and nline % 5000==0 and verbose: 
                print "## %d features loaded ..." % nline
            nline +=1
            bedtab = line.split("\t")
            chromosome, start, end = bedtab[:3]
            if len(bedtab)>3 & featuresOverCoord==True:
                name = bedtab[4]

            # BED files are zero-based, half-open as Intervals objects
            start = int(start) 
            end = int(end)
            if featuresOverCoord==True:
                featureNames.append(name.strip())
            else:
                featureNames.append(chromosome + ":" + str(start) + "-" + str(end))
            
            if chromosome in x:
                tree = x[chromosome]
                tree.add_interval(Interval(start, end, value={'pos' : nline - 1}))
            else:
                tree = Intersecter()
                tree.add_interval(Interval(start, end, value={'pos' : nline - 1}))
                x[chromosome] = tree
    bed_handle.close()
    return (x, featureNames)


def get_barcode_number_from_header(sam):
    """
    Read the BAM header, and extract the CO tag with barcode number that 
    is added during the flag step
    """
    barcode_number = None
    if 'CO' in sam.header:
        COitems = sam.header['CO']
        for com in COitems:
            if re.search("Barcodes", com):        
                barcode_number = int(com.split(":")[1])
    return barcode_number


def get_read_tag(read, tag):
    """
    Extract a flag from a read alignment
    """
    for t in read.tags:
        if t[0] == tag:
            return t[1]
    return None

def get_read_start(read):
    """
    Return the 5' end of the read
    Same as reference_start / reference_end in recent pysam version
    """
    if read.is_reverse:
        pos = read.pos + read.alen
    else:
        pos = read.pos
    return int(pos)


def get_chromosome_size_from_header(sam):
    """
    Extract chromosome size from header. 
    That way, we do not need any chromosome information from the user
    """
    chromSize = OrderedDict()
    SQitems = sam.header['SQ']
    for chrom in SQitems:
        chromSize[chrom['SN']] = chrom['LN']
    return(chromSize)


def get_chromosome_bins(chroms, bsize):
    """
    Get number of genomic bins per chromosome
    Is used to tranform genomic coordinates, into a bin number in 
    the count matrix
    """
    chrombins = OrderedDict()
    for chrname in chroms:
        x = float(chroms[chrname]) / float(bsize)
        chrombins[chrname] = int(math.ceil(x))
    return(chrombins)


def get_features_idx(intervals, chrom, read, useWholeRead=True, verbose=False):
    """
    Intersect a given read with the set of intervals
    intervals = the fragments [hash]
    chrom = the chromosome to look at [character]
    read = the read to intersect [AlignedRead]
    useWholeRead = True/False, use either the 5' end of the read or the full length
    """
    
    if useWholeRead:
        lpos = read.pos
        rpos = read.pos + read.alen
    else :
        if read.is_reverse:
            rpos = get_read_start(read)
            lpos = rpos - 1
        else:
            lpos = get_read_start(read)
            rpos = lpos + 1

    if chrom in intervals:

        # Overlap with the left/right positions of the read (zero-based)
        feat = intervals[chrom].find(lpos, rpos)

        if len(feat) == 0:
            if verbose: print >> sys.stderr, "Warning - no feature found for read at", chrom, ":", read.pos, "- skipped"
            return None
        else:
            feat_idx = []
            for i in range(len(feat)): 
                feat_idx.append(feat[i].value['pos'])
            return feat_idx
    else:
        if verbose: print >> sys.stderr, "Warning - no feature found for read at", chrom, ":", read.pos, "- skipped"
        return None


def get_bin_idx (chrname, read, chrom_idx, binSize, useWholeRead=True):
    """
    For a given chromosome and position, return the bin indice
    Return a zero-based index
    """
    try:
        chrpos = chrom_idx.keys().index(chrname)

    except ValueError:
        print >> sys.stderr, "Chromosome " + chrname + " not found !"
        sys.exit(-1)
    
    if useWholeRead:
        lpos = read.pos
        rpos = read.pos + read.alen
        if read.is_reverse:
            ## Require for reverse reads that start exactly at the end of a bin ...
            rpos = rpos - 1

        ### Sum of previous chromosome + current bin (0-based)
        idx = sum(chrom_idx.values()[:chrpos]) + int(math.floor(float(lpos) / float(binSize)))
        idx_end = sum(chrom_idx.values()[:chrpos]) + int(math.floor(float(rpos) / float(binSize)))
        
        if idx != idx_end:
            return range(idx, idx_end + 1)
        else:
            return [idx]
    else:
        lpos = get_read_start(read)
        if read.is_reverse:
            lpos = lpos - 1
        idx = sum(chrom_idx.values()[:chrpos]) + int(math.floor(float(lpos) / float(binSize)))
        
        return [idx]
    

def get_bins_coordinates(i, chromsize, chrom_idx, bsize):
    """
    Transform a bin indice into a genomic coordinate
    Need indice, chromosome size, number of indices per chromosome and binsize
    """
    cumidx = 0
  
    for k in chrom_idx.keys():
        if cumidx == 0:
            chrname = k
        if cumidx + chrom_idx[k] > i:
            chrname = k
            break
        else:
            cumidx += chrom_idx[k]
    
    start = i * bsize - cumidx * bsize
    end = int(start) + int(bsize)
    if end > chromsize[chrname]:
        end = chromsize[chrname]
    return np.array([str(chrname), str(start), str(end)])


def select_mat(x, nreads=500, verbose=False):
    """
    Select counts matrix columns
    """
    cx = x.tocsr()
    cols = np.array(cx.sum(axis=0))  # sum the columns
    idx = np.where(cols >= nreads)[1]
    if verbose:
        print "## Select " + str(len(idx)) + " columns with at least " + str(nreads) + " counts"
    return idx


def save2BinSparseMatrix(x, colnames, odir, chromsize, chrom_idx, bsize, filt, verbose=False, rmzeros=False):
    """
    Write the count table into a txt file without taking too much RAM
    Note that most of the args are used to convert bin coordinate into genomic coordinates
    """
    if verbose:
        print "## Writting output file ..."

    cx = x.tocsr()
    odir = odir + "_"+ str(bsize)
    if not os.path.exists(odir):
        os.mkdir(odir)

    mtx_file = odir + "/matrix.mtx"
    features_file = odir + "/features.tsv"
    barcodes_file = odir + "/barcodes.tsv"

    # Save matrix.mtx
    handle = open(mtx_file,'ab')
    sio.mmwrite(mtx_file, cx)
    handle.close()
    subprocess.Popen("gzip " + mtx_file, shell=True, stdin=subprocess.PIPE)
    
    # Save barcodes.tsv
    handle = open(barcodes_file,'ab')
    colnames = np.array([colnames])
    np.savetxt(handle, colnames, '%s', delimiter="\n")
    handle.close()
    subprocess.Popen("gzip " + barcodes_file,shell=True, stdin=subprocess.PIPE)

    # Save features.tsv
    handle = open(features_file,'ab')

    for i in range(cx.shape[0]):
        coord = get_bins_coordinates(i, chromsize, chrom_idx, bsize)
	coord = np.array([coord])
        np.savetxt(handle, coord, fmt='%s', delimiter="\t")
    
        if (i > 0 and i % 10000 == 0 and verbose):
            print "## Write line " + str(i)

    handle.close()
    subprocess.Popen("gzip " + features_file, shell=True, stdin=subprocess.PIPE)


def save2FeatSparseMatrix(x, colnames, odir, rownames, filt, verbose=False, rmzeros=False):
    """
    Write the count table into a txt file without taking too much RAM
    """
    if verbose:
        print "## Writting outpout file ..."

    
    cx = x.tocsr()
    if not os.path.exists(odir):
        os.mkdir(odir)

    mtx_file = odir + "/matrix.mtx"
    features_file = odir + "/features.tsv"
    barcodes_file = odir + "/barcodes.tsv"

    # Save matrix.mtx
    handle = open(mtx_file,'ab')
    sio.mmwrite(mtx_file, cx)
    handle.close()
    subprocess.Popen("gzip " + mtx_file, shell=True, stdin=subprocess.PIPE)

    # Save barcodes.tsv
    handle = open(barcodes_file,'ab')
    colnames = np.array([colnames])
    np.savetxt(handle, colnames, '%s', delimiter="\n")
    handle.close()
    subprocess.Popen("gzip " + barcodes_file,shell=True, stdin=subprocess.PIPE)

    # Save features.tsv
    handle = open(features_file,'ab')

    for i in range(cx.shape[0]):
        name = np.array([rownames[i]])
        np.savetxt(handle, name, '%s', delimiter="\t")
    
        if (i > 0 and i % 10000 == 0 and verbose):
            print "## Write line " + str(i)
 
    handle.close()
    subprocess.Popen("gzip " + features_file, shell=True, stdin=subprocess.PIPE)


if __name__ == "__main__":

    # Init variables
    reads_counter = 0
    non_overlapping_counter = 0
    allbarcodes = []
   
    # Reads args
    parser = argparse.ArgumentParser(prog='sc2counts.py', description='''
    Transform a BAM file to a count table based on barcode information and genomic bins/features
    ''', epilog='''
    Note that --bin and --bed options are exclusive ! Counts are generated either per bin (--bin) or per genomics features (--bed)
    ''')
    parser.add_argument('-i', '--input', help="BAM file with barcode tag", required=True)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-b', '--bin', help="Size of genomic bin size", default=None, type=int)
    group.add_argument('-B', '--bed', help="BED file of genomic features", default=None, type=str)
    parser.add_argument('-o', '--output', help="Output directory where to save matrix.mtx, features.tsv and barcodes.tsv files.", default="./sparse_matrix", type=str)
    parser.add_argument('-s', '--barcodes', help="Number of barcodes in the BAM file. Default: Extracted from the BAM 'CO' field", type=int)
    parser.add_argument('-t', '--tag', help="Barcode Tag. Default: XB", default="XB", type=str)
    parser.add_argument('-f', '--filt', help="Select barcodes with at least FILT counts. Default: None", default="1", type=str)
    parser.add_argument('-w', '--useWholeRead', help="Use the whole read in the count instead of the 5' end. Default: False", default=False, action="store_true")
    parser.add_argument('-r', '--rmZeros', help="Do not export bins/features with only zeros in the count table. Default: False", default=False, action="store_true")
    parser.add_argument('-F', '--featuresOverCoord', help="When counting on BED file, write feature name (column 4 of BED) as rownames of count matrix instead of coordinates. Default: False", default=False, action="store_true")
    parser.add_argument('-v', '--verbose', help="", action="store_true")
 
    args = parser.parse_args()
 
    # check args
    if args.bed is None and args.bin is None:
        print "Error: --bin or --bed must be defined"
        sys.exit(-1)
        
    # Verbose mode
    if args.verbose:
        print "## sc2counts.py"
        print "## input =", args.input
        print "## output =", args.output
        print "## tag =", args.tag
        print "## binSize =", args.bin
        print "## bedFile =", args.bed
        print "## barcodeNumber =", args.barcodes
        print "## minCounts =", args.filt
        print "## useWholeReads =", args.useWholeRead
        print "## rmZeros =", args.rmZeros
        print "## verbose =", args.verbose
        print "## featuresOverCoord =", args.featuresOverCoord
        print

    # Read the SAM/BAM file
    if args.verbose:
        print "## Opening SAM/BAM file '", args.input,"'..."

    samfile = pysam.Samfile(args.input, "rb")

    # Get info from header
    chromsize = get_chromosome_size_from_header(samfile)
    if len(chromsize)==0:
        print >> sys.stderr, "Error : chromosome lengths not available in BAM file. Exit"
        sys.exit(-1)

    # Get counts dimension
    if args.barcodes is None:
        N_barcodes = get_barcode_number_from_header(samfile)
        if N_barcodes is None :
            print >> sys.stderr, "Erreur : unable to find barcodes number. Exit"
            sys.exit(-1)
        elif args.verbose:
            print "## Barcodes Number: " + str(N_barcodes)
    elif args.barcodes is not None:
        N_barcodes = args.barcodes
        print "## Barcodes Number: " + str(N_barcodes)	
    if args.bin is not None:
        chromsize_bins = get_chromosome_bins(chromsize, args.bin)
        N_bins = sum(chromsize_bins.values())
    elif args.bed is not None:
        feat_bins = load_BED(args.bed, args.featuresOverCoord, args.verbose)
        N_bins = len(feat_bins[1]) 
 
    if args.verbose:
        print "## Bins/Features Number: " + str(N_bins)

    ## Init count table
    ## Note that lil matrix is a sparse (ie. RAM eficient) matrix
    ## design to speed incrementation
    counts = sparse.lil_matrix((N_bins, N_barcodes))

    for r1 in samfile.fetch(until_eof=True):
        reads_counter += 1

        r1_chrom = samfile.getrname(r1.tid)
  
        ## Get barcode
        barcode = str(get_read_tag(r1, args.tag))
        
        ## Get Barcode (ie col) indices
        try:
            j = allbarcodes.index(barcode)
        except ValueError:
            allbarcodes.append(barcode)
            j = len(allbarcodes)-1

        ## Get Bin (ie. rows) indice and increment bin matrix
        if args.bin is not None:
            i = get_bin_idx(r1_chrom, r1, chromsize_bins, args.bin, useWholeRead=args.useWholeRead)
            for ii in i: 
                counts[ii, j] += 1
        elif args.bed is not None:
            i = get_features_idx(feat_bins[0], r1_chrom, r1, useWholeRead=args.useWholeRead)
            if i is not None:
                for ii in i: 
                    counts[ii, j] += 1
            else:
                non_overlapping_counter += 1

        if (reads_counter % 1000000 == 0 and args.verbose):
            print "##", reads_counter
    
    samfile.close()

    if args.verbose and non_overlapping_counter > 0:
        print "## Warning:", non_overlapping_counter, "reads do not overlap any features !"

    ## filter mat
    if args.filt is not None:
        filters = map(int, args.filt.split(","))
        for filt in filters:
            sel_idx = select_mat(x=counts, nreads=filt, verbose=args.verbose)
            counts_reduced = counts[:, sel_idx]
            allbarcodes_reduced = np.array(allbarcodes)[sel_idx]
            allbarcodes_reduced = allbarcodes_reduced.tolist()

            ## save Matrix
            if args.bin is not None:
                save2BinSparseMatrix(counts_reduced, allbarcodes_reduced, args.output, chromsize, chromsize_bins, args.bin, filt, args.verbose, rmzeros=args.rmZeros)
            elif args.bed is not None:
                save2FeatSparseMatrix(counts_reduced, allbarcodes_reduced, args.output, feat_bins[1], filt ,args.verbose, rmzeros=args.rmZeros)



