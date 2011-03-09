'''
Created on Mar 5, 2011

@author: mkiyer
'''
import argparse
import logging
import os
import sys
import subprocess
import numpy as np
import pysam

from trackfactory.track import get_refs_from_sam, get_refs_from_bam, \
    get_refs_from_bowtie_index, parse_interval

from trackfactory.io.cwiggle import WiggleReader
from trackfactory.io.bed import parse_bed6
from trackfactory import TrackFactory
from trackfactory.sequencetrack import SequenceTrack
from trackfactory.arraytrack import ArrayTrack
from trackfactory.coveragetrack import CoverageTrack
from trackfactory.intervaltrack import IntervalTrack, get_bed_dtype_fields

ACTION_CREATE = "create"
ACTION_ADD = "add"
ACTION_REMOVE = "remove"
ACTION_VIEW = "view"

DESCRIPTION=("command-line tools for working with "
             "pytrackfactory tracks") 

def check_executable(filename):
    devnullfh = open(os.devnull, 'w')
    try:
        subprocess.call([filename], stdout=devnullfh, stderr=devnullfh)
    except OSError:
        return False
    devnullfh.close()
    return True

def create_trackfactory(parser, options):
    if options.file is None:
        parser.error("no filename specified")
    if os.path.exists(options.file):
        parser.error("file %s exists, cannot overwrite" % (options.file))
    # get references
    if options.refs_type == "chrom_sizes":
        if not os.path.exists(options.refs):
            parser.error("references file %s not found" % (options.refs))        
        refs = [tuple(line.strip().split(None,2)) for line in open(options.refs)]
    elif options.refs_type == "sam":
        if not os.path.exists(options.refs):
            parser.error("references file %s not found" % (options.refs))
        refs = get_refs_from_sam(options.refs)
    elif options.refs_type == "bam":
        if not os.path.exists(options.refs):
            parser.error("references file %s not found" % (options.refs))
        refs = get_refs_from_bam(options.refs)
    elif options.refs_type == "bowtie_index":
        if not check_executable("bowtie-inspect"):
            parser.error("'bowtie-inspect' executable not found")        
        refs = get_refs_from_bowtie_index(options.refs)
    tf = TrackFactory(options.file, 'w', refs=refs) 
    tf.close()
    logging.info("created trackfactory %s using refs from %s (%s)" % 
                 (options.file, options.refs, options.refs_type))

def remove_track(parser, options):
    tf = TrackFactory(options.file, "r+")
    if not tf.has_track(options.name):
        tf.close()
        parser.error("trackfactory '%s' does not contain track '%s'" %
                     (options.file, options.name))
    tf.delete_track(options.name)    
    tf.close()
    logging.info("removed track '%s' from trackfactory '%s'" %
                 (options.name, options.file))

def open_trackfactory(parser, options):
    tf = TrackFactory(options.file, "r+")
    if tf.has_track(options.name):
        tf.close()
        parser.error("trackfactory '%s' already has track named '%s'" %
                     (options.file, options.name))
    return tf

def check_datafile(parser, options):
    if options.data_file is not None:
        if not os.path.isfile(options.data_file):
            parser.error("data file '%s' not found or not a regular file" % 
                         (options.data_file))
    return options.data_file
    
def add_sequence_track(parser, options):
    tf = open_trackfactory(parser, options)
    t = tf.create_track(options.name, SequenceTrack, options.bpb)
    logging.info("added SequenceTrack %s to trackfactory %s" %
                 (options.name, options.file))
    if options.fasta_file is not None:
        logging.info("inserting fasta file %s" % (options.fasta_file))
        t.fromfasta(open(options.fasta_file))
    tf.close()

def add_array_track(parser, options):
    tf = open_trackfactory(parser, options)
    t = tf.create_track(options.name, ArrayTrack, np.dtype(options.dtype))
    logging.info("added %s %s to trackfactory %s" %
                 (ArrayTrack.__name__, options.name, options.file))
    datafile = check_datafile(parser, options)
    if datafile is not None:
        logging.info("inserting data file %s (type=%s)" % 
                     (options.data_file, options.file_type))
        if options.file_type == "wiggle":
            # add wiggle file
            t.fromintervals(WiggleReader(open(options.data_file)))
    tf.close()

def add_coverage_track(parser, options):
    tf = open_trackfactory(parser, options)
    t = tf.create_track(options.name, CoverageTrack)
    logging.info("added %s %s to trackfactory %s" %
                 (CoverageTrack.__name__, options.name, options.file))
    datafile = check_datafile(parser, options)
    if datafile is not None:
        logging.info("inserting data file %s (type=%s)" % 
                     (options.data_file, options.file_type))
        if options.file_type == "bam":
            # add BAM file
            bamfh = pysam.Samfile(options.data_file, "rb")
            t.frombam(bamfh, 
                      options.norm_rlen,
                      options.bam_nh,
                      options.bam_prob,
                      options.max_multihits,
                      options.keep_dup)
            bamfh.close()
    tf.close()

def add_interval_track(parser, options):
    tf = open_trackfactory(parser, options)
    datafile = check_datafile(parser, options)
    if options.file_type == "bed":
        dtype = np.dtype(get_bed_dtype_fields())
        t = tf.create_track(options.name, IntervalTrack, 
                            dtype=dtype)
        if datafile is not None:
            logging.info("inserting data file %s (type=%s)" % 
                         (options.data_file, options.file_type))
            t.fromintervals(parse_bed6(open(options.data_file)),
                            index=True)
    logging.info("added %s %s to trackfactory %s" %
                 (IntervalTrack.__name__, options.name, options.file))
    tf.close()

def view_track(parser, options):
    tf = TrackFactory(options.file, "r")
    if not tf.has_track(options.name):
        tf.close()
        parser.error("trackfactory '%s' does not contain track '%s'" %
                     (options.file, options.name))    
    region = parse_interval(options.region)
    t = tf.get_track(options.name)
    track_type = t.get_type()
    logging.debug("opened track '%s' type '%s'" % (options.name, track_type))        
    if track_type == CoverageTrack.__name__:
        if options.file_type == "bedgraph":
            t.tobedgraph(region, sys.stdout)
        else:
            print t[region]
    elif track_type == ArrayTrack.__name__:
        if options.file_type == "bedgraph":
            t.tobedgraph(region, sys.stdout)
        else:
            print t[region]
    logging.debug("done")
    tf.close()

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    subparsers = parser.add_subparsers(title="tracktools commands",
                                       description="valid subcommands",
                                       help="")
    # create the parser for the "new" command
    parser_new = subparsers.add_parser(ACTION_CREATE, 
                                       help='create a new trackfactory')
    parser_new.add_argument('file', help='file to create')
    parser_new.add_argument('refs', help="file with reference info")
    parser_new.add_argument('--chrom-sizes', dest="refs_type",
                            action="store_const", const="chrom_sizes", 
                            help="references file is a tab-delimited file "
                            "with each line containing a reference "
                            "(chromosome) name and size")
    parser_new.add_argument('--sam', dest="refs_type", 
                            action="store_const", const="sam", 
                            help="references derived from SAM header")
    parser_new.add_argument('--bam', dest="refs_type", 
                            action="store_const", const="bam", 
                            help="references derived from BAM header")
    parser_new.add_argument('--bowtie-index', dest="refs_type", 
                            action="store_const", const="bowtie_index", 
                            help="references derived from bowtie index "
                            "via the 'bowtie-inspect' command.  the "
                            "'bowtie-inspect' command must be available in "
                            "the current PATH")
    parser_new.set_defaults(func=create_trackfactory,
                            refs_type="chrom_sizes")
    #
    # create parser for the "remove" command
    #
    parser_rm = subparsers.add_parser(ACTION_REMOVE, 
                                      help="remove (delete) a track")
    parser_rm.add_argument('file', help='trackfactory file')
    parser_rm.add_argument('name', help="name of track to remove")
    parser_rm.set_defaults(func=remove_track)
    #
    # create the parser for the "add" command
    #
    parser_add = subparsers.add_parser(ACTION_ADD, help="add a track to "
                                       "an existing trackfactory")  
    parser_add.add_argument('file', help='trackfactory file')
    parser_add.add_argument('name',
                            help="name of track to add")
    # subparsers for each type of TrackFactory
    addsubparsers = parser_add.add_subparsers(title="track types",
                                              description="valid track types",
                                              help="")
    #
    # add a SequenceTrack
    # 
    parser_seq = addsubparsers.add_parser(SequenceTrack.__name__,
                                          help='sequence track')
    parser_seq.add_argument("bpb", type=int, choices=[2,3,4],
                            default=4,
                            help="Number of bits per base to store "
                            "(2bpb converts N -> A, 3bpb preserves 'N', "
                            "4bpb preserves 'N' and upper/lower case info) "
                            "[default=%(default)s]")    
    parser_seq.add_argument("fasta_file", default=None,
                            metavar="seq.fa",
                            help="(optional) if present, insert this fasta " 
                            "file into the track")
    parser_seq.set_defaults(func=add_sequence_track)    

    #
    # add an ArrayTrack
    # 
    parser_arr = addsubparsers.add_parser(ArrayTrack.__name__,
                                          help='array track')
    parser_arr.add_argument("dtype", default="f",
                            help="numpy dtype string specifying the "
                            "type of data to be stored in the track "
                            "[default=%(default)s]")
    parser_arr.add_argument("--wiggle", dest="file_type",
                            action="store_const", const="wiggle",
                            help="data file is in UCSC wiggle format")
    parser_arr.add_argument("--bed", dest="file_type",
                            action="store_const", const="bed",
                            help="data file is in UCSC BED format")
    parser_arr.add_argument("data_file", default=None,
                            help="(optional) file containing data to "
                            "insert into track")
    parser_arr.set_defaults(func=add_array_track,
                            file_type="wiggle")
    #
    # add a CoverageTrack
    #
    parser_cov = addsubparsers.add_parser(CoverageTrack.__name__,
                                          help='coverage track')
    parser_cov.add_argument("--bam", dest="file_type",
                            action="store_const", const="bam",
                            help="data file is in BAM format")
    parser_cov.add_argument("--norm-rlen", dest="norm_rlen",
                            action="store_true", default=False,
                            help="normalize reads by length")
    parser_cov.add_argument("--max-multihits", dest="max_multihits", 
                            default=None,
                            help="(BAM) maximum number of read hits "
                            "allowed during coverage calculation")
    parser_cov.add_argument("--dup", dest="keepdup", 
                            action="store_const", const=True, 
                            help="(BAM) count duplicate reads")
    parser_cov.add_argument("--nodup", dest="keep_dup", 
                            action="store_const", const=False, 
                            help="(BAM) count duplicate reads")
    parser_cov.add_argument("--bam-nh-tag", dest="bam_nh", default=None,
                            help="(BAM) tag containing number of hits per "
                            "read for use in calculate multimapping coverage")
    parser_cov.add_argument("--bam-prob-tag", dest="bam_prob", default=None,
                            help="(BAM) tag containing alignment "
                            "probability read for use in calculating "
                            "multimapping coverage")
    parser_cov.add_argument("data_file", default=None,
                            help="(optional) file containing data to "
                            "insert into track")
    parser_cov.set_defaults(func=add_coverage_track,
                            file_type="bam",
                            keepdup=True)
    #
    # add an IntervalTrack
    #
    subparser = addsubparsers.add_parser(IntervalTrack.__name__,
                                         help='interval track')
    subparser.add_argument("--bed", dest="file_type",
                           action="store_const", const="bed",
                           help="data file is in UCSC BED format")
    subparser.add_argument("data_file", default=None,
                           help="(optional) file containing data to "
                           "insert into track")    
    subparser.set_defaults(func=add_interval_track,
                           file_type="bed")
    #
    # create parser for "view" command
    #
    parser_view = subparsers.add_parser(ACTION_VIEW, 
                                       help="view data in a track")
    parser_view.add_argument('file', help='trackfactory file')
    parser_view.add_argument('name', help="track name")
    parser_view.add_argument('--bedgraph', action="store_const", 
                             const="bedgraph", dest="file_type", 
                             help="output data in bedgraph format")
    parser_view.add_argument('region', nargs="?", default=None, 
                             help="genomic region")
    parser_view.set_defaults(func=view_track,
                             file_type=None)
    #    
    # parse the args and call whatever function was selected
    #
    args = parser.parse_args()
    args.func(parser, args)


if __name__ == '__main__':
    main()
