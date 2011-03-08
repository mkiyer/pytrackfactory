'''
Created on Mar 8, 2011

@author: mkiyer
'''

#
#def profile_bam(bamfh):
#    for read in bamfh:
#        hits += 1
#        if read.mate_is_unmapped:
#            mate_unmapped_reads += 1
#        if read.is_unmapped or read.is_qcfail:
#            unmapped_reads += 1
#            num_read_hits = 0
#            if read.is_qcfail:
#                qcfail_reads += 1
#        else:            
#            if multimap:
#                num_read_hits = read.opt('NH')
#            else:
#                num_read_hits = 1
#            weighted_cov = 1.0 / num_read_hits
#            # count reads
#            reads += weighted_cov
#            if read.is_duplicate:
#                duplicate_reads += weighted_cov
#            read_lengths[read.rlen] += weighted_cov
#            # mismatch counts along read
#            mm_iter = mismatch_re.finditer(read.opt('MD'))
#            offset = 0
#            for mismatch in mm_iter:
#                skip, mmbase = mismatch.groups()
#                offset += int(skip)
#                if mmbase:
#                    pos = (read.rlen - offset - 1) if read.is_reverse else 0
#                    cycle_mismatch_counts[pos] += weighted_cov
#            # total mismatches per read
#            read_mismatch_counts[read.opt('NM')] += weighted_cov
#        multihit_counts[num_read_hits] += 1


def calc_cov(pileupcolumn, norm_rlen, nh_tag, prob_tag, mmap_max):
    cov = 0
    for pileupread in pileupcolumn.pileups:
        alignment = pileupread.alignment
        readcov = 1.0
        if prob_tag is not None:
            readcov *= alignment.opt(prob_tag)
        if nh_tag is not None:
            nh = alignment.opt(nh_tag)
            if nh > mmap_max:
                continue
            readcov /= nh
        if norm_rlen:
            # TODO: should actually use the genomic footprint of
            # the read rather than the read length, because indels
            # and other events can lead to extra coverage added
            readcov /= alignment.rlen
        cov += readcov
    return cov

def bam_pileup_to_array(pileup_iterator, arr,
                        norm_rlen=False,
                        num_hits_tag=None,
                        hit_prob_tag=None,
                        max_multimaps=None):
    """stores read pileups in a chromosome-sized array
    
    :param norm_rlen: whether to normalize pileup coverage by read length
    :param num_hits_tag: samtools tag that contains number of hits per read
    :param hit_prob_tag: (optional)
    """
    dtype = arr.dtype
    chunk_arr = array.array(dtype.char)
    start = None
    end = None
    for pileupcolumn in pileup_iterator:
        pos = pileupcolumn.pos        
        if start == None:
            start, end = pos, pos
        elif end != pos:
            arr[start:end] = np.array(chunk_arr, dtype=dtype)
            del chunk_arr
            chunk_arr = array.array(dtype.char)
            start, end = pos, pos
        end += 1
        cov = calc_cov(pileupcolumn, norm_rlen, num_hits_tag, hit_prob_tag, 
                       max_multimaps)
        chunk_arr.append(cov)
    # flush remaining array elements
    if len(chunk_arr) > 0:
        arr[start:(start + len(chunk_arr))] = np.array(chunk_arr, dtype=dtype)

def calc_cov_intervals(read_iterator,
                       norm_rlen=False,
                       num_hits_tag=None,
                       hit_prob_tag=None,
                       max_multimaps=None):
    """
    :param norm_rlen: whether to normalize pileup coverage by read length
    :param num_hits_tag: samtools tag that contains number of hits per read
    :param hit_prob_tag: (optional)
    """
    for read in read_iterator:
        strand = "-" if read.is_reverse else "+"
        # compute read coverage
        cov = 1.0
        if hit_prob_tag is not None:
            cov *= read.opt(hit_prob_tag)
        elif num_hits_tag is not None:
            nh = read.opt(num_hits_tag)
            if (max_multimaps is not None) and (nh > max_multimaps):
                continue
            cov /= nh
        # find genomic intervals of read alignment
        for start,end,seq in get_genomic_intervals(read):        
            if norm_rlen:
                cov /= (end - start)
            yield start, end, strand, cov, seq
            
#    def frombam(self, bamfh, 
#                norm_rlen=False,
#                num_hits_tag=None,
#                hit_prob_tag=None,
#                max_multimaps=None):
#        for rname in self.get_rnames():
#            if rname not in bamfh.references:
#                logging.debug("Reference %s not found in BAM header" % (rname))
#                continue
#            logging.debug("Adding ref '%s'" % (rname))            
#            arr = self._get_array(rname)
#            # insert the reads from the BAM file into the array
#            bam_pileup_to_array(bamfh.pileup(rname), arr,
#                                norm_rlen=norm_rlen,
#                                num_hits_tag=num_hits_tag,
#                                hit_prob_tag=hit_prob_tag,
#                                max_multimaps=max_multimaps)
#    def insert_sam(self, filename, genome=False, fragment_length=-1, rmdup=False,
#                   multimap=False):
#        bamfh = pysam.Samfile(filename, "rb")
#        # count reads and get other info from BAM file
#        logging.debug("%s: Computing track statistics" % (filename))
#        hits = 0
#        reads = 0.0
#        duplicate_reads = 0.0
#        unmapped_reads = 0.0
#        mate_unmapped_reads = 0.0
#        qcfail_reads = 0.0
#        read_lengths = collections.defaultdict(lambda: 0.0)
#        multihit_counts = collections.defaultdict(lambda: 0.0)
#        cycle_mismatch_counts = collections.defaultdict(lambda: 0.0)
#        read_mismatch_counts = collections.defaultdict(lambda: 0.0)
#        mismatch_re = re.compile(r"(\d+)([agtcAGTC]*)")
#
#        for read in bamfh:
#            hits += 1
#            if read.mate_is_unmapped:
#                mate_unmapped_reads += 1
#            if read.is_unmapped or read.is_qcfail:
#                unmapped_reads += 1
#                num_read_hits = 0
#                if read.is_qcfail:
#                    qcfail_reads += 1
#            else:            
#                if multimap:
#                    num_read_hits = read.opt('NH')
#                else:
#                    num_read_hits = 1
#                weighted_cov = 1.0 / num_read_hits
#                # count reads
#                reads += weighted_cov
#                if read.is_duplicate:
#                    duplicate_reads += weighted_cov
#                read_lengths[read.rlen] += weighted_cov
#                # mismatch counts along read
#                mm_iter = mismatch_re.finditer(read.opt('MD'))
#                offset = 0
#                for mismatch in mm_iter:
#                    skip, mmbase = mismatch.groups()
#                    offset += int(skip)
#                    if mmbase:
#                        pos = (read.rlen - offset - 1) if read.is_reverse else 0
#                        cycle_mismatch_counts[pos] += weighted_cov
#                # total mismatches per read
#                read_mismatch_counts[read.opt('NM')] += weighted_cov
#            multihit_counts[num_read_hits] += 1
#
#        # number of reads
#        parent_group = self.hdf_group
#        parent_group._v_attrs.hits = hits
#        logging.info("%s: Number of alignment hits=%d" % (filename,hits))
#        parent_group._v_attrs.reads = reads
#        logging.info("%s: Number of reads=%d" % (filename,int(reads)))
#        parent_group._v_attrs.duplicate_reads = duplicate_reads
#        logging.info("%s: Number of qc fail reads=%d" % (filename,int(qcfail_reads)))
#        parent_group._v_attrs.qcfail_reads = qcfail_reads           
#        logging.info("%s: Number of read flagged as duplicates=%d" % (filename,int(duplicate_reads)))
#        parent_group._v_attrs.unmapped_reads = unmapped_reads
#        logging.info("%s: Number of mapped reads=%d" % (filename,int(reads-unmapped_reads)))
#        logging.info("%s: Number of unmapped reads=%d" % (filename,unmapped_reads))
#        parent_group._v_attrs.mate_unmapped_reads = mate_unmapped_reads
#        logging.info("%s: Number of reads with unmapped mate=%d" % (filename,mate_unmapped_reads))
#        parent_group._v_attrs.read_length_variable = True if len(read_lengths) > 1 else False
#        # set to the most abundant read length in the case of variable read length runs
#        parent_group._v_attrs.read_length = int(sorted(read_lengths.items(), key=operator.itemgetter(1), reverse=True)[0][0])
#        logging.info("%s: Most common read length=%d" % (filename, parent_group._v_attrs.read_length))
#        parent_group._v_attrs.read_length_distribution = dict(read_lengths)
#        logging.info("%s: Number of different read lengths=%d" % (filename, len(read_lengths)))
#        parent_group._v_attrs.fragment_length = fragment_length
#        logging.info("%s: Fragment length (specified as input)=%d" % (filename, parent_group._v_attrs.fragment_length))
#
#        max_multihits = max(multihit_counts.keys())
#        arr = np.zeros(max_multihits, dtype=np.int)
#        for i in xrange(max_multihits):
#            if i in multihit_counts:
#                arr[i] = int(multihit_counts[i])
#        parent_group._v_attrs.multihits = list(arr)
#        logging.info("%s: Multimapping read counts=%s" % (filename, parent_group._v_attrs.multihits))
#
#        max_rlen = max(read_lengths.keys())
#        arr = np.zeros(max_rlen, dtype=np.int)
#        for i in xrange(max_rlen):
#            if i in cycle_mismatch_counts:
#                arr[i] = int(cycle_mismatch_counts[i])
#        parent_group._v_attrs.cycle_mismatches = list(arr)
#        logging.info("%s: Mismatches per cycle=%s" % (filename, parent_group._v_attrs.cycle_mismatches))
#
#        max_mismatches = max(read_mismatch_counts.keys())
#        arr = np.zeros(max_mismatches, dtype=np.int)
#        for i in xrange(max_mismatches):
#            if i in read_mismatch_counts:
#                arr[i] = int(read_mismatch_counts[i])
#        parent_group._v_attrs.read_mismatches = list(arr)
#        logging.info("%s: Mismatches per read=%s" % (filename, parent_group._v_attrs.read_mismatches))
#
#        # insert data
#        for ref,length in zip(self._get_references(), self._get_lengths()):
#            if ref not in bamfh.references:
#                logging.debug("Skipping reference %s in BAM file that is not part of the track" % (ref))
#                continue
#            logging.debug("Adding file %s ref %s to database" % (filename, ref))
#            ca = self.hdf_group._f_getChild(ref)
#            # insert the reads from the BAM file into the array
#            if genome:
#                bam_to_array_chunked(bamfh.fetch(ref), ca, fragment_length=fragment_length, rmdup=rmdup)
#            elif multimap:
#                bam_pileup_to_array_multihit(bamfh.pileup(ref), ca, dtype=ca.atom.dtype)
#            else:
#                bam_pileup_to_array(bamfh.pileup(ref), ca, dtype=ca.atom.dtype)
#        bamfh.close()