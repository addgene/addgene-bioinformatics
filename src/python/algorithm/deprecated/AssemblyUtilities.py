def count_k_mers_in_seq(seq, seq_id=0, k_mer_len=25, k_mers=None):
    """Count k-mers of a specified length in a sequence with a
    specified identifier. Optionally update dictionary returned by
    this method.

    Parameters
    ----------
    seq : Bio.Seq.Seq
        The sequence in which to count
    seq_id : int
        The source sequence identifer
    k_mer_len : int
        The length of k-mers to count
    k_mers : dct
        Emtpy, or returned by this method

    Returns
    -------
    k_mers : dct
        Dictionay containing k-mer keys and counts and source sequence
        identifier values

    """
    if k_mers is None:
        k_mers = {}
    seq_len = len(seq)
    for i_seq in range(seq_len - k_mer_len + 1):
        k_mer = str(seq[i_seq : i_seq + k_mer_len])
        k_mer_rc = Seq(k_mer).reverse_complement()
        if k_mer not in k_mers and k_mer_rc not in k_mers:
            k_mers[k_mer] = {}
            k_mers[k_mer]["src"] = set([seq_id])
            k_mers[k_mer]["cnt"] = 1
        elif k_mer in k_mers:
            k_mers[k_mer]["src"].add(seq_id)
            k_mers[k_mer]["cnt"] += 1
        elif k_mer_rc in k_mers:
            k_mers[k_mer_rc]["src"].add(seq_id)
            k_mers[k_mer_rc]["cnt"] += 1
        else:
            raise Exception("Should not get here")

    return k_mers


def count_k_mers_in_rds(rd_fnm, k_mer_len=25, k_mers=None, seq_rcds=None):
    """Count k-mers of a specified length in a gzipped file of reads
    in the specified format. Optionally update dictionary and list
    returned by this method.

    Parameters
    ----------
    rd_fnm : str
        The name of the gzipped read file
    k_mers : dict
        Empty, or returned by count_k_mers_in_seq()
    seq_rcds : list(Bio.SeqRecord.SeqRecord)
        Empty, or returned by count_k_mers_in_rds()
    k_mer_len : int
        The length of k-mers to count

    Returns
    -------
    k_mers : dict
        Containing k-mer keys and count and source sequence record
        index values
    seq_rcds : list(Bio.SeqRecord.SeqRecord)
        Contains sequence records in which k-mers were counted

    """
    if k_mers is None:
        k_mers = {}
    if seq_rcds is None:
        seq_rcds = []
    read_format, is_gzipped = ru.get_bio_read_format(rd_fnm)
    _open = open
    if is_gzipped:
        _open = gzip.open
    with _open(rd_fnm, "rt") as f:
        seq_rcd_gen = SeqIO.parse(f, format=read_format)
        i_seq = -1
        for seq_rcd in seq_rcd_gen:
            i_seq += 1
            k_mers = count_k_mers_in_seq(seq_rcd.seq, seq_id=i_seq, k_mers=k_mers)
            seq_rcds.append(seq_rcd)

    return k_mers, seq_rcds


def write_k_mer_counts_in_rds(k_mers_in_rd1, k_mers_in_rd2, k_mer_counts_fnm):
    """Write the k-mer counts corresponding to paired reads to the
    specified file.

    Parameters
    ----------
    k_mers_in_rd1 : dct
        Dictionary returned by count_k_mers_in_seq()
    k_mers_in_rd2 : dct
        Dictionary returned by count_k_mers_in_seq()
    k_mer_counts_fnm
        File name to which to write k-mers and their counts

    Returns
    -------
    None

    """
    with open(k_mer_counts_fnm, "w") as f:
        k_mer_set_rd1 = set(k_mers_in_rd1.keys())
        k_mer_set_rd2 = set(k_mers_in_rd2.keys())
        k_mer_set = k_mer_set_rd1.union(k_mer_set_rd2)
        for k_mer in sorted(k_mer_set):
            k_mer_cnt_rd1 = 0
            if k_mer in k_mer_set_rd1:
                k_mer_cnt_rd1 = k_mers_in_rd1[k_mer]["cnt"]
            k_mer_cnt_rd2 = 0
            if k_mer in k_mer_set_rd2:
                k_mer_cnt_rd2 = k_mers_in_rd2[k_mer]["cnt"]
            f.write(
                "{0} {1:10d} {2:10d}\n".format(
                    k_mer, int(k_mer_cnt_rd1), int(k_mer_cnt_rd2)
                )
            )


def read_k_mer_counts(k_mer_counts_fnm, seq_id=0, k_mers=None):
    """Read the k-mer counts corresponding to paired reads from the
    specified file. Optionally update dictionary returned by
    count_k_mers_in_seq().

    Parameters
    ----------
    k_mers_counts_fnm : str
        The name of the file containing k-mers and their counts
    seq_id : int
        The sequence identifer
    k_mers : None or dct
        None or dictionary returned by count_k_mers_in_seq()

    Returns
    -------
    k_mers : dict
        Dictionay containing k-mer keys and counts, and source
        sequence identifiers values
    """
    if k_mers is None:
        k_mers = {}
    with open(k_mer_counts_fnm) as f:
        for ln in f:
            flds = ln.split()
            k_mer = flds[0]
            cnt_rd1 = int(flds[1])
            cnt_rd2 = int(flds[2])
            if k_mer in k_mers:
                logger.ingof(f"Second occurance of k-mer: {k_mer} unexpected")
            else:
                k_mers[k_mer] = {}
                k_mers[k_mer]["src"] = set([seq_id])
                k_mers[k_mer]["rd1"] = cnt_rd1
                k_mers[k_mer]["rd2"] = cnt_rd2

    return k_mers


def collect_k_mer_cnt(k_mers, seq_cnt=2):
    """Collect k-mer counts, assuming that we counted in a doubled
    sequence. Set sequence count otherwise.

    Parameters
    ----------
    k_mers : dict
        Dictionary returned by count_k_mers_in_seq()
    seq_cnt : int
        Number of times the sequence was repeated for counting
        (default: 2)

    Returns
    -------
    numpy.ndarray
        Counts of k_mers

    """
    return np.array([val["cnt"] / seq_cnt for val in k_mers.values()])


def write_paired_reads_for_cnt(
    k_mer_cnt_rd1,
    coverage_rd1,
    k_mers_rd1,
    seq_rcds_rd1,
    k_mer_cnt_rd2,
    coverage_rd2,
    k_mers_rd2,
    seq_rcds_rd2,
):
    """Write paired reads for counts.

    Parameters
    ----------
    k_mer_cnt_rd1 : numpy.ndarray
        Counts of k_mers
    coverage_rd1 : int
        Read one coverage
    k_mers_rd1 : dict
        Contains k-mer keys and count and source sequence record index
        values
    seq_rcds_rd1 : list(Bio.SeqRecord.SeqRecord)
        Contains sequence records in which k-mers were counted
    k_mer_cnt_rd2  : numpy.ndarray
        Counts of k_mers
    coverage_rd2 : int
        Read two coverage
    k_mers_rd2 : dict
        Contains k-mer keys and count and source sequence record index
        values
    seq_rcds_rd2 : list(Bio.SeqRecord.SeqRecord)
        Contains sequence records in which k-mers were counted

    Returns
    -------
    rd1_wr_file_name : str
        File name containing reads one with repeats
    rd1_wo_file_name : str
        File name containing reads one without repeats
    rd2_wr_file_name : str
        File name containing reads two with repeats
    rd2_wo_file_name : str
        File name containing reads two without repeats

    """
    # Find common sequence records
    if len(seq_rcds_rd1) != len(seq_rcds_rd2):
        raise (Exception("Number of sequence records for each pair must agree"))
    i_rcds_rd1, _ = find_seq_rcds_for_cnt(k_mer_cnt_rd1, coverage_rd1, k_mers_rd1)
    i_rcds_rd2, _ = find_seq_rcds_for_cnt(k_mer_cnt_rd2, coverage_rd2, k_mers_rd2)
    i_rcds = i_rcds_rd1.intersection(i_rcds_rd2)

    # Write read one and two files with and without repeats
    rd1_wr_file_name, rd1_wo_file_name = write_reads_for_cnt(
        i_rcds, seq_rcds_rd1, "rd1"
    )
    rd2_wr_file_name, rd2_wo_file_name = write_reads_for_cnt(
        i_rcds, seq_rcds_rd2, "rd2"
    )
    return rd1_wr_file_name, rd1_wo_file_name, rd2_wr_file_name, rd2_wo_file_name


def find_seq_rcds_for_cnt(k_mer_cnt_rds, min_n_clusters=2, max_n_clusters=8):
    pass


def write_reads_for_cnt(i_rcds, seq_rcds, case):
    """Write a FASTQ file containing the sequence records (reads)
    specified, and a file containing all other sequence records
    (reads).

    Parameters
    ----------
    i_rcds : set(int)
        Index of sequence records in a cluster
    seq_rcds : list(Bio.SeqRecord.SeqRecord)
        Contains sequence records
    case : str
        Label for read files

    Returns
    -------
    read_wr_file_name : str
        File name of reads with repeats
    read_wo_file_name : str
        File name of reads without repeats

    """
    # TODO: Settle
    BASE_FILE_NAME = "random_seq"
    NUMBER_PAIRS = 25000

    # Write a FASTQ file containing the specified sequence records
    # (reads)
    read_wr_file_name = BASE_FILE_NAME + "_wr_" + case + ".fastq"
    with open(read_wr_file_name, "w") as f:
        for i_rcd in i_rcds:
            SeqIO.write(seq_rcds[i_rcd], f, "fastq")

    # Write a FASTQ file containing all other sequence records (reads)
    read_wo_file_name = BASE_FILE_NAME + "_wo_" + case + ".fastq"
    with open(read_wo_file_name, "w") as f:
        for i_rcd in set(range(NUMBER_PAIRS)).difference(i_rcds):
            SeqIO.write(seq_rcds[i_rcd], f, "fastq")

    return read_wr_file_name, read_wo_file_name


# TODO: Name and review
def ghi(initial_seq, rd1_file_name, rd2_file_name):
    # TODO: Settle
    BASE_FILE_NAME = "random_seq"
    EXP_CNT = 16
    K_MER_CNT_REP = [EXP_CNT]
    OUTER_DISTANCE = 500
    FRAGMENT_LEN = OUTER_DISTANCE

    # Count k-mers in the random sequence, doubled to represent
    # circular DNA
    # TODO: Fix
    with timing("Counting k_mers in initial sequence"):
        k_mers_seq = count_k_mers_in_seq(initial_seq + initial_seq)

    # Count k-mers in the paired reads, and write the result to a
    # file
    with timing("Counting k_mers in reads"):
        k_mers_rd1, seq_rcds_rd1 = count_k_mers_in_rds(rd1_file_name)
        k_mers_rd2, seq_rcds_rd2 = count_k_mers_in_rds(rd2_file_name)
    with timing("Writing k_mers and counts in reads"):
        write_k_mer_counts_in_rds(k_mers_rd1, k_mers_rd2, BASE_FILE_NAME + "_cnt.txt")

    # Collect k-mer counts, and compute coverage
    k_mer_cnt_seq = collect_k_mer_cnt(k_mers_seq, seq_cnt=2)
    k_mer_cnt_rd1 = collect_k_mer_cnt(k_mers_rd1, seq_cnt=1)
    k_mer_cnt_rd2 = collect_k_mer_cnt(k_mers_rd2, seq_cnt=1)
    coverage_rd1 = int(np.sum(k_mer_cnt_rd1) / np.sum(k_mer_cnt_seq))
    coverage_rd2 = int(np.sum(k_mer_cnt_rd2) / np.sum(k_mer_cnt_seq))

    # Separate reads into those containing a k-mer with the
    # expected count, and all others
    with timing("Writing reads based on read count"):
        (
            rd1_wr_file_name,
            rd1_wo_file_name,
            rd2_wr_file_name,
            rd2_wo_file_name,
        ) = write_paired_reads_for_cnt(
            k_mer_cnt_rd1,
            coverage_rd1,
            k_mers_rd1,
            seq_rcds_rd1,
            k_mer_cnt_rd2,
            coverage_rd2,
            k_mers_rd2,
            seq_rcds_rd2,
            K_MER_CNT_REP,
            EXP_CNT,
        )

    # Assemble paired reads with repeats using SSAKE
    with timing("Assembling paired reads with repeats using SSAKE"):
        ssake_out_base_name = BASE_FILE_NAME + "_ssake"
        if Path(ssake_out_base_name).exists():
            shutil.rmtree(ssake_out_base_name)
            ru.ssake(
                rd1_wr_file_name,
                rd2_wr_file_name,
                ssake_out_base_name,
                FRAGMENT_LEN,
                phred_threshold=20,  # -x
                n_consec_bases=70,  # -n
                ascii_offset=33,  # -d
                min_coverage=5,  # -w
                n_ovrlap_bases=20,  # -m
                n_reads_to_call=2,  # -o
                base_ratio=0.7,  # -r
                n_bases_to_trim=0,  # -t
                contig_size=100,  # -z
                do_track_cvrg=0,  # -c
                do_ignore_mppng=0,  # -y
                do_ignore_headr=0,  # -k
                do_break_ties=0,  # -q
                do_run_verbose=0,
            )  # -v

    # Assemble paired reads with trusted contigs using SPAdes
    with timing(
        "Assembling paired reads repeats and with trusted contigs using SPAdes"
    ):
        spades_wc_out_dir = BASE_FILE_NAME + "_spades_wc"
        trusted_contigs_fnm = str(Path(ssake_out_base_name) / ssake_out_base_name + "_scaffolds.fa")
        if Path(spades_wc_out_dir).exists():
            shutil.rmtree(spades_wc_out_dir)
        ru.spades(
            rd1_file_name,
            rd2_file_name,
            spades_wc_out_dir,
            trusted_contigs_fnm=trusted_contigs_fnm,
        )

