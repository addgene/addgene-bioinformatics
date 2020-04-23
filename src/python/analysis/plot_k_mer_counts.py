import os
import pickle
import random
import re
import time

from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import numpy as np

import utilities


def create_r_seq_w_rep(k_mer_len=25, k_mer_cnt=[2**(n+1) for n in range(8)]):
    r_seq = utilities.create_r_seq(k_mer_len)
    for kMC in k_mer_cnt:
        t_seq = utilities.create_r_seq(k_mer_len)
        for iKM in range(kMC):
            r_seq += t_seq
    return r_seq


# TODO: Add error rate, outer distance standard deviation, indel
# fraction, indel extension probablility, and random seed:
# error_rate=0.020, standard_deviation=50, indel_fraction=0.15,
# indel_extended_prob=0.30, random_seed=0,
def simulate_paired_reads_clean(seq, seq_nm, number_pairs=50000,
                                len_first_read=250, len_second_read=250,
                                outer_distance=500):
    rd1_seq_rcds = []
    rd2_seq_rcds = []
    seq_len = len(seq)
    for iP in range(number_pairs):

        i_rd1 = random.randint(0, seq_len - outer_distance)
        rd1_seq_rcd = SeqRecord(
            seq[i_rd1:i_rd1 + len_first_read],
            id=str(iP), name="one", description="first read of read pair")
        rd1_seq_rcd.letter_annotations["phred_quality"] = (
            40 * np.ones(len_first_read, dtype=int)).tolist()
        rd1_seq_rcds.append(rd1_seq_rcd)

        i_rd2 = i_rd1 + outer_distance - len_second_read
        rd2_seq_rcd = SeqRecord(
            seq[i_rd2:i_rd2 + len_first_read],
            id=str(iP), name="one", description="second read of read pair")
        rd2_seq_rcd.letter_annotations["phred_quality"] = (
            40 * np.ones(len_second_read, dtype=int)).tolist()
        rd2_seq_rcds.append(rd2_seq_rcd)

    rd1_fNm = seq_nm + "_r1.fastq"
    rd2_fNm = seq_nm + "_r2.fastq"
    SeqIO.write(rd1_seq_rcds, rd1_fNm, "fastq")
    SeqIO.write(rd2_seq_rcds, rd2_fNm, "fastq")
    return rd1_fNm, rd2_fNm


# TODO: Validate
def simulate_paired_reads_wgsim(seq, seq_nm, number_pairs=50000,
                                len_first_read=250, len_second_read=250,
                                outer_distance=500):
    seq_rcd = SeqRecord(seq, id="0", name="base", description="reference")
    seq_fNm = seq_nm + ".fasta"
    rd1_fNm = seq_nm + "_r1.fastq"
    rd2_fNm = seq_nm + "_r2.fastq"
    SeqIO.write(seq_rcd, seq_fNm, "fasta")
    utilities.wgsim(seq_fNm,
                    rd1_fNm,
                    rd2_fNm,
                    error_rate=0.0,
                    outer_distance=outer_distance,
                    standard_deviation=0,
                    number_pairs=number_pairs,
                    len_first_read=len_first_read,
                    len_second_read=len_second_read,
                    mutation_rate=0.0,
                    indel_fraction=0.0,
                    indel_extended_prob=0.0,
                    random_seed=0,
                    ambiguous_base_frac=1.0,
                    haplotype_mode="")
    return rd1_fNm, rd2_fNm


def count_k_mers_in_seq(seq, k_mers, seq_id=0, k_mer_len=25):
    seq_len = len(seq)
    for i_seq in range(seq_len - k_mer_len + 1):
        k_mer = str(seq[i_seq:i_seq + k_mer_len])
        if k_mer not in k_mers:
            k_mers[k_mer] = {}
            k_mers[k_mer]["src"] = set([seq_id])
            k_mers[k_mer]["cnt"] = 1
        else:
            k_mers[k_mer]["src"].add(seq_id)
            k_mers[k_mer]["cnt"] += 1
    return k_mers


def count_k_mers_in_rds(rd_fNm, k_mers):
    # TODO: Try Biopython read? Seemed very, very slow before
    p_ss = re.compile(r"^@")
    found_identifier = False
    rds = []
    with open(rd_fNm) as f:
        seq_id = -1
        for line in f:
            if p_ss.search(line) is not None:
                found_identifier = True
            elif found_identifier:
                found_identifier = False
                seq_str = line[:-1]
                seq_id += 1
                k_mers = count_k_mers_in_seq(seq_str, k_mers, seq_id)
                # TODO: Will need quality later
                rds.append(seq_str)
    return k_mers, rds


def write_k_mer_counts_in_rds(k_mers_in_rd1, k_mers_in_rd2, k_mer_counts_fNm):
    with open(k_mer_counts_fNm, "w") as f:
        k_mer_set_rd1 = set(k_mers_in_rd1.keys())
        k_mer_set_rd2 = set(k_mers_in_rd2.keys())
        k_mer_set = k_mer_set_rd1.union(k_mer_set_rd2)
        for k_mer in sorted(k_mer_set):
            # TODO: Be careful when errors present
            k_mer_cnt = 0
            if k_mer in k_mer_set_rd1:
                k_mer_cnt += k_mers_in_rd1[k_mer]['cnt']
            if k_mer in k_mer_set_rd2:
                k_mer_cnt += k_mers_in_rd2[k_mer]['cnt']
            f.write("{0} {1:10d}\n".format(k_mer, int(k_mer_cnt / 2)))


def read_k_mer_counts(k_mer_counts_fNm, k_mers, seq_id=0):
    with open(k_mer_counts_fNm) as f:
        for ln in f:
            flds = ln.split()
            k_mer = flds[0]
            cnt = int(flds[1])
            if k_mer in k_mers:
                print("Second occurance of k_mer: {0} unexpected".format(k_mer))
            else:
                k_mers[k_mer] = {}
                k_mers[k_mer]["src"] = set([seq_id])
                k_mers[k_mer]["cnt"] = cnt
    return k_mers


if __name__ == "__main__":

    start_time = time.time()
    print("creating random sequence with repeats ...", end=" ", flush=True)
    r_seq = create_r_seq_w_rep()
    print("done in {0} s".format(time.time() - start_time), flush=True)
    print("random sequence length: {0}".format(len(r_seq)))
    print("r_seq: ", r_seq)
    print("r_seq + r_seq: ", r_seq + r_seq)

    # TODO: Double or single?
    start_time = time.time()
    print("counting k_mers in sequence ...", end=" ", flush=True)
    k_mers_in_seq = count_k_mers_in_seq(r_seq + r_seq, {})
    print("done in {0} s".format(time.time() - start_time), flush=True)

    # === clean ===

    # TODO: Double or single?
    start_time = time.time()
    print("simulating paired reads ...", end=" ", flush=True)
    rd1_fNm, rd2_fNm = simulate_paired_reads_clean(r_seq + r_seq, "ccc_clean")
    print("done in {0} s".format(time.time() - start_time), flush=True)

    start_time = time.time()
    print("counting k_mers in reads ...", end=" ", flush=True)
    k_mers_in_rd1_c, _ = count_k_mers_in_rds(rd1_fNm, {})
    k_mers_in_rd2_c, _ = count_k_mers_in_rds(rd2_fNm, {})
    print("done in {0} s".format(time.time() - start_time), flush=True)

    start_time = time.time()
    print("writing k_mers and counts in reads ...", end=" ", flush=True)
    write_k_mer_counts_in_rds(k_mers_in_rd1_c, k_mers_in_rd2_c, "ccc_clean_cnt.txt")
    print("done in {0} s".format(time.time() - start_time), flush=True)

    start_time = time.time()
    print("counting k_mers in reads using kmc ...", end=" ", flush=True)
    read_file_names = [rd1_fNm, rd2_fNm]
    inp_database_file_name = "ccc_clean"
    out_database_file_name = "ccc_clean_kmc.txt"
    utilities.kmc(read_file_names, inp_database_file_name,
                  count_min=0, max_count=1e9, canonical_form=False)
    utilities.kmc_transform(inp_database_file_name, "dump",
                            out_database_file_name)
    print("done in {0} s".format(time.time() - start_time), flush=True)

    start_time = time.time()
    print("reading k_mers counted using kmc ...", end=" ", flush=True)
    k_mers_by_kmc_c = read_k_mer_counts(out_database_file_name, {})
    print("done in {0} s".format(time.time() - start_time), flush=True)

    # === wgsim ===

    start_time = time.time()
    print("simulating paired reads ...", end=" ", flush=True)
    rd1_fNm, rd2_fNm = simulate_paired_reads_wgsim(r_seq + r_seq, "ccc_wgsim")
    print("done in {0} s".format(time.time() - start_time), flush=True)

    start_time = time.time()
    print("writing k_mers and counts in reads ...", end=" ", flush=True)
    k_mers_in_rd1_n, _ = count_k_mers_in_rds(rd1_fNm, {})
    k_mers_in_rd2_n, _ = count_k_mers_in_rds(rd2_fNm, {})
    print("done in {0} s".format(time.time() - start_time), flush=True)

    start_time = time.time()
    print("writing k_mers and counts in reads ...", end=" ", flush=True)
    write_k_mer_counts_in_rds(k_mers_in_rd1_n, k_mers_in_rd2_n, "ccc_wgsim_cnt.txt")
    print("done in {0} s".format(time.time() - start_time), flush=True)

    start_time = time.time()
    print("counting k_mers in reads using kmc ...", end=" ", flush=True)
    read_file_names = [rd1_fNm, rd2_fNm]
    inp_database_file_name = "ccc_wgsim"
    out_database_file_name = "ccc_wgsim_kmc.txt"
    utilities.kmc(read_file_names, inp_database_file_name,
                  count_min=0, max_count=1e9, canonical_form=False)
    utilities.kmc_transform(inp_database_file_name, "dump",
                            out_database_file_name)
    print("done in {0} s".format(time.time() - start_time), flush=True)

    start_time = time.time()
    print("reading k_mers counted using kmc ...", end=" ", flush=True)
    k_mers_by_kmc_n = read_k_mer_counts(out_database_file_name, {})
    print("done in {0} s".format(time.time() - start_time), flush=True)

    for k_mer in k_mers_in_seq.keys():

        ln = "k_mer: {0} in_seq: {1:6d}".format(k_mer, k_mers_in_seq[k_mer]["cnt"])

        if k_mer in k_mers_in_rd1_c.keys():
            ln += ", in_rd1_c: {0:6d}".format(k_mers_in_rd1_c[k_mer]["cnt"])
        else:
            ln += ", in_rd1_c: {0:6d}".format(0)

        if k_mer in k_mers_in_rd2_c.keys():
            ln += ", in_rd2_c: {0:6d}".format(k_mers_in_rd2_c[k_mer]["cnt"])
        else:
            ln += ", in_rd2_c: {0:6d}".format(0)

        if k_mer in k_mers_by_kmc_c.keys():
            ln += ", by_kmc_c: {0:6d}".format(k_mers_by_kmc_c[k_mer]["cnt"])
        else:
            ln += ", by_kmc_c: {0:6d}".format(0)

        if k_mer in k_mers_in_rd1_n.keys():
            ln += ", in_rd1_n: {0:6d}".format(k_mers_in_rd1_n[k_mer]["cnt"])
        else:
            ln += ", in_rd1_n: {0:6d}".format(0)

        if k_mer in k_mers_in_rd2_n.keys():
            ln += ", in_rd2_n: {0:6d}".format(k_mers_in_rd2_n[k_mer]["cnt"])
        else:
            ln += ", in_rd2_n: {0:6d}".format(0)

        if k_mer in k_mers_by_kmc_n.keys():
            ln += ", by_kmc_n: {0:6d}".format(k_mers_by_kmc_n[k_mer]["cnt"])
        else:
            ln += ", by_kmc_n: {0:6d}".format(0)

        print(ln)
