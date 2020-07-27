from Bio import Seq, SeqIO
import random
import primer3
import pandas as pd
import sys

input_file = "input_sequences.fasta"
insert1 = "ACTACTTACCTCGGGC"
#insert2 = "TTTTCCCCTTTTCCCCTTTTCCCCTTTTCCCC"
insert2 = "CCCCTTTTCCCCTTTTCCCCTTTTCCCC"
#insert2 = "CCAACCCGCCCTACCCAC"
target_length = 80
offset_range = range(2, 25)
design_temperature = 62
monovalent = 25
divalent = 10
dntps = 0
dna = 50

pd.options.display.width = 0

def get_probe(sequence, offset):
    try:
        temporary_sequence = ""
        previous_temp = 0
        for i in sequence[offset:]:
            temporary_sequence += i
            melting_temp = get_tm(temporary_sequence)
            if melting_temp >= design_temperature:
                if abs(design_temperature - melting_temp) <= abs(design_temperature - previous_temp):
                    return(temporary_sequence)
                else:
                    return(temporary_sequence[:-1])
            previous_temp = melting_temp
    except:
        return ""

def get_probes(sequence):
    probe_sequences = []
    probe_lengths = []
    probe_offsets = []
    for offset in offset_range:
        probe1 = ""
        probe2 = ""
        probe1 = get_probe(sequence, offset)
        if probe1 is not "" and probe1 is not None:
            probe2 = get_probe(sequence, offset + len(probe1))
            if probe2 is not "" and probe2 is not None:
                probe_sequences.append((probe1, probe2))
                probe_lengths.append((len(probe1) + len(probe2)) / 2)
                probe_offsets.append(offset)
    if len(probe_sequences) == 0:
        return False
    else:
        index_min = min(range(len(probe_lengths)), key=probe_lengths.__getitem__)
        return probe_sequences[index_min], probe_offsets[index_min]

def design_probes(sequence):
    try:
        (probe1, probe2), offset = get_probes(sequence)
    except:
        (probe1, probe2), offset = ("", ""), ""
    try:
        (probe1_rc, probe2_rc), offset_rc = get_probes(rev_com(sequence))
    except:
        (probe1_rc, probe2_rc), offset_rc = ("", ""), ""
    if probe1 is not "" and probe2 is not "" and len(probe1+probe2) < len(probe1_rc + probe2_rc):
        return (probe1, probe2), offset, 1
    else:
        return (probe1_rc, probe2_rc), offset_rc, -1

def get_tm(sequence):
    return float(primer3.calcTm(sequence, mv_conc = monovalent, dv_conc = divalent, dntp_conc = dntps, dna_conc = dna))

def rev_com(sequence):
    return(str(Seq.Seq(sequence).reverse_complement()))

def fill_sequence(actual_length, target_length):
    fill_len = target_length - actual_length
    if fill_len <= 0:
        return ""
    else:
        return "".join(random.choice("GATC") for i in range(fill_len))

results = pd.DataFrame(columns = ["Name", "Sequence1", "Insert1", "Filler", "Insert2", "Sequence2", "Tm1", "Tm2", "Length", "Design_strand", "5_Offset", "Design_temperature", "Monovalent", "Divalent", "dNTPs", "DNA"])
input_sequences = SeqIO.parse(input_file, "fasta")
for input_sequence in input_sequences:
    (probe1, probe2), offset, strand = design_probes(str(input_sequence.seq))
    pp_length = len(probe1)+len(probe2)+len(insert1)+len(insert2)
    fill_seq = fill_sequence(pp_length, target_length)
    pp_length_filled = pp_length + len(fill_seq)
    results.loc[len(results)] = [str(input_sequence.id), rev_com(probe1), insert1, fill_seq, insert2, rev_com(probe2), get_tm(probe1), get_tm(probe2), pp_length_filled, strand, offset, design_temperature, monovalent, divalent, dntps, dna]

results.to_csv("padlock_probes.csv", sep = ";", decimal = ",")
