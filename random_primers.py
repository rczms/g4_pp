import random
import primer3
import pandas as pd

target_temp = 55
monovalent = 17.5
divalent = 10
dntps = 0.5
dna = 50

gc_lower_limit = 40
gc_upper_limit = 60
assay_temp = 37

loops = 10000

def get_random_nucleotide():
    return "".join(random.choice("GATC") for i in range(1))

def get_tm(sequence, mv, dv, dntps, dna):
    return float(primer3.calcTm(sequence, mv_conc = mv, dv_conc = dv, dntp_conc = dntps, dna_conc = dna))

def design_primer():
    random_sequence = get_random_nucleotide()
    melting_temp = get_tm(random_sequence, monovalent, divalent, dntps, dna)
    previous_temp = 0
    while True:
        random_sequence += get_random_nucleotide()
        melting_temp = get_tm(random_sequence, monovalent, divalent, dntps, dna)
        if melting_temp >= target_temp:
            if abs(target_temp - melting_temp) <= abs(target_temp - previous_temp):
                return(random_sequence, melting_temp)
            else:
                return(random_sequence[:-1], previous_temp)
        previous_temp = melting_temp

def gc_content(sequence):
    return float((sequence.count("G") + sequence.count("C"))) / len(sequence) * 100

def is_hairpin(sequence, mv, dv, dntps, dna, temp):
    return primer3.calcHairpin(sequence, mv_conc = mv, dv_conc = dv, dntp_conc = dntps, dna_conc = dna, temp_c = temp).structure_found

results = pd.DataFrame(columns = ["Sequence", "Length", "Melting_temp", "GC"])
for i in range(0, loops):
    primer, primer_tm = design_primer()
    gc = gc_content(primer)
    if gc >= gc_lower_limit and gc <= gc_upper_limit:
        if is_hairpin(primer, monovalent, divalent, dntps, dna, assay_temp) == False:
            results.loc[len(results)] = [primer, len(primer), primer_tm, gc]

results.drop_duplicates(inplace = True)
results.to_csv("primers.csv", sep = ";", decimal = ",")
