from Bio import SeqIO
from Bio.Seq import UndefinedSequenceError
import re
from mitos import mito
from mitos.gb.unify import unify_name_protein, unify_name_rrna
import subprocess


leu1AntiCodons = [["tag", "aag", "cag", "gag", "nag"], ["CUN", "CUA", "CUC", "CUG", "CUA"]]
leu2AntiCodons = [["taa", "caa", "yaa"], ["UUR", "UUA", "UUG"]]
ser1AntiCodons = [["act", "cct", "gct", "tct", "nct"], ["AGN", "AGA", "AGC", "AGG", "AGU"]]
ser2AntiCodons = [["nga", "aga", "cga", "gga", "tga"], ["UCN", "UCA", "UCC", "UCG", "UCT"]]

def leucin_serin_determiner(feature, product):

    gene = feature.qualifiers.get("gene", ["unknown"])[0].lower()
    note = feature.qualifiers.get("note", ["unknown"])[0].lower()
    codon_specif = feature.qualifiers.get("anticodon", ["Unknown"])[0][-4:-1]
    codon_specif = feature.qualifiers.get("codon_recognized", [codon_specif])[0]
    if product == "trnL":
        if any(codon_specif in l for l in leu1AntiCodons) or "1" in note or "1" in gene or any((anticodon in gene or anticodon in note) for anticodon in leu1AntiCodons[0]):
            product = product+"1"
        elif any(codon_specif in l for l in leu2AntiCodons) or "2" in note or "2" in gene or any((anticodon in gene or anticodon in note) for anticodon in leu2AntiCodons[0]):
            product = product+"2"
    if product == "trnS":
        if any(codon_specif in l for l in ser1AntiCodons) or "1" in note or "1" in gene or any((anticodon in gene or anticodon in note) for anticodon in ser1AntiCodons[0]):
            product = product+"1"
        elif any(codon_specif in l for l in ser2AntiCodons) or "2" in note or "2" in gene or any((anticodon in gene or anticodon in note) for anticodon in ser2AntiCodons[0]):
            product = product+"2"
    return product


def append_suffix_to_duplicates(input_list):
    counts = {} 
    output_list = []
    occurrences = {item: input_list.count(item) for item in input_list}
    for item in input_list:
        if occurrences[item] > 1:
            if item in counts:
                counts[item] += 1
            else:
                counts[item] = 0
            new_item = f"{item}_{counts[item]}"
        else:
            new_item = item
        output_list.append(new_item)
    return output_list


def my_parser(gb_file):

    rnas = []
    rnas_seq = []
    rnas_seq_transc = []
    proteins = []
    proteins_seq = []
    proteins_seq_transl =[]
    others = []
    others_seq = []
    with open(gb_file, "r") as input_handle:
        for record in SeqIO.parse(input_handle, "genbank"):

            rRNAs = set()
            for feature in record.features:
                if "rRNA" in feature.type:
                    name = unify_name_rrna(feature.qualifiers.get("product", ["Unknown"])[0], None)
                    if name is None or (name not in ["12S", "15S", "16S", "21S", "23S"]):
                        continue
                    rRNAs.add(name)

            try:
                str(record.seq)
                seq_exists = True
            except UndefinedSequenceError:
                seq_exists = False
            for feature in record.features:
                if "RNA" in feature.type and seq_exists:
                    product = feature.qualifiers.get("product", ["Unknown"])[0]
                    if product == "Unknown":
                        product = product+"_"+feature.type
                    elif "tRNA" == feature.type:
                        if product.strip().split("-")[1].upper() in mito.trnamap.keys():
                            product = mito.trnamap[product.strip().split("-")[1].upper()] 
                            if product in ["trnL","trnS"]:
                                product = leucin_serin_determiner(feature, product)
                    else:
                        product = unify_name_rrna(product, rRNAs) or product.replace(" ", "_").lower()
                    rnas.append(product)
                    rnas_seq.append(feature.extract(record.seq))
                    rnas_seq_transc.append(feature.extract(record.seq).transcribe())

                if feature.type == "CDS" and ("translation" in feature.qualifiers or seq_exists):
                    product = feature.qualifiers["product"][0]
                    gene = feature.qualifiers.get("gene", [product])[0]
                    gene = unify_name_protein(gene) or gene.replace(" ", "_").lower()
                    if "orf" in gene.lower() and gene.lower()[0] in ['f', 'h', 'm']:
                        gene = gene.lower()[0]+"orf"
                    elif "LAGLIDADG" in product.upper():
                        gene = "lagli"
                    elif "GIY" in product.upper():
                        gene = "giy"
                    elif gene == "hypothetical_protein":
                        gene = "hyp"
                    proteins.append(gene)
                    if seq_exists:
                        proteins_seq.append(feature.extract(record.seq))
                    else:
                        proteins_seq.append("")
                    proteins_seq_transl.append(feature.qualifiers.get("translation", [""])[0].replace('J','L'))

                if feature.type == "intron" and seq_exists:
                    gene = feature.qualifiers.get("gene", [""])[0]
                    intron_type = feature.qualifiers.get("note", ["Unknown"])[0]
                    if "group" in intron_type:
                        if "II" in intron_type or "2" in intron_type:
                            intron = "gpII"
                        elif "I" in intron_type or "1" in intron_type:
                            intron = "gpI"
                    else:
                        intron = "intron_"+gene.lower()
                    others.append(intron)
                    others_seq.append(feature.extract(record.seq))
                if feature.type == "misc_feature" and "note" in feature.qualifiers and seq_exists:
                    note = feature.qualifiers.get("note", ["Unknown"])[0]
                    if "control region" in note:
                        others.append("OH")
                        others_seq.append(feature.extract(record.seq))

                    
    return [append_suffix_to_duplicates(rnas), rnas_seq, rnas_seq_transc], [append_suffix_to_duplicates(proteins), proteins_seq, proteins_seq_transl], [append_suffix_to_duplicates(others), others_seq]

    
def mitos_parser(path):

    rnas = []
    proteins = []
    
    output = subprocess.run(["getfeatures.py", "-f", "\"%%bed\"", path], capture_output=True, text=True).stdout.split('.')[:-1]
    for gene in output:
        gene = gene.split()[3]

        if re.search(r"-[a-z]", gene):
            if "-a" in gene:
                gene = re.sub(r"-a.*", "", gene)
            else:
                continue
        if "rn" in gene:
            rnas.append(gene.split("(")[0])
        else:
            proteins.append(gene)

    return rnas, proteins