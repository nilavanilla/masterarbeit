from Bio import SeqIO
from Bio.Seq import Seq
from parser import my_parser, mitos_parser
import os

def write_seq_records(genes, sequences, species):

    seq_records = []
    for i, gene in enumerate(genes):
        if len(sequences[i]) > 0:
            seq_rec = SeqIO.SeqRecord(
                        Seq(sequences[i]),
                        id=species+"|"+gene,
                        description="")
            seq_records.append(seq_rec)

    return seq_records


def collect_data_from_dir(category):

    rna_seq = []
    prot_seq = []
    all_gene_seq = []

    my_parser_out = open('/home/nila/Dokumente/masterarbeit/my_data_collection/parser_results/my_parser_gene_list_'+category+'.txt',"w")
    mitos_parser_out = open('/home/nila/Dokumente/masterarbeit/my_data_collection/parser_results/mitos_parser_gene_list_'+category+'.txt',"w")
    for gb_file in os.listdir("/home/nila/Dokumente/masterarbeit/datenbank/fullRefSeq/"+category+"/"):
        if gb_file.endswith("gb"):
            gb_file_path = os.path.join("/home/nila/Dokumente/masterarbeit/datenbank/fullRefSeq/"+category+"/",gb_file)
            species = gb_file[:-3]
            rnas_my_parser, proteins_my_parser, others_my_parser = my_parser(gb_file_path)
            rnas_mitos_parser, proteins_mitos_parser = mitos_parser(gb_file_path)
            rna_seq = rna_seq + write_seq_records(rnas_my_parser[0], rnas_my_parser[2], species)
            prot_seq = prot_seq + write_seq_records(proteins_my_parser[0], proteins_my_parser[2], species)
            all_gene_seq = all_gene_seq + write_seq_records(rnas_my_parser[0]+proteins_my_parser[0]+others_my_parser[0], rnas_my_parser[1]+proteins_my_parser[1]+others_my_parser[1], species)
            my_parser_out.write(species+'\t[{0}]\t'.format(', '.join(gene for gene in rnas_my_parser[0]))+'[{0}]\t'.format(', '.join(gene for gene in proteins_my_parser[0]))+'[{0}]\t'.format(', '.join(gene for gene in others_my_parser[0]))+'\n')
            mitos_parser_out.write(species+'\t[{0}]\t'.format(', '.join(gene for gene in rnas_mitos_parser))+'[{0}]\t'.format(', '.join(gene for gene in proteins_mitos_parser))+'\n')
    my_parser_out.close()
    mitos_parser_out.close()

    SeqIO.write(rna_seq, '/home/nila/Dokumente/masterarbeit/my_data_collection/parser_results/my_parser_rna_seq_'+category+'.fa', "fasta")
    SeqIO.write(prot_seq, '/home/nila/Dokumente/masterarbeit/my_data_collection/parser_results/my_parser_prot_seq_'+category+'.fa', "fasta")
    SeqIO.write(all_gene_seq, '/home/nila/Dokumente/masterarbeit/my_data_collection/parser_results/my_parser_gene_seq_'+category+'.fa', "fasta")

    print("Parsing "+category+" finished")


def read_gene_list_from_file(path, seperated, with_others):

    proteins = []
    rnas = []
    others = []
    with open(path, 'r') as file:
        for line in file:
            rnas_str = line.split("\t")[1].strip("[]")
            proteins_str = line.strip().split("\t")[2].strip("[]")
            proteins.append([p.strip() for p in proteins_str.split(', ')])
            rnas.append([r.strip() for r in rnas_str.split(', ')])
            if with_others:
                others_str = line.split("\t")[3].strip("[]")
                others.append([o.strip() for o in others_str.split(', ')])
    
    if seperated:
        if with_others:
            return rnas, proteins, others
        return rnas, proteins
    if with_others:
        return [rnas+proteins[i]+others[i] for i, rnas in enumerate(rnas)]
    return [rnas+proteins[i] for i, rnas in enumerate(rnas)]

