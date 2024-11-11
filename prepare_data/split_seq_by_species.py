from Bio import SeqIO
import os

# Eingabe-FASTA-Datei
input_fasta = "/home/nila/Dokumente/masterarbeit/my_data_collection/parser_results/my_parser_prot_seq_mollusca.fa"
output_dir = "/home/nila/Dokumente/masterarbeit/my_data_collection/parser_results/sequences_by_species/mollusca/"

# Erstelle das Ausgabeverzeichnis, falls es noch nicht existiert
os.makedirs(output_dir, exist_ok=True)

# Dictionary, um die Sequenzen nach Art zu speichern
species_sequences = {}

# Parsen der FASTA-Datei
for record in SeqIO.parse(input_fasta, "fasta"):
    # Extrahiere den Art-Namen aus der Beschreibung, z.B. den ersten Teil des Headers
    # Annahme: Der Artname ist in der Beschreibung und wird durch ein Leerzeichen oder einen Unterstrich getrennt
    species_name = record.description.split("|")[0]
    
    # Füge die Sequenz zur entsprechenden Art hinzu
    if species_name not in species_sequences:
        species_sequences[species_name] = []
    species_sequences[species_name].append(record)

# Schreibe die Sequenzen jeder Art in separate FASTA-Dateien
for species, sequences in species_sequences.items():
    output_file = os.path.join(output_dir, f"{species}.fasta")
    SeqIO.write(sequences, output_file, "fasta")
    print(f"{len(sequences)} Sequenzen für {species} in {output_file} geschrieben.")

