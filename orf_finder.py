import customtkinter as ctk
from tkinter import filedialog, messagebox
import csv



STOP_CODONS = ['TAA', 'TAG', 'TGA']

CODON_TABLE = {
    'ATA':'I','ATC':'I','ATT':'I','ATG':'M',
    'ACA':'T','ACC':'T','ACG':'T','ACT':'T',
    'AAC':'N','AAT':'N','AAA':'K','AAG':'K',
    'AGC':'S','AGT':'S','AGA':'R','AGG':'R',                 
    'CTA':'L','CTC':'L','CTG':'L','CTT':'L',
    'CCA':'P','CCC':'P','CCG':'P','CCT':'P',
    'CAC':'H','CAT':'H','CAA':'Q','CAG':'Q',
    'CGA':'R','CGC':'R','CGG':'R','CGT':'R',
    'GTA':'V','GTC':'V','GTG':'V','GTT':'V',
    'GCA':'A','GCC':'A','GCG':'A','GCT':'A',
    'GAC':'D','GAT':'D','GAA':'E','GAG':'E',
    'GGA':'G','GGC':'G','GGG':'G','GGT':'G',
    'TCA':'S','TCC':'S','TCG':'S','TCT':'S',
    'TTC':'F','TTT':'F','TTA':'L','TTG':'L',
    'TAC':'Y','TAT':'Y','TAA':'_','TAG':'_','TGA':'_',
    'TGC':'C','TGT':'C','TGG':'W'
}

def reverse_complement(seq):
    complement = {'A':'T','T':'A','C':'G','G':'C'}
    return ''.join(complement.get(base, base) for base in reversed(seq))

def translate_dna(dna):
    return ''.join(CODON_TABLE.get(dna[i:i+3], '?') for i in range(0, len(dna)-2, 3))

def find_orfs(sequence):
    sequence = sequence.upper().replace(" ", "").replace("\n", "")
    rev_seq = reverse_complement(sequence)
    seq_len = len(sequence)
    orfs = []

    
    for frame in range(3):
        i = frame
        while i + 3 <= seq_len:
            if sequence[i:i+3] == 'ATG':
                for j in range(i+3, seq_len, 3):
                    if sequence[j:j+3] in STOP_CODONS:
                        orf_seq = sequence[i:j+3]
                        orfs.append((f'+{frame+1}', i+1, j+3, len(orf_seq), orf_seq, translate_dna(orf_seq)))
                        break
            i += 3

    
    for frame in range(3):
        i = frame
        while i + 3 <= seq_len:
            if rev_seq[i:i+3] == 'ATG':
                for j in range(i+3, seq_len, 3):
                    if rev_seq[j:j+3] in STOP_CODONS:
                        orf_seq = rev_seq[i:j+3]
                        start = seq_len - j
                        end = seq_len - i + 1
                        orfs.append((f'-{frame+1}', start, end, len(orf_seq), orf_seq, translate_dna(orf_seq)))
                        break
            i += 3

    return orfs














ctk.set_appearance_mode("System")  
ctk.set_default_color_theme("blue")  

app = ctk.CTk()
app.title("ORF Finder")
app.geometry("820x720")


input_label = ctk.CTkLabel(app, text="Enter DNA Sequence:", font=ctk.CTkFont(size=14, weight="bold"))
input_label.pack(pady=(10, 0))

input_box = ctk.CTkTextbox(app, height=100, font=("Consolas", 12))
input_box.pack(padx=20, pady=10, fill="x")


def run_orf_finder():
    dna = input_box.get("1.0", "end").strip()
    result_box.delete("1.0", "end")
    if not dna:
        messagebox.showerror("Error", "DNA sequence is empty.")
        return
    orfs = find_orfs(dna)
    if not orfs:
        result_box.insert("end", "No ORFs found.\n")
        return
    for idx, (frame, start, end, length, seq, protein) in enumerate(orfs, 1):
        result_box.insert("end", f"ORF {idx} - Frame {frame} | Start: {start}, End: {end}\n")
        result_box.insert("end", f"Length: {length} bp ({length//3} codons)\n")
        result_box.insert("end", f"Sequence: {seq}\n")
        result_box.insert("end", f"Protein: {protein}\n\n")

def load_fasta():
    path = filedialog.askopenfilename(filetypes=[("FASTA files", "*.fasta *.fa")])
    if path:
        with open(path) as f:
            lines = f.readlines()
            sequence = ''.join(line.strip() for line in lines if not line.startswith('>'))
            input_box.delete("1.0", "end")
            input_box.insert("end", sequence)

def export_results():
    content = result_box.get("1.0", "end").strip()
    if not content:
        messagebox.showwarning("No output", "Nothing to export.")
        return
    path = filedialog.asksaveasfilename(defaultextension=".txt", filetypes=[("Text files", "*.txt"), ("CSV", "*.csv")])
    if path:
        with open(path, "w", newline='') as f:
            f.write(content)


btn_frame = ctk.CTkFrame(app, fg_color="transparent")
btn_frame.pack(pady=5)

ctk.CTkButton(btn_frame, text="Find ORFs", command=run_orf_finder, width=120).grid(row=0, column=0, padx=10)
ctk.CTkButton(btn_frame, text="Load FASTA", command=load_fasta, width=120).grid(row=0, column=1, padx=10)
ctk.CTkButton(btn_frame, text="Export Results", command=export_results, width=120).grid(row=0, column=2, padx=10)


output_label = ctk.CTkLabel(app, text="ORF Results:", font=ctk.CTkFont(size=14, weight="bold"))
output_label.pack(pady=(15, 0))

result_box = ctk.CTkTextbox(app, height=400, font=("Consolas", 11))
result_box.pack(padx=20, pady=10, fill="both", expand=True)


footer = ctk.CTkLabel(app, text="© Abrar's Modern ORF Finder", font=("Segoe UI", 10), text_color="gray")
footer.pack(pady=5)

app.mainloop()













 






#ATGCGATACGCTTGAATGTAGATAGCCATGCTGAAATGTGAATGCGTGA
