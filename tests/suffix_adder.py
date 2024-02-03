import os
import csv
import tkinter as tk
from tkinter import filedialog, messagebox
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

class SuffixAdderGUI:
    """Simple GUI for the Suffix Adder Tool"""

    def __init__(self, master):
        """Initialize the GUI components"""
        self.master = master
        master.title("Suffix Adder Tool")

        self.proj_dir_label = tk.Label(master, text="Enter the path to your project folder:")
        self.proj_dir_label.pack()

        self.proj_dir_entry = tk.Entry(master)
        self.proj_dir_entry.pack()

        self.browse_button = tk.Button(master, text="Browse", command=self.browse_folder)
        self.browse_button.pack()

        self.run_button = tk.Button(master, text="Run Suffix Adder", command=self.run_suffix_adder)
        self.run_button.pack()

    def browse_folder(self):
        """Open a file dialog to browse and select the project folder"""
        folder_path = filedialog.askdirectory()
        self.proj_dir_entry.delete(0, tk.END)
        self.proj_dir_entry.insert(0, folder_path)

    def run_suffix_adder(self):
        """Run the Suffix Adder tool using the provided project folder path"""
        proj_dir = self.proj_dir_entry.get()
        if proj_dir and os.path.exists(proj_dir):
            suffix_adder = SuffixAdder(proj_dir)
            suffix_adder.run()
            messagebox.showinfo("Success", "Suffixes Added!")
        else:
            messagebox.showerror("Error", "Invalid project folder path.")

class SuffixAdder:
    """Core functionality for the Suffix Adder Tool"""

    def __init__(self, proj_dir):
        """Initialize with the project directory"""
        self.proj_dir = proj_dir

    def load_suffix_mapping(self):
        """Load suffix mapping from 'occurrences.csv' file"""
        suffix_mapping = {}
        try:
            with open(os.path.join(self.proj_dir, "occurrences.csv"), 'r') as csv_file:
                reader = csv.DictReader(csv_file)
                for row in reader:
                    occurrence_id = row['occurrence_id']
                    country = row['country']
                    suffix_mapping[occurrence_id] = country
            print("Loaded suffix mapping:", suffix_mapping)
            return suffix_mapping
        except Exception as e:
            print(f"An error occurred while loading suffix mapping: {e}")
            return {}

    def map_suffix(self, suffix_mapping):
        """Map suffixes to sequence headers and create a new FASTA file"""
        try:
            with open(os.path.join(self.proj_dir, "seq_master.fasta"), 'r') as input_fasta, \
                 open(os.path.join(self.proj_dir, "suffix_added_seqs.fasta"), 'w') as output_fasta:
                for line in input_fasta:
                    if line.startswith('>'):
                        header = line.strip()[1:]
                        header_parts = header.split("_")
                        sequence_name = header_parts[0]
                        suffix = suffix_mapping.get(sequence_name, '')
                        revised_header_parts = [sequence_name]

                        if len(header_parts) > 1:
                            revised_header_parts.extend(header_parts[1:])

                        if suffix:
                            revised_header_parts.append(suffix)

                        revised_header = ">" + "_".join(revised_header_parts) + "\n"
                        output_fasta.write(revised_header)
                    else:
                        output_fasta.write(line)

            print("Mapping completed successfully!")
        except Exception as e:
            print(f"An error occurred: {e}")

    def run(self):
        """Run the Suffix Adder tool"""
        suffix_map = self.load_suffix_mapping()
        self.map_suffix(suffix_map)

if __name__ == "__main__":
    root = tk.Tk()
    app = SuffixAdderGUI(root)
    root.mainloop()
