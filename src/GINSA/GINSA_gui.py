import tkinter as tk
from tkinter import filedialog, messagebox
from GINSA import GINSAClass

class GINSAGUI:
        
    """
    GINSAGUI - GUI for GINSA Analysis

    This class defines a simple GUI for performing analysis using GINSA.

    Parameters:
    - master (tk.Tk): The root window for the GUI.

    Note: Requires the GINSAClass.
    """
    def __init__(self, master):
        """
        Initializes the GINSAGUI instance.
        """
        self.master = master
        self.master.title("GINSA GUI")

        self.project_dir_label = tk.Label(master, text="Project Directory:")
        self.project_dir_label.grid(row=0, column=0, padx=10, pady=5, sticky="e")

        self.project_dir_entry = tk.Entry(master, width=50)
        self.project_dir_entry.grid(row=0, column=1, padx=10, pady=5, columnspan=2, sticky="w")

        self.browse_project_dir_button = tk.Button(master, text="Browse", command=self.browse_project_dir)
        self.browse_project_dir_button.grid(row=0, column=3, padx=5, pady=5, sticky="w")

        self.species_name_label = tk.Label(master, text="Species Name:")
        self.species_name_label.grid(row=1, column=0, padx=10, pady=5, sticky="e")

        self.species_name_entry = tk.Entry(master, width=50)
        self.species_name_entry.grid(row=1, column=1, padx=10, pady=5, columnspan=2, sticky="w")

        self.execute_button = tk.Button(master, text="Execute Analysis", command=self.execute_analysis)
        self.execute_button.grid(row=2, column=1, columnspan=2, pady=10)

    def browse_project_dir(self):
        """
        Opens a file dialog to browse and select the project directory.
        Also updates project directory entry in the GUI.
        """
        directory = filedialog.askdirectory()
        if directory:
            self.project_dir_entry.delete(0, tk.END)
            self.project_dir_entry.insert(0, directory)

    def execute_analysis(self):
        """
        Executes the GINSA analysis using the provided project directory and species name.
        Also displays a message box with the analysis result.
        """
        project_dir = self.project_dir_entry.get()
        species_name = self.species_name_entry.get()

        if not project_dir or not species_name:
            messagebox.showerror("Error", "Please enter both project directory and species name.")
            return

        ginsa_instance = GINSAClass()  # Initialize GINSAClass without parameters; default arguments = None
        try:
            ginsa_instance.main(project_dir, species_name)
            messagebox.showinfo("Analysis Complete", "GINSA analysis completed successfully!")
        except Exception as e:
            messagebox.showerror("Error", f"An error occurred during analysis: {str(e)}")

if __name__ == "__main__":
    root = tk.Tk()
    app = GINSAGUI(root)
    root.mainloop()

