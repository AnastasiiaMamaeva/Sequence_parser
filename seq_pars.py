import os
import sys
import re
import tkinter as tk
from tkinter import filedialog
from Bio import Entrez


Entrez.email = "mamaeva.anastas@gmail.com"

def is_ncbi_accession(text):
    accession_pattern = r"^[A-Z]{1,2}_\d+(\.\d+)?$"
    return bool(re.match(accession_pattern, text))

def get_acc_numbers(text):
    splited_text = re.split(r"[;,\| .\t\n]", text)
    acc = [i for i in splited_text if is_ncbi_accession(i)]
    print(f"Accession numbers:\n{acc}")
    return set(acc)

def fetch_fasta(accession, db):
    try:
        with Entrez.efetch(db=db, id=accession, rettype="fasta", retmode="text") as h:
            fasta_data = h.read()
            return fasta_data
    except Exception as e:
        print(f"Error fetching data for {accession}: {e}")
        return None

def create_new_folder(folder_path, iter=0):
    if iter > 100:
        print("Error: Unable to create a new folder after 100 attempts.")
        sys.exit(1)
    name = f"{folder_path}_{iter}" if iter > 0 else folder_path
    try:
        os.mkdir(name)
        print(f"Folder '{name}' created successfully.")
        return name
    except FileExistsError:
        return create_new_folder(folder_path, iter + 1)

def save_FASTA(FASTA, folder_path, name):
    file_path = os.path.join(folder_path, f"{name}.fasta")
    try:
        with open(file_path, "w") as file:
            if isinstance(FASTA, list):
                for text in FASTA:
                    if text:
                        file.write(text + "\n")
            else:
                if FASTA:
                    file.write(f"{FASTA}\n")
        print(f"FASTA file saved as {file_path}")
    except Exception as e:
        print(f"Error saving FASTA data: {e}")

#a class to get file from tkinter
class g_file_store():
    def __init__(self):
        self.f = "" #variable to store file path to a txt file with seq accessions
        self.t = "" #variable to store the text that will be used to perse seq
        self.db = "" #variable to store data base 
        self.sl = "" #variable to store wether user wants one file or separat
    #function to get the path to the file with sequences and save it in self.f
    def g_file(self):
        self.f = filedialog.askopenfilenames(title="Select Text Files", filetypes=[("Text files", "*.txt")])
        if self.f:
            lbl_2 = tk.Label(root, text="The File is located at : " + str(self.f), font=('Aerial 11'))
            lbl_2.grid(row = 2, column = 0)
    
    #function to get keybord or file input into text as self.t
    def r_text(self):
        line = ""
        line = entry.get()
        if line:
            self.t = line

        elif self.f:
            with open(self.f, "r") as file:
                self.t = file.read()
        else:
            lbl_3 = tk.Label(root, text="please input file or text", font=('Aerial 11'))
            lbl_3.grid(row = 5, column = 0)
        self.db = opt_db.get()
        self.sl = opt_sl.get()




if __name__ == "__main__":
    
    #create a storage for file
    
    
    #create the main tk interface
    
    root = tk.Tk()
    root.title("sequence parser")
    new_path = g_file_store()


    lbl = tk.Label(root, text="Select files containing sequences to parse:")
    lbl.grid(row=0, column=0)

    open_button = tk.Button(root, text="Open Files", command=new_path.g_file)
    open_button.grid(row=0, column=0)

    entry = tk.Entry(root)
    entry.grid(row=3, column=0)
    

    opt_db_all = ["protein", "nucleotide"]
    opt_db = tk.OptionMenu(root, tk.StringVar().set("protein"), *opt_db_all)
    opt_db.grid(row =1, column = 1 )

    opt_sl_all = ["one", "split"]
    opt_sl = tk.OptionMenu(root, tk.StringVar().set("one"), *opt_sl_all)
    opt_sl.grid(row = 2, column = 1)
    submit_button = tk.Button(root, text="Submit", command=new_path.r_text)
    submit_button.grid(row=4, column=0)

    root.mainloop()


    db = sys.argv[2]
    fnum = sys.argv[3]

    if fnum not in ["split", "one"]:
        print("Error: <split/one> argument must be 'split' or 'one'.")
        sys.exit(1)


    acc_all = get_acc_numbers(text)

    if not acc_all:
        print("No accession numbers found in the input file.")
        sys.exit(1)

    folder_path = create_new_folder("FASTA_files")

    if fnum == "split":
        for acc in acc_all:
            FASTA = fetch_fasta(acc, db)
            if FASTA:
                save_FASTA(FASTA, folder_path, acc)
    else:
        FASTA_all = []
        for acc in acc_all:
            FASTA = fetch_fasta(acc, db)
            if FASTA:
                FASTA_all.append(FASTA)
        save_FASTA(FASTA_all, folder_path, "search_results")

    print("All done!")



#TODO; write code for saving parameters from csv 
#TODO: make variabels in inmput non mandatory, describe defoult parameters in read_me
#TODO: check fnum wether it is correct or not

    




