import os
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
        return None
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


def main(text, db, sl, folder_path):
    acc_all = get_acc_numbers(text)
    folder_path_updated = create_new_folder(f"{folder_path}/output")

    if sl == "split":
        for acc in acc_all:
            FASTA = fetch_fasta(acc, db)
            if FASTA:
                save_FASTA(FASTA, folder_path_updated, acc)
    else:
        FASTA_all = []
        for acc in acc_all:
            FASTA = fetch_fasta(acc, db)
            if FASTA:
                FASTA_all.append(FASTA)
        save_FASTA(FASTA_all, folder_path_updated, "search_results")
    
    root.destroy()
    root_2 = tk.Tk()

    lb_last = tk.Label(root_2, text = "All done!")
    lb_last.grid(row=0, column=0)

    Q_btn = tk.Button(root_2, text="Great! Thank you!", command = root_2.quit)
    Q_btn.grid(row=1, column=0)
    root_2.mainloop()

    print("All done!")



#a class to get file from tkinter
class g_file_store():
    def __init__(self):
        self.f = "" #variable to store file path to a txt file with seq accessions
        self.t = "" #variable to store the text that will be used to perse seq
        self.db = "" #variable to store data base 
        self.sl = ""
        self.out_f = ""  #variable to store wether user wants one file or separat
    #function to get the path to the file with sequences and save it in self.f
    def g_file(self):
        self.f = filedialog.askopenfilenames(title="Select Text Files", filetypes=[("Text files", "*.txt")])
        if self.f:
            lbl_2 = tk.Label(root, text="The File is located at : " + str(self.f))
            lbl_2.grid(row=0, column=0)
            self.f = self.f[0]

    def o_file(self):
        self.out_f = filedialog.askdirectory(
            title="Save File As")
        if self.out_f:
            lbl_3 = tk.Label(root, text="The output will be saved in " + str(self.out_f))
            lbl_3.grid(row=4, column=0)

    
    #function to get keybord or file input into text as self.t
    def r_text(self):
        tests = [False, False]
        line = ""
        line = entry.get()
        if line:
            self.t = line
            tests[0] = True
        elif self.f:
            try:
                with open(self.f, "r") as file:
                    self.t = file.read()
                    tests[0] = True
            except Exception as e:
                lbl_3 = tk.Label(root, text=f"Error reading file: {e}", fg="red")
                lbl_3.grid(row=4, column=0)
        else:
            lbl = tk.Label(root, text="please input file or text", fg="red")
            lbl.grid(row = 0, column = 0)
        self.db = opt_db_var.get()
        self.sl = opt_sl_var.get()
        if not self.out_f:
            lbl_out = tk.Label(root, text="please, provide a folder to store results", fg="red")
            lbl_out.grid(row=4, column=0)
        else:
            tests[1] = True
        if all(tests):
            main(self.t, self.db, self.sl, self.out_f)
        else:
            print("you didn't provide one of the arguments, please, do")




#create the main tk interface

root = tk.Tk()
root.title("sequence parser")
new_path = g_file_store()

lbl = tk.Label(root, text="Select files containing sequences to parse:")
lbl.grid(row=0, column=0)

lbl_2 = tk.Label(root, text="or input text with sequences to parse:")
lbl_2.grid(row=1, column=0)

lbl_db = tk.Label(root, text="selecte an NCBI database to search in")
lbl_db.grid(row=2, column=0)

lbl_sl = tk.Label(root, text="do you want output as a one fasta file or several (split)?")
lbl_sl.grid(row=3, column=0)

open_button = tk.Button(root, text="Open File", command=new_path.g_file)
open_button.grid(row=0, column=1)

entry = tk.Entry(root)
entry.grid(row=1, column=1)


opt_db_all = ["protein", "nucleotide"]
opt_sl_all = ["one", "split"]


opt_db_var = tk.StringVar(value="protein")
opt_sl_var = tk.StringVar(value="one")

opt_db = tk.OptionMenu(root, opt_db_var, *opt_db_all)
opt_db.grid(row=2, column=1)

opt_sl = tk.OptionMenu(root, opt_sl_var, *opt_sl_all)
opt_sl.grid(row=3, column=1)


#TODO make it work:
lbl_out = tk.Label(root, text="save output file in:")
lbl_out.grid(row=4, column=0)

btn_out = tk.Button(root, text="Select Folder", command=new_path.o_file)
btn_out.grid(row=4, column=1)

submit_button = tk.Button(root, text="Submit", command=new_path.r_text)
submit_button.grid(row=5, column=0)

root.mainloop()

#parsing





#TODO; write code for saving parameters from csv 
#TODO: make variabels in inmput non mandatory, describe defoult parameters in read_me
#TODO: check fnum wether it is correct or not
#TODO: write in read me that if output file is not selectet, the deafult path is current derectory

    




