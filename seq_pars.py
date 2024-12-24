import os
import sys
import re
from Bio import Entrez

# Set email for Entrez
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

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python script.py <file_path> <db_name> <split/one>")
        sys.exit(1)

    file_path = sys.argv[1]
    db = sys.argv[2]
    fnum = sys.argv[3]

    if fnum not in ["split", "one"]:
        print("Error: <split/one> argument must be 'split' or 'one'.")
        sys.exit(1)

    try:
        with open(file_path, "r") as file:
            text = file.read()
    except FileNotFoundError:
        print(f"Error: File '{file_path}' not found.")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading file '{file_path}': {e}")
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

    




