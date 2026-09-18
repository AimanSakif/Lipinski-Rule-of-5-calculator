import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import pandas as pd
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski
import requests
from bs4 import BeautifulSoup
import time

# Global list to store data for CSV export
data_list = []
sn_counter = 1

def get_smiles_from_imppat(imphy_id):
    """Fetches the SMILES string from the IMPPAT database with enhanced browser mimicking and parsing."""
    # Ensure the ID starts with IMPHY
    if not imphy_id.upper().startswith("IMPHY"):
        imphy_id = "IMPHY" + imphy_id
        
    url = f"https://cb.imsc.res.in/imppat/phytochemical-detailedpage/{imphy_id}"
    
    # Pretend to be a real browser
    headers = {
        'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/115.0.0.0 Safari/537.36',
        'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,*/*;q=0.8',
        'Accept-Language': 'en-US,en;q=0.5'
    }
    
    try:
        print(f"[IMPPAT] Fetching: {url}")
        response = requests.get(url, headers=headers, timeout=25)
        
        if response.status_code != 200:
            print(f"[IMPPAT] HTTP Error {response.status_code}")
            return None, f"HTTP Error {response.status_code}"
            
        soup = BeautifulSoup(response.text, 'html.parser')
        
        # 1. Extract Name
        name = imphy_id
        name_tag = soup.find('strong', string=lambda t: t and 'Phytochemical name' in t)
        if name_tag and name_tag.next_sibling:
            name = name_tag.next_sibling.strip()
            
        # 2. Extract SMILES using multiple fallback strategies
        smiles = None
        
        # Strategy A: Look for the <strong>SMILES:</strong> tag
        smiles_label = soup.find('strong', string=lambda t: t and 'SMILES' in t)
        if smiles_label:
            # Check the immediate text sibling
            if isinstance(smiles_label.next_sibling, str):
                smiles = smiles_label.next_sibling.strip()
            # Check for a <text> tag sibling
            elif smiles_label.find_next_sibling('text'):
                smiles = smiles_label.find_next_sibling('text').get_text(strip=True)
            # Check the parent's text
            else:
                parent_text = smiles_label.parent.get_text(strip=True)
                smiles = parent_text.replace('SMILES:', '').strip()
                
        # Strategy B: If Strategy A fails, search all <text> tags for typical SMILES characters
        if not smiles:
            print("[IMPPAT] Standard SMILES tag not found. Searching all text tags...")
            text_tags = soup.find_all('text')
            for tag in text_tags:
                txt = tag.get_text(strip=True)
                # SMILES strings typically contain C, c, O, N, brackets, etc.
                if any(char in txt for char in ['C', 'c', 'O', 'N', '(', ')']) and len(txt) > 5:
                    smiles = txt
                    print(f"[IMPPAT] Found candidate SMILES in text tag: {smiles[:20]}...")
                    break
                    
        if not smiles:
            print(f"[IMPPAT] Failed to extract SMILES for {imphy_id}.")
            return name, "SMILES Not Found on Page"
            
        print(f"[IMPPAT] Successfully extracted SMILES for {imphy_id}")
        return name, smiles
        
    except requests.exceptions.Timeout:
        print(f"[IMPPAT] Timeout error for {imphy_id}")
        return None, "Connection Timeout"
    except requests.exceptions.ConnectionError:
        print(f"[IMPPAT] Connection error for {imphy_id}")
        return None, "Connection Error (Check Internet)"
    except Exception as e:
        print(f"[IMPPAT] Unexpected Exception: {e}")
        return None, f"Error: {str(e)[:20]}"

def calculate_batch():
    """Reads all inputs from the text box and processes them in a loop."""
    global sn_counter
    input_text = text_input.get("1.0", tk.END).strip()
    
    if not input_text:
        messagebox.showwarning("Input Error", "Please enter at least one identifier.")
        return
        
    identifiers = [line.strip() for line in input_text.split('\n') if line.strip()]
    selected_type = combo_type.get()
    
    btn_add.config(state=tk.DISABLED)
    btn_clear.config(state=tk.DISABLED)
    root.config(cursor="watch")
    
    for identifier in identifiers:
        if selected_type == "Auto-Detect":
            if identifier.upper().startswith("IMPHY"):
                process_type = "IMPHY"
            elif identifier.isdigit():
                process_type = "CID"
            else:
                process_type = "Name"
        else:
            process_type = selected_type
            
        name = None
        smiles = None
        error_msg = ""
        
        try:
            if process_type == "IMPHY":
                name, result = get_smiles_from_imppat(identifier.upper())
                if result and result.startswith("Error") or result in ["Connection Timeout", "Connection Error (Check Internet)", "SMILES Not Found on Page"]:
                    error_msg = result
                else:
                    smiles = result
            elif process_type == "CID":
                compounds = pcp.get_compounds(identifier, 'cid')
                if compounds:
                    name = compounds[0].iupac_name or (compounds[0].synonyms[0] if compounds[0].synonyms else identifier)
                    smiles = compounds[0].canonical_smiles
                else:
                    error_msg = "CID Not Found"
            else: # Name
                compounds = pcp.get_compounds(identifier, 'name')
                if compounds:
                    name = compounds[0].iupac_name or (compounds[0].synonyms[0] if compounds[0].synonyms else identifier)
                    smiles = compounds[0].canonical_smiles
                else:
                    error_msg = "Name Not Found"

            if error_msg or not smiles:
                row = [sn_counter, identifier, error_msg or "No SMILES returned", "", "", "", "", "", "Error"]
                tree.insert("", "end", values=row)
                data_list.append(row)
                sn_counter += 1
                root.update()
                continue

            # Calculate Lipinski properties with RDKit
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                row = [sn_counter, name, smiles, "RDKit Parse Error", "", "", "", "", "Error"]
                tree.insert("", "end", values=row)
                data_list.append(row)
                sn_counter += 1
                root.update()
                continue
                
            mw = Descriptors.MolWt(mol)
            logp = Descriptors.MolLogP(mol)
            hbd = Lipinski.NumHDonors(mol)
            hba = Lipinski.NumHAcceptors(mol)
            
            violations = 0
            if mw > 500: violations += 1
            if logp > 5: violations += 1
            if hbd > 5: violations += 1
            if hba > 10: violations += 1
            
            rule = "Pass" if violations <= 1 else "Fail"
            
            row = [
                sn_counter, name, smiles,
                round(mw, 3), round(logp, 5),
                hbd, hba, violations, rule
            ]
            
            data_list.append(row)
            tree.insert("", "end", values=row)
            sn_counter += 1
            
        except Exception as e:
            row = [sn_counter, identifier, f"Critical Error: {str(e)[:20]}", "", "", "", "", "", "Error"]
            tree.insert("", "end", values=row)
            data_list.append(row)
            sn_counter += 1
            
        root.update()
        time.sleep(0.5) # Be polite to servers to prevent getting blocked
        
    text_input.delete("1.0", tk.END)
    btn_add.config(state=tk.NORMAL)
    btn_clear.config(state=tk.NORMAL)
    root.config(cursor="")
    messagebox.showinfo("Complete", "Batch processing complete! Check the table for results.")

def clear_table():
    global sn_counter, data_list
    for item in tree.get_children():
        tree.delete(item)
    data_list = []
    sn_counter = 1

def save_to_csv():
    if not data_list:
        messagebox.showwarning("No Data", "The table is empty.")
        return
        
    file_path = filedialog.asksaveasfilename(
        defaultextension=".csv",
        filetypes=[("CSV files", "*.csv"), ("All files", "*.*")],
        title="Save results as CSV"
    )
    
    if file_path:
        columns = ["S/N", "Compound", "Canonical SMILES", "Molecular Weight", "LogP", "HBD", "HBA", "Violations", "Lipinski Rule"]
        df = pd.DataFrame(data_list, columns=columns)
        df.to_csv(file_path, index=False)
        messagebox.showinfo("Success", f"Data saved to {file_path}")

# ==========================================
# GUI Setup
# ==========================================
root = tk.Tk()
root.title("Lipinski Rule of 5 Batch Calculator")
root.geometry("1200x750")

input_frame = tk.Frame(root, pady=10, padx=10)
input_frame.pack(fill=tk.X)

text_frame = tk.Frame(input_frame)
text_frame.pack(side=tk.LEFT, fill=tk.Y)

tk.Label(text_frame, text="Enter Identifiers (One per line):", font=("Arial", 10, "bold")).pack(anchor=tk.W)
text_input = tk.Text(text_frame, height=6, width=45, font=("Arial", 10))
text_input.pack(pady=5)

control_frame = tk.Frame(input_frame, padx=20)
control_frame.pack(side=tk.LEFT, fill=tk.Y)

tk.Label(control_frame, text="Input Type:", font=("Arial", 10, "bold")).pack(anchor=tk.W, pady=(0, 5))
combo_type = ttk.Combobox(control_frame, values=["Auto-Detect", "Name", "CID", "IMPHY"], state="readonly", width=18, font=("Arial", 10))
combo_type.current(0)
combo_type.pack(pady=(0, 15))

btn_add = tk.Button(control_frame, text="Calculate All", command=calculate_batch, bg="#4CAF50", fg="white", font=("Arial", 10, "bold"), width=15)
btn_add.pack(pady=5)

btn_clear = tk.Button(control_frame, text="Clear Table", command=clear_table, bg="#f44336", fg="white", font=("Arial", 10, "bold"), width=15)
btn_clear.pack(pady=5)

btn_save = tk.Button(control_frame, text="Save to CSV", command=save_to_csv, bg="#008CBA", fg="white", font=("Arial", 10, "bold"), width=15)
btn_save.pack(pady=5)

table_frame = tk.Frame(root, padx=10, pady=10)
table_frame.pack(fill=tk.BOTH, expand=True)

columns = ("S/N", "Compound", "Canonical SMILES", "Molecular Weight", "LogP", "HBD", "HBA", "Violations", "Lipinski Rule")
tree = ttk.Treeview(table_frame, columns=columns, show="headings")

for col in columns:
    tree.heading(col, text=col)
    if col == "Canonical SMILES":
        tree.column(col, width=400)
    elif col == "Compound":
        tree.column(col, width=200)
    elif col == "Molecular Weight" or col == "LogP":
        tree.column(col, width=120, anchor=tk.CENTER)
    else:
        tree.column(col, width=80, anchor=tk.CENTER)

vsb = ttk.Scrollbar(table_frame, orient="vertical", command=tree.yview)
hsb = ttk.Scrollbar(table_frame, orient="horizontal", command=tree.xview)
tree.configure(yscrollcommand=vsb.set, xscrollcommand=hsb.set)

tree.grid(row=0, column=0, sticky="nsew")
vsb.grid(row=0, column=1, sticky="ns")
hsb.grid(row=1, column=0, sticky="ew")

table_frame.grid_rowconfigure(0, weight=1)
table_frame.grid_columnconfigure(0, weight=1)

root.mainloop()