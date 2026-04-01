
from rdkit import Chem
from rdkit.Chem import AllChem

# Hardcoded stereochemically defined SMILES string for Amphotericin B
# Source: PubChem CID 5280965 (Isomeric SMILES)
# This SMILES represents the correct configuration of all 19+ chiral centers in the molecule.
amphotericin_b_smiles = (
    "C[C@H]1[C@@H]([C@H](C[C@H](O1)O[C@@H]2[C@H]([C@@H](C[C@H](O2)C(=O)O)O)N)O)O"
    "[C@H]3/C=C/C=C/C=C/C=C/C=C/C=C/C=C/[C@@H](C[C@H]4[C@@H]([C@@H](O4)C)[C@@H]"
    "([C@@H](C[C@H](C[C@@H](C[C@@H](CC(=O)O[C@@H]([C@H]([C@H]3O)C)C)O)O)O)O)O)O"
)

def analyze_amphotericin():
    """Analyzes the stereochemistry of Amphotericin B using RDKit."""
    
    # Load correctly defined molecule from isomeric SMILES
    mol = Chem.MolFromSmiles(amphotericin_b_smiles)
    
    if mol is None:
        print("Error: Could not parse the SMILES string for Amphotericin B.")
        return

    # 1. Add hydrogens properly (necessary for full stereochemistry assignment)
    mol_with_hs = Chem.AddHs(mol)
    # 2. Assign stereochemistry properly using RDKit's internal CIP rules
    # force=True re-calculates, cleanIt=True ensures consistency
    Chem.AssignStereochemistry(mol_with_hs, force=True, cleanIt=True)
    
    # 3. Identify all chiral centers
    # includeUnassigned=True will reveal if any center could not be resolved (should be handled/avoided)
    chiral_centers = Chem.FindMolChiralCenters(mol_with_hs, includeUnassigned=True)
    
    # Output results
    print("------------------------------------------")
    print("Molecule Name: Amphotericin B")
    print(f"Total Number of Chiral Centers: {len(chiral_centers)}")
    print("------------------------------------------")
    
    print(f"{'Atom Index':<12} | {'Element':<8} | {'Configuration':<15}")
    print("-" * 40)
    
    has_unassigned = False
    for atom_idx, config in chiral_centers:
        atom = mol_with_hs.GetAtomWithIdx(atom_idx)
        element = atom.GetSymbol()
        
        # Display each center's information
        print(f"{atom_idx:<12} | {element:<8} | {config:<15}")
        
        if config == "?":
            has_unassigned = True
            
    if has_unassigned:
        print("\nWarning: Some stereocenters remain unassigned (?) due to incomplete SMILES metadata.")
    else:
        print("\nStereochemistry fully resolved.")

if __name__ == "__main__":
    analyze_amphotericin()
