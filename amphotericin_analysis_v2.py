import streamlit as st
from rdkit import Chem
import pandas as pd
from rdkit.Chem.Draw import rdMolDraw2D
import base64

def draw_molecule_svg(mol):
    d2d = rdMolDraw2D.MolDraw2DSVG(400, 300)
    d2d.DrawMolecule(mol)
    d2d.FinishDrawing()
    return d2d.GetDrawingText()

def analyze_amphotericin(smiles_string):
    """Analyzes the stereochemistry of a molecule using RDKit."""
    # Load correctly defined molecule from isomeric SMILES
    mol = Chem.MolFromSmiles(smiles_string)
    
    if mol is None:
        return None, "Error: Could not parse the SMILES string."

    # 1. Add hydrogens properly (necessary for full stereochemistry assignment)
    mol_with_hs = Chem.AddHs(mol)
    # 2. Assign stereochemistry properly using RDKit's internal CIP rules
    Chem.AssignStereochemistry(mol_with_hs, force=True, cleanIt=True)
    
    # 3. Identify all chiral centers
    chiral_centers = Chem.FindMolChiralCenters(mol_with_hs, includeUnassigned=True)
    
    # Create Dataframe
    data = []
    has_unassigned = False
    
    for atom_idx, config in chiral_centers:
        atom = mol_with_hs.GetAtomWithIdx(atom_idx)
        element = atom.GetSymbol()
        data.append({"Atom Index": atom_idx, "Element": element, "Configuration": config})
        
        if config == "?":
            has_unassigned = True
            
    df = pd.DataFrame(data)
    
    return mol, df, has_unassigned

def main():
    st.set_page_config(page_title="Amphotericin B Stereochemistry Analysis", page_icon="🧬", layout="centered")
    
    st.markdown(
        """
        <style>
        .student-info {
            position: fixed;
            top: 70px;
            right: 20px;
            text-align: right;
            background-color: rgba(255, 255, 255, 0.1);
            padding: 10px 15px;
            border-radius: 8px;
            box-shadow: 0px 2px 5px rgba(0,0,0,0.1);
            z-index: 9999;
            font-size: 14px;
        }
        /* Handle dark mode gracefully */
        @media (prefers-color-scheme: dark) {
            .student-info {
                background-color: rgba(0, 0, 0, 0.3);
                box-shadow: 0px 2px 5px rgba(255,255,255,0.05);
            }
        }
        </style>
        <div class="student-info">
            <b>Name:</b> Vijay Thomas<br>
            <b>Class:</b> AIML-A<br>
            <b>Roll No:</b> RA2511026050006
        </div>
        """,
        unsafe_allow_html=True
    )
    
    st.title("🧬 Amphotericin B Stereochemistry Analyzer")
    st.markdown("""
    This app analyzes the stereochemistry of Amphotericin B (or any other molecule) using [RDKit](https://www.rdkit.org/).
    It identifies all local chiral centers and outputs their CIP classifications.
    """)
    
    st.markdown("---")
    st.markdown("## 📚 Educational Guide:\n## Stereochemistry Basics")
    
    with st.expander("What is Stereochemistry & Chirality?"):
        st.markdown("""
        **Stereochemistry** is the study of the three-dimensional arrangement of atoms in molecules.
        
        **Chirality** refers to the property of a molecule that makes it non-superimposable on its mirror image (like your left and right hands). A carbon atom bonded to four different groups is a common chiral center.
        """)
        
    with st.expander("Determining R/S Configuration (CIP Rules)"):
        st.markdown("""
        The **Cahn-Ingold-Prelog (CIP) sequence rules** are used to assign R/S configurations:
        1. **Assign priorities** to the four groups attached to the stereocenter based on atomic number (highest atomic number = 1).
        2. **Orient the molecule** so the lowest priority group (4) points away.
        3. **Determine the direction** of groups 1 → 2 → 3:
           - **Clockwise** = **R** (Rectus)
           - **Counterclockwise** = **S** (Sinister)
        """)
        
    with st.expander("About the Example Compound"):
        st.markdown("""
        **Amphotericin B** is a potent antifungal medication. It's a polyene macrolide containing numerous chiral centers along its rigid macrocyclic ring. The precise 3D configuration of all these stereocenters is critical to how the molecule inserts into fungal cell membranes to form pores, leading to cell death. Identifying and correctly assigning every single center is a major challenge in chemical analysis!
        """)
    
    st.markdown("---")

    # Hardcoded stereochemically defined SMILES string for Amphotericin B
    default_smiles = (
        "C[C@H]1[C@@H]([C@H](C[C@H](O1)O[C@@H]2[C@H]([C@@H](C[C@H](O2)C(=O)O)O)N)O)O"
        "[C@H]3/C=C/C=C/C=C/C=C/C=C/C=C/C=C/[C@@H](C[C@H]4[C@@H]([C@@H](O4)C)[C@@H]"
        "([C@@H](C[C@H](C[C@@H](C[C@@H](CC(=O)O[C@@H]([C@H]([C@H]3O)C)C)O)O)O)O)O)O"
    )

    smiles_input = st.text_area("Input SMILES String:", default_smiles, height=100)
    
    if st.button("Analyze Stereochemistry", type="primary"):
        with st.spinner("Analyzing..."):
            mol, df, has_unassigned = analyze_amphotericin(smiles_input)
            
            if mol is None:
                st.error(df) # Here df is the error message
            else:
                st.success("Analysis Complete!")
                
                col1, col2 = st.columns([1, 1])
                
                with col1:
                    st.subheader("2D Structure")
                    svg_string = draw_molecule_svg(mol)
                    b64 = base64.b64encode(svg_string.encode('utf-8')).decode('utf-8')
                    html = f'<img src="data:image/svg+xml;base64,{b64}" width="100%">'
                    st.markdown(html, unsafe_allow_html=True)
                
                with col2:
                    st.subheader("Stereocenters")
                    st.metric("Total Chiral Centers", len(df))
                    st.dataframe(df, use_container_width=True)
                
                if has_unassigned:
                    st.warning("⚠️ Some stereocenters remain unassigned (?) due to incomplete SMILES metadata.")
                else:
                    st.info("✅ Stereochemistry fully resolved.")
                    
if __name__ == "__main__":
    main()
