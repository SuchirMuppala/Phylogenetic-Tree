import streamlit as st
import matplotlib.pyplot as plt
from Bio import Entrez, Phylo
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceMatrix
Entrez.email = "25smuppala@gmail.com"

@st.cache_data(show_spinner=False, max_entries=5)
def get_sequence(organism, gene = 'Cytochrome c'):
    query = f"{organism}[Organism] AND {gene}[Protein]"

    with Entrez.esearch(db = 'protein', term = query, retmax = 1) as handle:
        results = Entrez.read(handle)
    results_list = results.get("IdList", [])
    if not results_list:
        return None
    protein_id = results_list[0]
    with Entrez.efetch(db = 'protein', id = protein_id, rettype = 'fasta', retmode = 'text') as handle:
        fasta_results = handle.read()
    lines = [line.strip() for line in fasta_results.splitlines() if not line.startswith('>')]
    return "".join(lines)

st.set_page_config(page_title = 'Phylogenetic Tree Generator', layout = 'wide')
st.title('Phylogenetic Tree Generator')
st.sidebar.header('Input Organisms')
preset = 'Homo sapiens\nGallus gallus\nPan paniscus\nVaranus komodoensis'
with st.sidebar.form(key = 'tree_generator'):
    target = st.text_input('Marker Gene: ', value = 'Cytochrome c')
    organism_input = st.text_area(
        'Enter scientific name of organism (one per line):', 
        value = preset,
        height = 150)
    submit = st.form_submit_button(label = 'Generate Tree')
if submit:
    organisms = [line.strip() for line in organism_input.splitlines() if line.strip()]
    if len(organisms) < 3:
        st.error('You must input 3 or more organisms to form a coherent tree')
    else:
        with st.spinner('Forming tree...'):
            sequences = {}
            fails = []
            for organism in organisms:
                sequence = get_sequence(organism, gene = target)
                if sequence:
                    sequences[organism] = sequence
                else:
                    fails.append(organism)
        if fails:
            st.warning(f'{target} sequence not found for {', '.join(fails)}')
        if len(sequences) >= 3:
            max_len = max(len(s) for s in sequences.values())
            padded_sequences = {org: seq.ljust(max_len, '_') for org, seq in sequences.items()}
            names = list(padded_sequences.keys())
            matrix = []
            for i, name in enumerate(names):
                row = []
                sequence1 = padded_sequences[name]
                for j in range(i + 1):
                    sequence2 = padded_sequences[names[j]]
                    mismatches = sum(1 for a, b in zip(sequence1, sequence2) if a != b)
                    row.append(mismatches / max_len)
                matrix.append(row)

            distance_matrix = DistanceMatrix(names = names, matrix = matrix)
            constructor = DistanceTreeConstructor()
            tree = constructor.upgma(distance_matrix)

            for clade in tree.find_clades():
                if clade.name and clade.name.startswith('Inner'):
                    clade.name = ''
            
            st.subheader(f'Phylogenetic tree for target {target}')
            fig, ax = plt.subplots(figsize = (10, 6))
            Phylo.draw(tree, axes = ax, do_show = False)
            ax.xaxis.set_visible(False)
            ax.yaxis.set_visible(False)
            for spine in ax.spines.values():
                spine.set_visible(False)
            plt.tight_layout()
            st.pyplot(fig)