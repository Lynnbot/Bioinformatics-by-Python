import Bio.PDB, numpy as np, sys


class PairFrag (object):
    """A non-contiguous pair of contiguous residues, such that the middle residues are in contact.
    E.g., residues 1-5 and 36-40, if 3 and 38 are in contact.
    Instance variables:
        chain: Chain, from which the fragment is extracted
        resis1, resis2: [ Residue ], the discontiguous pair of contiguous residues
        cas: [ Atom ], a list of all the CA atoms in the fragment (useful for alignment purposes)
    """
    
    psup = Bio.PDB.Superimposer() # to be used in best_rmsd
    
    def __init__(self, chain, resis1, resis2):
        self.resis1 = resis1
        self.resis2 = resis2
        self.cas = [r['CA'] for r in resis1+resis2]
        self.name = '%s %d*%d' % (chain.get_full_id()[0], resis1[len(resis1)//2].get_id()[1], resis2[len(resis2)//2].get_id()[1])
        
    def __repr__(self):
        return self.name
    
    def __str__(self):
        return self.name

    def best_rmsd(self, frag2):
        """A convenience function: returns the best RMSD when self is aligned with frag2"""
        PairFrag.psup.set_atoms(self.cas, frag2.cas)
        return PairFrag.psup.rms

    def write(self, out, rot=None, tran=None):
        """Writes the ATOM records for this backbone of the fragment
        Args:
            out: a file handle
            rot, tran: possibly a transformation to apply to each atom
        """
        atom_num = 1; res_num = 1
        for resi in self.resis1 + self.resis2:
            for atom in resi.get_atoms():
                if atom.get_name() in ['N','CA','C']:
                    # hacked for our special case, ignoring bfactor etc.
                    coord = atom.coord
                    # determine transformed version within actually modifying the atom (which can be in multiple frags)
                    if rot is not None: coord = np.dot(atom.coord, rot)
                    if tran is not None: coord += tran
                    out.write('ATOM   %4d  %-3s %3s A %3d    %8.3f%8.3f%8.3f\n' %(atom_num, atom.get_name(), resi.get_resname(), res_num, round(coord[0],3), round(coord[1],3), round(coord[2],3)))
                    atom_num += 1
            res_num += 1
        out.write('END\n')

def make_pair_frags(chain, contact_thresh=8):
    """Generates all the legit PairFrag instances for the chain, such that the middle residue CAs are within the contact threshold
    Args:
        chain: Chain
        contact_thresh: max distance to be considered a contact
    Returns: 
        [ PairFrag ]
    """    
    # TODO your code here
    residues = [r for r in chain if r.get_resname() in Bio.PDB.Polypeptide.d3_to_index]
    atom_ca = [atom for atom in chain.get_atoms() if atom.get_name() == 'CA']   

    start_id = residues[0].get_id()[1]
    end_id = residues[-1].get_id()[1]

    # use the efficient contact finder
    searcher = Bio.PDB.NeighborSearch(atom_ca)
    atom_pairs = searcher.search_all(contact_thresh) # threshold of 8

    pairfrags = []
    for (a,b) in atom_pairs:
        i = a.get_parent().get_id()[1] - 1 # converting to 0 indexed
        j = b.get_parent().get_id()[1] - 1  # converting to 0 indexed
        i_loc = i - start_id
        j_loc = j - start_id

        if i_loc<2 or i > end_id-2 or j_loc<2 or j > end_id-2:
            continue
        if abs(i_loc-j_loc)<5:
            continue
        if j_loc+3 > len(residues) or i_loc+3 > len(residues):
            continue

        pairfrag = PairFrag(chain=chain, resis1 = residues[i_loc-2: i_loc+3], resis2= residues[j_loc-2:j_loc+3])
        pairfrags.append(pairfrag)
    
    return pairfrags



# simplistic command-line driver
# call with arguments
#   in_filename -- will use the first chain from model 0 in this pdb file
#   out_filename -- store the frags in this pdb file
#   contact_thresh
# ex: python frag.py data/1ubq.pdb out/1ubq_frags.pdb 8

if __name__ == '__main__':
    in_filename = sys.argv[1]
    out_filename = sys.argv[2]
    contact_thresh = float(sys.argv[3])

    parser = Bio.PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('test',in_filename)
    chain = next(structure[0].get_chains()) # get the first chain of the first model
        
    frags = make_pair_frags(chain, contact_thresh)
    print(len(frags), 'fragments')

    with open(out_filename,'w') as out:
        for frag in frags:
            frag.write(out)

    for frag in sorted(frags, key=str):
        print(frag)

