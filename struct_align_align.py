import Bio.PDB
import collections, numpy as np, sys

from frag import make_pair_frags
import copy

def matching_frags(frags1, frags2, max_rmsd=2):
    """All pairs of fragments, one from each list, such that their CAs align to within the given maximum RMSD.
    Args:
        frags1, frags2: [ PairFrag ]
        max_rmsd: threshold
    Returns:
        [ (PairFrag, PairFrag) ], all such good pairs, one from frags1 and the other from frags2
    """
    # TODO your code here
    matching_pairs = []

    for frag1 in frags1:
        for frag2 in frags2:
            if frag1.best_rmsd(frag2) <= max_rmsd:   
                matching_pairs.append((frag1, frag2))
    print(matching_pairs)
    return matching_pairs

def corresponding_atoms(atoms1, atoms2, max_rmsd=2):
    """Identify pairs of atoms, one from each list, such that they are within the given maximum RMSD, and each is the closest neighbor of the other.
    For example, if the distance of atoms1[0] and atoms2[10] is d<=2, and no other atoms2 is within d of atoms1[0] and no other atoms1 is within d of atoms2[10], then they correspond.
    Thus only some atoms will have partners, and each atom will have at most 1 partner.
    Args:
        atoms1, atoms2: [ Atom ]
    Returns:
        [ (Atom, Atom) ] satisfying the above criteria
    """
    # TODO your code here
    best_rmsd = 0
    corr_atoms_pairs = []

    for atom1 in atoms1:
        closest_d = max_rmsd+1
        atom1_pair = 0

        for atom2 in atoms2:
            d = atom1.__sub__(atom2) # distance between atom1 and atom2
            if d < closest_d:
                closest_d = d
                atom1_pair = atom2 # find the closest atom2 to atom1
        if closest_d < max_rmsd:
            corr_atoms_pairs.append((atom1, atom1_pair))
    
    return corr_atoms_pairs


    

def rmsd(atoms1, atoms2):
    """RMSD between corresponding lists of atoms (i.e., atoms1[0] vs. atoms2[0], etc.)
    Args:
        atoms1, atoms2: [ Atom ]
    Returns:
        rmsd
    """
    
    # TODO your code here
    #corr_pairs = corresponding_atoms(atoms1, atoms2) #find all the corresponding pairs
    #print(len(corr_pairs))
    sum=0
    for pair in zip(atoms1, atoms2):
        crd1 = pair[0].get_coord()
        crd2 = pair[1].get_coord()
        d = (crd1[0]-crd2[0])**2 + (crd1[1]-crd2[1])**2 + (crd1[2]-crd2[2])**2
        sum += d
    rmsd = (sum/len(atoms1))**0.5
    
    return rmsd



def adjust(moving_chain, fixed_bb, moving_bb):
    """Update moving_chain to get as good a correspondence as possible between its backbone atoms, moving_bb, and those in fixed_bb.
    Correspondence quality is defined in terms of (n,rmsd), where n is the number of corresponding atoms between the two bb atom sets,
    and rmsd is their RMSD. Larger n is always better, and among the largest n, go for smallest RMSD.
    Args:
        moving_chain: Chain
        fixed_bb, moving_bb: [ Atom ], the backbone Atoms from the fixed chain and the moving chain
                             Note that the two lists need not be of the same length.
    Returns:
        best (n,rmsd), and moving_chain has been updated
    """
    
    # TODO your code here
    psup = Bio.PDB.Superimposer()
    psup.set_atoms(fixed_bb, moving_bb)
    psup.apply([a for a in moving_chain.get_atoms()])
    
    return (len(corresponding_atoms(fixed_bb, moving_bb)), rmsd(fixed_bb, moving_bb))


def align(fixed_chain, moving_chain, max_frag_rmsd=2, max_corr_rmsd=2):
    """Update moving_chain to align as best as possible to fixed_chain.
    
    Args:
        fixed_chain, moving_chain: Chain
        max_frag_rmsd: threshold passed to matching_frags
        max_corr_rmsd: threshold passed to corresponding_atoms
    Returns:   
        best (n,rmsd). and moving_chain has been updated
    """
    
    # TODO your code here
    # create paired fragments
    pair_fixed = make_pair_frags(chain = fixed_chain)
    pair_moving = make_pair_frags(chain = moving_chain)
    # find all corresponding pairs between moving and fixed chains
    matching_pairs = matching_frags(frags1=pair_fixed, frags2=pair_moving, max_rmsd=max_frag_rmsd)
    
    # 
    best_n = None  
    best_rmsd = None
    best_pair = matching_pairs[0]
    moving_chain_copy = copy.deepcopy(moving_chain)

    for (frag_fixed, frag_moving) in matching_pairs:
        cur_moving_chain = moving_chain_copy
        imp = Bio.PDB.Superimposer()
        atoms_f=[]; atoms_m = []

        for res in frag_fixed.resis1 + frag_fixed.resis2:
            atoms_f += [a for a in res.get_atoms() if a.get_name() in ['C', 'CA', 'N']]
        for res in frag_moving.resis1 + frag_moving.resis2:
            atoms_m += [a for a in res.get_atoms() if a.get_name() in ['C', 'CA', 'N']]
        
        imp.set_atoms(atoms_f, atoms_m)
        imp.apply(cur_moving_chain)
        
        n = len(corresponding_atoms(atoms_f, atoms_m, max_corr_rmsd))
        r = rmsd(atoms_f, atoms_m)

        new_corr_qual = (n,r)
        
        # Stop when n starts to drop or rm starts to go up
        while n < new_corr_qual[0] or (n == new_corr_qual[0] and r > new_corr_qual[1]):
            n = new_corr_qual[0]
            r = new_corr_qual[1]
            new_corr_qual = adjust(cur_moving_chain, atoms_f, atoms_m)
            
        # Update the best (n, rmsd)
        if not best_n and not best_rmsd:
            best_n = n
            best_rmsd = r
        if n>best_n:
            best_n = n
            best_rmsd = r
            best_pair = (frag_fixed, frag_moving)
        elif n==best_n and r<best_rmsd:
            best_n = n
            best_rmsd = r
            best_pair = (frag_fixed, frag_moving)


    # Update moving chain
    imp = Bio.PDB.Superimposer()

    atoms_f=[]; atoms_m = []

    for res in best_pair[0].resis1 + best_pair[1].resis2:
        atoms_f += [a for a in res.get_atoms() if a.get_name() in ['C', 'CA', 'N']]
    for res in best_pair[0].resis1 + best_pair[1].resis2:
        atoms_m += [a for a in res.get_atoms() if a.get_name() in ['C', 'CA', 'N']]
    
    imp.set_atoms(atoms_f, atoms_m)
    imp.apply([a for a in moving_chain.get_atoms()])

    # Stop when n starts to drop or rm starts to go up
    while n < new_corr_qual[0] or (n == new_corr_qual[0] and r > new_corr_qual[1]):
        n = new_corr_qual[0]
        r = new_corr_qual[1]
        new_corr_qual = adjust(cur_moving_chain, atoms_f, atoms_m)
        
    n = len(corresponding_atoms(atoms_f, atoms_m, max_corr_rmsd))
    r = rmsd(atoms_f, atoms_m)

    new_corr_qual = (n,r)
    # Stop when n starts to drop or rm starts to go up
    while n < new_corr_qual[0] or (n == new_corr_qual[0] and r > new_corr_qual[1]):
        n = new_corr_qual[0]
        r = new_corr_qual[1]
        new_corr_qual = adjust(cur_moving_chain, atoms_f, atoms_m)

    return best_n,best_rmsd

# simplistic command-line driver
# call with arguments
#   fixed.pdb -- filename
#   moving.pdb -- filename
#   moved.pdb -- filename for resulting aligned version of moving.pdb
#   max_frag_rmsd -- parameter to align
#   max_corr_rmsd -- parameter to align
# ex: python align.py data/1ubq.pdb data/3r3mA.pdb out/3r3mA_aligned.pdb 1 1
#     python align.py 1ubq.pdb structs/1a6j.pdb out/1a6j_aligned.pdb 1 1
#     python align.py 1ubq.pdb structs/1a2o.pdb out/1a20_aligned.pdb 1 1


if __name__ == '__main__':
    fixed_fn = sys.argv[1]
    moving_fn = sys.argv[2]
    moved_fn = sys.argv[3]
    max_frag_rmsd = float(sys.argv[4])
    max_corr_rmsd = float(sys.argv[5])

    parser = Bio.PDB.PDBParser(QUIET=True)
    fixed_struct = parser.get_structure('fixed',fixed_fn)
    fixed_chain = next(fixed_struct[0].get_chains())   
    moving_struct = parser.get_structure('moving',moving_fn)
    moving_chain = next(moving_struct[0].get_chains())

    # first couple pieces by themselves
    fixed_frags = make_pair_frags(fixed_chain)
    moving_frags = make_pair_frags(moving_chain)
    print('fragments', len(fixed_frags), len(moving_frags))
    matches = matching_frags(fixed_frags, moving_frags, max_frag_rmsd)
    print('matches',len(matches))

    # align
    print('align', align(fixed_chain, moving_chain, max_frag_rmsd, max_corr_rmsd))

    io = Bio.PDB.PDBIO()
    io.set_structure(moving_struct)
    io.save(moved_fn)
