import Bio.PDB
import collections, numpy as np, glob, random, os, sys

from frag import make_pair_frags

def generate_frags(filenames):
    """Generate fragments for the pdb files.
    Args:
        filenames: [ string ], each a pdb filename
    Returns:
        [ PairFragment ]
    """
    parser = Bio.PDB.PDBParser(QUIET=True)
    frags = []
    for filename in filenames:
        name = os.path.splitext(os.path.basename(filename))[0]
        struct = parser.get_structure(name, filename)
        chain = next(struct[0].get_chains())
        frags.extend(make_pair_frags(chain))
    return frags

def add_frags(frags, frag_library, max_rmsd):
    """Add the fragments to the fragment library, forming landmarks separated by max_rmsd
    Args:
        frags: [ PairFragment ]
        frag_library: { PairFragment: [ PairFragment ] }, mapping landmark fragments to other fragments similar enough to it
        max_rmsd: threshold for adding a fragment to a landmark's list
    Returns:
        nothing; updates frag_library
    """
    
    # TODO your code here
    for frag in frags:
        sim_found = False
        for key in frag_library:
            if frag.best_rmsd(key) < max_rmsd:
                frag_library[key].append(frag)
                sim_found = True
        if not sim_found:
            frag_library[frag] = []
    

def select_reps(frag_library, min_rmsd=4):
    """Select some representative, diverse landmark fragments
    Args:
        frag_library: { PairFragment: [ PairFragment ] }, mapping landmark fragments to other fragments similar enough to it
        min_rmsd: representatives should be at least this different from each other
    Returns:
        [ PairFragment ]
    """

    # TODO your code here
    # Find the first rep that have the longest list of similar fragments
    longest_rep = 0
    for key, v in frag_library.items():
        if len(v) > longest_rep:
            first_rep = key
            longest_rep = len(v)
    
    ext_reps = [first_rep]
    largest_rmsd = 0
    next_rep = None
    while largest_rmsd > min_rmsd:
        if next_rep:
            ext_reps.append(key)
        
        for key in frag_library.keys():
            # find the largest rmsd between the new potential frag with the existing frags in [reps]
            local_min_rmsd = 999
            
            # closest rep to the key
            for ext_rep in ext_reps:
                cur_rmsd = ext_rep.best_rmsd(key)
                if cur_rmsd < local_min_rmsd: # find the smallest rmsd btw the key and any rep
                    local_min_rmsd = cur_rmsd
                
            # find the greatest rmsd among keys - rmsd is the closest distance btw the key and all the reps
            if local_min_rmsd > largest_rmsd:
                largest_rmsd = local_min_rmsd
                next_rep = key
    
    return ext_reps
            
            



def chimera(target_frags, reps, frag_library, lib_rmsd, max_sub_rmsd):
    """Identify the best library fragment to substitute for each target fragment, when possible
    (i.e. when there is one within the max_sub_rmsd)
    Args:
        target_frags: [ PairFragment ], a list of the fragments from the target protein, for which we want to find substitute fragments from our library
        reps: [ PairFragment ], the list of representatives selected
        frag_library: { PairFragment: [ PairFragment ] }, mapping landmark fragments (including reps) to other fragments similar enough to it
        lib_rmsd: the max_rmsd value used in constructing the library
        max_sub_rmsd: substitutes need to be at most this far away from a selected substitute
    Returns:
        [ (PairFragment, PairFragment) ], the original and the selected subtitute, for cases where that's is possible
                          so the list may be shorter than target_frags
    """

    # TODO your code here
    subt_pairs = []
    for target_frag in target_frags:
        for rep in reps:
            if rep.best_rmsd(target_frag) <= max_sub_rmsd + lib_rmsd: # the visited rep should be close enough w/ the target frag
                cloesest_rmsd = rep.best_rmsd(target_frag)
                subt = rep
                for rep_nbr in frag_library[rep]:
                    if rep_nbr.best_rmsd(target_frag) < cloesest_rmsd:
                        cloesest_rmsd = rep_nbr.best_rmsd(target_frag)
                        subt = rep_nbr
        if cloesest_rmsd <= max_sub_rmsd:
            subt_pairs.append((target_frag, subt))
    return subt_pairs



def write_aligned_frags(filename, frags):
    """Save the frags in a pdb format file, all aligned to the 0th frag
    Args:
        filename: String
        frags: [ PairFragment ]
    """
    psup = Bio.PDB.Superimposer()
    with open(filename,'w') as out:
        for frag in frags:
            psup.set_atoms(frags[0].cas, frag.cas)
            (rot, tran) = psup.rotran
            frag.write(out, rot, tran)

def write_chimera(filename, matches):
    """Save the substitute fragments in a pdb format file, each aligned to its partner from the target
    Args:
        filename: String
        matches: [ (PairFragment, PairFragment) ], listing (target, substitute) matches
    """
    psup = Bio.PDB.Superimposer()
    with open(filename,'w') as out:
        for (orig, sub) in matches:
            psup.set_atoms(orig.cas, sub.cas)
            (rot, tran) = psup.rotran
            sub.write(out, rot, tran)
     
# simplistic command-line driver
# call with arguments
#    struct_dir -- path to the database of structures
#    num_structs -- how many to extract fragments from
#    target_filename -- pdb file to substitute with fragments
#    lib_rmsd, sub_rmsd, rep_rmsd -- thresholds for the various functions
#    out_dir -- where to put output pdb files
# example
# > python fraglib.py structs 10 1ubq.pdb 3 2 3 out
# > python fraglib.py structs 25 1ubq.pdb 3 2 3 out

if __name__=='__main__':
    struct_dir = sys.argv[1]
    num_structs = int(sys.argv[2])
    target_filename = sys.argv[3]
    lib_rmsd = float(sys.argv[4])
    sub_rmsd = float(sys.argv[5])
    rep_rmsd = float(sys.argv[6])
    out_dir = sys.argv[7]

    random.seed(12345) # comment this out if you want different results every time
    struct_filenames = random.sample(glob.glob(struct_dir+'/*.pdb'), num_structs)

    frags = generate_frags(struct_filenames)
    frag_library = collections.defaultdict(list)
    add_frags(frags, frag_library, lib_rmsd)

    reps = select_reps(frag_library, rep_rmsd)
    if len(reps)>100: 
        print('bailing on saving reps because too many')
    else:
        for (i,rep) in enumerate(reps):
            example_frags = [rep]+frag_library[rep]
            if len(example_frags)>50: example_frags = example_frags[:50]
            write_aligned_frags(out_dir+'/r'+str(i)+'.pdb', example_frags)

    parser = Bio.PDB.PDBParser(QUIET=True)
    target_struct = parser.get_structure('target',target_filename)
    target_chain = next(target_struct[0].get_chains())
    target_frags = make_pair_frags(target_chain)

    chimera_frags = chimera(target_frags, reps, frag_library, lib_rmsd, sub_rmsd)
    print('chimera covers',len(chimera_frags),'of',len(target_frags))
    write_chimera(out_dir+'/chimera_'+os.path.basename(target_filename), chimera_frags)

