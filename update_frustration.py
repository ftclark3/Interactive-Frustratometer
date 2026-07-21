import sys
import traceback
import timeit

import numpy as np
from scipy.spatial.distance import pdist, squareform
import numba

import MDAnalysis as mda
from MDAnalysis.analysis import distances

from frustratometer import SparseMatrix, Structure, AWSEM 


# Function to validate input structure
def validate_structure(pdb_filename):
    sys.stdout.write('PYTHON: Validating structure\n')
    sys.stdout.flush()
    u = mda.Universe(pdb_filename)
    chains = [] 
    chain_breaks = [] 
    resids = [] #{}
    zeroindexed_residue_index_counter = -1
    for seg in u.segments: # segments are typically equivalent to chains
        if seg.segid == '':
            # empty segid makes name lookup weird in tcl, so we'll ask the user
            # to use nonempty chain names
            sys.stdout.write(f'# PYTHON ERROR:\n')
            sys.stdout.write('    empty chain names not supported')
            raise Exception('    empty chain names not supported')
        else:
            chains.append(seg.segid)
        resids_for_chain = []
        for residue in seg.residues:
            try:
                assert residue.resid == int(residue.resid) # check that it's not a float
                resid = int(residue.resid) # check that it's not anything that can't be interpreted as an int
            except:
                sys.stdout.write(f'# PYTHON ERROR:\n')
                sys.stdout.write(f"    don't know how to parse resid {residue.resid}")
                raise Exception(f"    don't know how to parse resid {residue.resid}")
            zeroindexed_residue_index_counter += 1
            resids_for_chain.append(resid) # mda's resid equivalent to VMD's resid
        assert len(resids_for_chain) > 1
        for counter in range(len(resids_for_chain)-1):
            if resids_for_chain[counter+1] - resids_for_chain[counter] != 1:
                sys.stdout.write(f'# PYTHON ERROR:\n')
                sys.stdout.write(f"    non-contiguous resids in chain {seg.segid}: {resids_for_chain}")
                raise Exception(f"    non-contiguous resids in chain {seg.segid}: {resids_for_chain}")
        #resids.update({seg.segid: resids_for_chain})
        for resid_for_chain in resids_for_chain:
            resids.append(resid_for_chain)
        chain_breaks.append(zeroindexed_residue_index_counter+1) # chain_breaks is supposed to be the start of the next chain, so we add 1
    chain_breaks = chain_breaks[:-1] # the final index of chain_breaks is the zero-indexed starting residue index of a chain that doesn't exist
    return chains, resids, chain_breaks


# Function to format contacts like a SparseMatrix,
# returning lists of row indices, column indices, and data
@numba.njit(parallel=False) 
def get_sparse_contacts(coords, cutoff=40.0):
    # parallelization would be tricky because we don't know how many elements
    # to pre-allocate to the array
    m = coords.shape[0]
    rows = []
    cols = []
    data = []
    # Numba runs this double loop "instantly"
    for i in range(m):
        for j in range(i + 1, m):
            # Unrolling the math is faster than using np.linalg
            dx = coords[i, 0] - coords[j, 0]
            dy = coords[i, 1] - coords[j, 1]
            dz = coords[i, 2] - coords[j, 2]
            dist = (dx*dx + dy*dy + dz*dz) ** 0.5
            if dist <= cutoff:
                rows.append(i)
                cols.append(j)
                data.append(dist)
    # Return as NumPy arrays for the Frustratometer
    return np.array(rows), np.array(cols), np.array(data)


# Function to build frustratometer server-style output string
# from lists of residue pairs, distances, frustration indices, and sigma_wats 
def get_return_string(Res1s, Res2s, distances, frustration_indices, sigma_wats,
                      resids, chain_breaks):
    # this function uses native python str.join,
    # which is supposed to be the fastest way to concatenate all these strings
    chain_breaks = [0] + chain_breaks # now chain_breaks has the 0-indexed start index of each chain
    N = len(Res1s)
    def loop_helper(counter):
        # determine the style of line that should be drawn, based on the distance
        distance = distances[counter]
        if distance <= 9.5:
            if distance >= 6.5:
                if sigma_wats[counter] >= 0.5:
                    style = 'water-mediated'
                else:
                    style = 'long' # protein-mediated
            elif distance >= 3.5:
                style = 'short' # direct / short-range
            else:
                return '' # no line should be drawn
        else:
            return '' # no line should be drawn
        # determine whether to drawn line or the color of the line, based on the frustration index
        frustration_index = frustration_indices[counter]
        if frustration_index > .78:
            descriptor = 'minimally'
        elif frustration_index < -1.0:
            descriptor = 'highly'        
        else:
            return '' # no line should be drawn
        # convert 0-indexed residue indices to (chain,resid) pairs
        residue1 = Res1s[counter]
        residue2 = Res2s[counter]
        chain1 = None
        chain2 = None
        for chain_index in range(len(chains)-1):
            if chain_breaks[chain_index] <= residue1 < chain_breaks[chain_index+1]:
                chain1 = chains[chain_index]
            if chain_breaks[chain_index] <= residue2 < chain_breaks[chain_index+1]:
                chain2 = chains[chain_index]
        if chain1 is None:  # if the if statements in the above loop are never activated, we're in the last chain
            chain1 = chains[-1]
        if chain2 is None:
            chain2 = chains[-1]
        resid1 = resids[residue1]
        resid2 = resids[residue2]
        return f'{resid1} {resid2} {chain1} {chain2} foo foo foo foo foo foo foo foo {style} {descriptor}\n'
    return ''.join([loop_helper(counter) for counter in range(N)])


# get the desired pdb from tcl
pdb_filename = sys.argv[1]


# validate structure
chains, resids, chain_breaks = validate_structure(pdb_filename)


# get arguments needed to initialize AWSEM class
try:
    # this way of setting up from a Structure won't work if it's 
    # an openawsem structure with IGL/IPR/NGP resnames instead of real resnames
    s = Structure(pdb_filename, repair_pdb=False, sparse=True)
    #s = Structure('2026716185347387546.pdb', repair_pdb=False, sparse=True)
    #s = Structure('aligned_first-openmmawsem.pdb', repair_pdb=False, sparse=True)
    dists = s.distance_matrix 
    #chain_breaks = s.chain_breaks
    sequence = s.sequence
except Exception as e:
    # ASSUME IGL/IPR/NGP resnames
    # 1. Ask Tcl for the Sequence
    sys.stdout.write("__PROMPT__\n")
    sys.stdout.write(f"Topology parser failed ({e}). Please enter the full amino acid sequence in standard one-letter code:\n")
    sys.stdout.flush()
    sequence = sys.stdin.readline().strip()  # Wait for the pop-up answer
    # 2. Ask Tcl for the Chain Breaks
    #sys.stdout.write("__PROMPT__\n")
    #sys.stdout.write("Please enter space-separated chain breaks (0-indexed start residue index of each chain except the first; leave blank if just one chain):\n")
    #sys.stdout.flush()
    #cb_input = sys.stdin.readline().strip()  # Wait for the pop-up answer
    #chain_breaks = [int(x) for x in cb_input.split()] if cb_input else []
    # 3. Build the distance matrix manually
    u = mda.Universe(pdb_filename)
    CACB_sele = u.select_atoms('(name CA and resname IGL) or name CB')
    coords = CACB_sele.positions
    rows, cols, data = get_sparse_contacts(coords, cutoff=40.0) # 40 used for electrostatics
    dists = SparseMatrix(rows, cols, data=data, shape=coords.shape[0])


# initialize AWSEM class 
try:
    f = AWSEM.from_distance_matrix(
        # structural parameters
        dists,
        sequence, #'AYVEIIEQPKQRGMRFRYKCEGRSAGSIPGERSTDTTKTHPTIKINGYTGPGTVRISLVTKDPPHRPHPHELVGKDCRDGYYEADLCPDRSIHSFQNLGIQCVKKRDLEQAISQRIQTNNNPFHVPIEEQRGDYDLNAVRLCFQVTVRDPAGRPLLLTPVLSHPIFDNRAPNTAELKICRVNRNSGSCLGGDEIFLLCDKVQKEDIEVYFTGPGWEARGSFSQADVHRQVAIVFRTPPYADPSLQAPVRVSMQLRRPSDRELSEPMEFQYLPDTDDRHRIEEKRKRTYETFKSIMKKSPFSG',
        chain_breaks=chain_breaks, #chain_breaks= [302], # 1-indexed index of beginning of chain B is 303, so 0-indexed beginning is 302
        full_pdb_distance_matrix = None,
        z_coords=None, # for membrane calculations
        # object initialization parameters
        expose_indicator_functions = False, sparse=True, fast=True, backend='numba', 
        # AWSEM hamiltonian parameters
        min_sequence_separation_contact=7, min_sequence_separation_electrostatics=7,
        min_sequence_separation_rho=2, minimum_sequence_separation=7,
        k_electrostatics=5*4.184)
except Exception as e:
    # The user may try to modify the AWSEM initialization parameters,
    # so AWSEM.from_distance_matrix is one of the most likely points
    # of failure in this script. We'll want to print it to the console
    error_msg = traceback.format_exc() # nicely formatted exception string
    sys.stdout.write(f'# PYTHON ERROR:\n')
    for err_line in error_msg.split('\n'):
        sys.stdout.write(f"#  {err_line}\n")
    raise 


# TELL TCL WE ARE DONE INITIALIZING
sys.stdout.write("__READY__\n")
sys.stdout.flush()


# main loop
while True:
    # This blocks until Tcl sends a line (ending in \n)
    line = sys.stdin.readline()
    if not line: 
        break # Tcl closed the pipe
    try:
        def frustratometry_logic():
            # Process
            floats = [float(x) for x in line.split()]
            #with open('saved_floats.txt','w') as f:
            #    for float_ in floats:
            #        f.write(f'{float_}\n')
            #exit()
            coords = np.array(floats).reshape((-1,3)) # some number of atoms x 3 coordinates each
            dists = pdist(coords) # note this is the flattened upper triangular part of the matrix
            rows, cols, data = get_sparse_contacts(coords, cutoff=40.0) # 40 used for electrostatics
            sm = SparseMatrix(rows, cols, data=data, shape=coords.shape[0])
            f.change_conformation(sm, sm, sm, sm)
            frustration_indices = f.configurational_frustration()[sm.row, sm.col]
            rho1 = f.rho_r[sm.row]
            rho2 = f.rho_r[sm.col]
            sigma_wats = 0.25 * (1 - np.tanh(f.eta_sigma * (rho1 - f.rho_0))) * (1 - np.tanh(f.eta_sigma * (rho2 - f.rho_0)))
            return sm, frustration_indices, sigma_wats
        sm_of_distances, frustration_indices, sigma_wats = frustratometry_logic()
        assert len(sm_of_distances.row) == len(sm_of_distances.col) == len(sm_of_distances.data) == len(frustration_indices) == len(sigma_wats)
        to_return = get_return_string(sm_of_distances.row, sm_of_distances.col, sm_of_distances.data, frustration_indices, sigma_wats,
                                      resids, chain_breaks)
        sys.stdout.write(to_return)
    except Exception:
        error_msg = traceback.format_exc() # nice string format of error
        sys.stdout.write(f'# PYTHON ERROR:\n')
        for err_line in error_msg.split('\n'):
            sys.stdout.write(f"#  {err_line}\n")
    finally:
        sys.stdout.write("__END_OF_CALC__\n")
        sys.stdout.flush()
