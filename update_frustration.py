import sys

import numpy as np

from frustratometer import Structure, AWSEM 


s = Structure('2026716185347387546.pdb', repair_pdb=False)
f = AWSEM(s)


while True:
    # This blocks until Tcl sends a line (ending in \n)
    line = sys.stdin.readline()
    if not line: 
        break # Tcl closed the pipe
        
    # Process
    floats = [float(x) for x in line.split()]
    #result = sum(floats)
    
    # Reply back to Tcl -- need to update this to actually compute frustration
    #                      instead of loading precomputed file
    #sys.stdout.write(f"{result}\n")
    with open('2026716185347387546.pdb_configurational', 'r') as f:
        stuff = ''.join(f.readlines())
    sys.stdout.write(stuff)
    sys.stdout.write("__END_OF_CALC__\n")
    sys.stdout.flush()
