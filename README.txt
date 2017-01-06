RNAMoIP 1.1

	This software runs in command line and is written in python2.7
	with the GUROBI5.1 python API.

DEPENDANCIES:
	1) Gurobi (implemented with Gurobi 4.6 free academic licence)
		http://gurobi.com/
		To make sure that installation work, from a terminal,
		you should be able to type:
			gurobi.sh
		and enter in the Gurobi Interactive Shell

INSTALLATION:
    1) The RNAMoIP.py script is located in Src/RNAMoIP.py

	If you don't have a database of motifs in a ".desc" formatm, 
	as outputed bu RNA3DMotifs (http://rna3dmotif.lri.fr/) we 
	provided the one we used. Else jump to execution.

    2) Extract the database i.e.
            tar -xzf CATALOGUE.tgz
    3) The folder No_Redondance_DESC contains the description files,
		as .desc, of the motifs used by RNAMoIP
	   The fold No_Redondance_VIEW3D contains, for each description
		file, all motifs .pdb having the same sequence as 
		described in the .desc.


EXECUTION:

	RNAMoIP.py is only a script that need to be launch with gurobi.sh.
    RNAMoIP.py arguments are as follows:
      REQUIRED:
        -s The rna Sequence 
        -ss The rna Secondary Structure (no pseudoknots) OR
            path to file with list of secondary structures. In that
            case all optimal secondary stuctures will be outputed. 
        -d the path to the .Desc files
      OPTIONAL
        -r (default 0.3) 
           the percentage of basepairs that can be Removed 
        -c (default 4, can be 3)
            the max nb of Components in motifs 
        -m_sols (default 1)
            maximal number of optimal solutions to output
    e.g.

        gurobi.sh RNAMoIP.py -s 'GCCAGGGUGGCAGAGGGGCUUUGCGGCGGACUGCAGAUCCGCUUUACCCCGGUUCGAAUCCGGGCCCUGGC' -ss '(((((((..((((((...)))))).(((((.......))))).....(((((.......))))))))))))' -d 'No_Redondance_DESC' > my_output.txt

        NB: This is a suboptimal structure for the given sequence. It can be obtained with "RNAsubopt -e 3"

        gurobi.sh RNAMoIP.py -s 'GCCAGGGUGGCAGAGGGGCUUUGCGGCGGACUGCAGAUCCGCUUUACCCCGGUUCGAAUCCGGGCCCUGGC' -ss "/path/to/sec/struct/list.tt" -d 'No_Redondance_DESC' -r 0.4 -m_sols 5  > my_output.txt
 
	The file my_output.txt will contain all the usual output 
	of Gurobi. Followed by
	The values of the inserted motifs 
	and deleted base pairs:
    And the last line contains the value of the objective 
    function.
   
    Output exemple (omitting gurobi's) :

################################################################################
Solution for the secondary structure:
    (((((((..((((((...)))))).(((((.......))))).....(((((.......))))))))))))

Optimal solution nb:  1
Corrected secondary structure:
    ((((((...(((.........)))..((((.......))))......((((.........)))).))))))
    C-2DU6.D.2-12-22-1
    C-3CUL.D.6-51-61-1
    C-2DU3.D.3-30-38-1
    C-2DU5.D.1-6-10-1
    C-2DU5.D.1-24-27-2
    C-2DU5.D.1-41-48-3
    C-2DU5.D.1-64-66-4
    D-15-19
    D-14-20
    D-13-21
    D-26-42
    D-52-60
    D-7-65

The optimal solutions has as value:
    -663.0
