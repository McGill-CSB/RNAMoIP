RNAMoIP 
	This software runs in command line and is written in python2.7
	with the GUROBI4.6 python API.

DEPENDANCIES:
	1) Python2.7 	
	2) Gurobi (implemented with Gurobi 4.6 free academic licence)
		http://gurobi.com/
		To make sure that installation work, from a terminal,
		you should be able to type:
			gurobi.sh
		and enter in the Gurobi Interactive Shell

INSTALLATION:
	1) Download the RNAMoIP.py script

	If you don't have a database of motifs in a ".desc" formatm, 
	as outputed bu RNA3DMotifs (http://rna3dmotif.lri.fr/) we 
	provided the one we used. Else jump to execution.
	4) The folder No_Redondance_DESC contains the description files,
		as .desc, of the motifs used by RNAMoIP
	   The fold No_Redondance_VIEW3D contains, for each description
		file, all motifs .pdb having the same sequence as 
		described in the .desc.

	RNAMoIP.py is only a script that need to be launch with gurobi.sh.
		 ii) The rna secondary structure (without pseudoknots, 
			only use the characters: '(', ')', '.'
RNAMoIP.py takes 5 arguments:
                   i) The rna sequence
                  ii) The rna secondary structure (without pseudoknots,
                         only use the characters: '(', ')', '.'
                 iii) The path to the motifs database ".desc" files.
                         If you downloaded the one provides, it is the
                         full path of the folder "No_Redondance_DESC"
                  iv) The percentage  base pairs that can be removed.
                         This must be a number between 0 and 1.
			 We used "0.3".
                   v) The max number of components in the motif inserted.
                         The implementation has been done up to '4'.
 
 EXAMPLE:
         If you have the sequence / secondary structure:
		GGGCGGCCUUCGGGCUAGACGGUGGGAGAGGCUUCGGCUGGUCCACCCGUGACGCUC
		((((((((....))))..((((..(((..(((....)))..)))..))))...))))
	And you wish to insert up to 4-way junction allowing at most 30% 
	the base pairs to be removed. You need the complete path to your
	motifs folder (if you downloaded the one we provide, it is:
		"No_Redondance_DESC"
		"/path/to/motifs/database/No_Redondance_DESC" \
	You should in a terminal do as follows:
		gurobi.sh MODEL_4_RECOMB.py 'GGGCGGCCUUCGGGCUAGACGGUGGGAGAGGCUUCGGCUGGUCCACCCGUGACGCUC' '((((((((....))))..((((..(((..(((....)))..)))..))))...))))' '/data/pkgs/gurobi/vreinh/MERGED_CATALOGUE/No_Redondance_DESC' 0.3 4 > my_output.txt

	The file my_output.txt will contain all the usual output 
	of Gurobi. The line starting with "Best objective" 
	contains the value of the objective function.
	At the end are the values of the inserted motifs 
	and deleted base pairs:

		 C-1F7Y.B.6-31-38-1
		 C-2KMJ.A.2-7-14-1
		 C-1KUQ.B.5-26-30-1
		 C-2ZJR.X.176-23-25-1
		 C-1FKA.A.44-4-6-1
		 C-1KUQ.B.5-39-42-2
		 C-2ZJR.X.176-44-46-2
		 C-1FKA.A.44-15-19-2
		 C-1FKA.A.44-50-54-3
		 D-8-13
		 D-5-16
		 D-32-37
		 D-27-42

	The rows starting with a "D" are the positions of the deleted base pairs.
	In this case [(8,13), (5,16), (32,37), (27,42)]. The resulting secondary
	structure is thus:
		((((.((......))...((((..((...((......))...))..))))...))))
	The rows with a "C" are the components of the motifs inserted. They are 
	composed of 4 values:
			C-<1>-<2>-<3>-<4>
		<1>: The name of the motifs
		<2>: The first position where the component is inserted
		<3>: The last position where the component is inserted
		<4>: The number of the component in the motifs (starting at 1)
	Therefore in this case we have 5 motifs inserted as follow:
		A) 	1F7Y.B.6 	[(31,38)]
		B) 	2KMJ.A.2 	[(7,14)]
		C) 	1KUQ.B.5 	[(26,30), (39,42)]
		D) 	2ZJR.X.176 	[(23,25), (44,46)]
		E) 	1FKA.A.44 	[(4,6,), (15, 19), (50,54)]
	Where the motifs 1FKA.A.44 has its first components on positions 4-6, 
	its second on 15-19 and its third on 50-54.`
	We can visualize the result as follow:
		GGGCGGCCUUCGGGCUAGACGGUGGGAGAGGCUUCGGCUGGUCCACCCGUGACGCUC
		((((.((......))...((((..((...((......))...))..))))...))))	
		   EEEBBBBBBBBEEEEE   DDDCCCCCAAAAAAAACCCC DDD   EEEEE                      
	Now, e.g. in the folder "No_Redondance_VIEW3D/1F7Y.B.6", all the
	".pdb" files of the motifs having a sequence equal to the sequence
	of 1F7Y.B.6 are deposited. You can use them directly with MC-Sym.
