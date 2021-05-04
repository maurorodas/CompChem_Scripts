#!/bin/bash

adtToolsPath='/usr/local/mgltools/MGLToolsPckgs/AutoDockTools/Utilities24'


centerX=50.3948
centerY=-37.9373
centerZ=17.3044

sizeX=15
sizeY=17
sizeZ=14

spacingVal=1.000
gaEvals=750000

# file with proteins and chain to perform a redocking calculation
filename='proteins.dat'
# Name in the pdb file of the ligand to use for the redocking
ligandName='BRL'
# Number of the proteins to perform the redocking calculation
protNumber=5

# Directory to save all ligands 
if [ -d ligandsDir ]
then
    rm -r ligandsDir
    mkdir ligandsDir
else
    mkdir ligandsDir
fi

for line in `seq 1 $protNumber`
do
    # To obtain the pdbcode of the file with proteins
    pdbcode=`awk -v var=$line '{if(NR==var) print $1}' $filename`
    # To obtain the chain letter of the file with proteins
    chain=`awk -v var=$line '{if(NR==var) print $2}' $filename`
    # To contruct a name with pdbcode and chain to fetch only de chain of interest
    my_query=`echo "$pdbcode$chain"`
    # The script creates a directory for each protein with the name protein1, protein2, and so on
    fileProt=`echo "protein$line"`
    # The script extract all ligands for each protein and names them as ligand1, ligand2, and so on 
    ligand=`echo "ligand$line"`

    # Checks if the protein directory exists, if true remove it and make it a new directory
    # if false make a new directory
    if [ -d $fileProt ]
    then
	rm -r $fileProt
	mkdir $fileProt
    else
	mkdir $fileProt
    fi


    # These variables are necessary to use the pdb code of reference protein for the alignment step
    # with all the other proteins, the reference protein is the first protein in the file
    if [ $line -eq 1 ]
    then
	protRef=$pdbcode
	chainRef=$chain
	query_ref=`echo "$protRef$chainRef"`
    fi

    if [ $line -eq 1 ]
    then
	# This is the reference protein, only fetch his pdb file, extract the ligand and performs a basic cleaning
	# note that this protein does not need and alignment, is the reference!
	cd $fileProt
	pymol -c -d 'fetch '$my_query', type=pdb, async=0; save '$pdbcode'.pdb'
	grep $ligandName $pdbcode.pdb > $ligand.pdb
	grep -v -e HOH -e $ligandName -e CONECT -e HETATM $pdbcode.pdb > temp_$pdbcode.pdb
	mv temp_$pdbcode.pdb $pdbcode.pdb
	# Using reduce from ambertools to add hydrogens
	ligand_h=`echo $ligand"_H"`
	protein_h=`echo $pdbcode"_H"`
	reduce $ligand.pdb > $ligand_h.pdb
	reduce $pdbcode.pdb > $protein_h.pdb
	cp $ligand_h.pdb ../ligandsDir/.
	cd ..
    else
	# Performs the aligment using the "align" method of Pymol, saves the protein with the new coordinates,
	# extract the ligand and make a basic cleaning
	cd $fileProt
	pymol -c -d 'fetch '$my_query' '$query_ref', type=pdb, async=0; align '$my_query', '$query_ref'; save '$pdbcode'.pdb, '$my_query''
	grep $ligandName $pdbcode.pdb > $ligand.pdb
	grep -v -e HOH -e $ligandName -e CONECT -e HETATM $pdbcode.pdb > temp_$pdbcode.pdb
	mv temp_$pdbcode.pdb $pdbcode.pdb
	# Using reduce from ambertools to add hydrogens
	ligand_h=`echo $ligand"_H"`
	protein_h=`echo $pdbcode"_H"`
	reduce $ligand.pdb > $ligand_h.pdb
	reduce $pdbcode.pdb > $protein_h.pdb
	cp $ligand_h.pdb ../ligandsDir/.
	cd ..
    fi      
done

line=''

for line in `seq 1 $protNumber`
do
    # To obtain the pdbcode of the file with proteins
    pdbcode=`awk -v var=$line '{if(NR==var) print $1}' $filename`
    # The script creates a directory for each protein with the name protein1, protein2, and so on
    fileProt=`echo "protein$line"`
    # The script extract all ligands for each protein and names them as ligand1, ligand2, and so on 
    protein_h=`echo $pdbcode"_H"`

    cp ligandsDir/ligand* $fileProt/.

    cd $fileProt
    pythonsh $adtToolsPath/prepare_receptor4.py -r $protein_h.pdb -A hydrogens -o $protein_h.pdbqt
    
    ligand=''
    for lig in `seq 1 $protNumber`
    do
	ligand=`echo "ligand$lig"`	
	if [ -d $ligand ]
	then
	    rm -r $ligand
	    mkdir $ligand
	else
	    mkdir $ligand
	fi

	ligand_h=`echo $ligand"_H"`
	
	mv $ligand_h.pdb $ligand/.
	cp $protein_h.pdbqt $ligand/.

	cd $ligand

	pythonsh $adtToolsPath/prepare_ligand4.py -l $ligand_h.pdb -o $ligand_h.pdbqt
	ligandRandom=`echo $ligand"_random"`
	pythonsh $adtToolsPath/write_random_state_ligand.py -l $ligand_h.pdbqt -o $ligandRandom.pdbqt
	
	pythonsh $adtToolsPath/prepare_gpf4.py -l $ligandRandom.pdbqt -r $protein_h.pdbqt -p npts="$sizeX,$sizeY,$sizeZ" -p gridcenter="$centerX,$centerY,$centerZ" -p spacing="$spacingVal" -o $pdbcode$ligand.gpf
	pythonsh $adtToolsPath/prepare_dpf4.py -l $ligandRandom.pdbqt -r $protein_h.pdbqt -p ga_num_evals=$gaEvals -o $pdbcode$ligand.dpf

	autogrid4 -p $pdbcode$ligand.gpf -l $pdbcode$ligand.glg
	autodock4 -p $pdbcode$ligand.dpf -l $pdbcode$ligand.dlg

	pythonsh $adtToolsPath/write_conformations_from_dlg.py -d $pdbcode$ligand.dlg
	mkdir poses

	touch rmsd.dat
	
	for pose in `seq 1 10`
	do
	    poseName=`echo $ligandRandom"_"$pose`
	    mv $poseName.pdbqt poses/.

	    cd poses

	    pythonsh $adtToolsPath/pdbqt_to_pdb.py -f $poseName.pdbqt
	    poseRed=`echo $poseName"_H"`
	    reduce $poseName.pdb > $poseRed.pdb
	    
	    rmsd=`pymol ../$ligand_h.pdb $poseRed.pdb -c -d 'align '$ligand_h', '$poseName'' | grep 'RMSD = ' | awk '{print $4}'`

	    echo $ligand" pose "$pose": "$rmsd >> ../rmsd.dat
	    
	    cd ..
	    
	done
	
	cd ..
    done
        
    cd ..
    
done
