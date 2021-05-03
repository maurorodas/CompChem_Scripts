#!/bin/bash

# file with proteins and chain to perform a redocking calculation
filename='proteins.dat'
# Name in the pdb file of the ligand to use for the redocking
ligandName='BRL'
# Number of the proteins to perform the redocking calculation
protNumber=5

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
	pymol -c -d 'fetch '$my_query', async=0; save '$pdbcode'.pdb'
	grep $ligandName $pdbcode.pdb > $ligand.pdb
	grep -v -e HOH -e $ligandName -e CONECT $pdbcode.pdb > temp_$pdbcode.pdb
	mv temp_$pdbcode.pdb $pdbcode.pdb
	cd ..
    else
	# Performs the aligment using the "align" method of Pymol, saves the protein with the new coordinates,
	# extract the ligand and make a basic cleaning
	cd $fileProt
	pymol -c -d 'fetch '$my_query' '$query_ref', async=0; align '$my_query', '$query_ref'; save '$pdbcode'.pdb, '$my_query''
	grep $ligandName $pdbcode.pdb > $ligand.pdb
	grep -v -e HOH -e $ligandName -e CONECT $pdbcode.pdb > temp_$pdbcode.pdb
	mv temp_$pdbcode.pdb $pdbcode.pdb
	cd ..
    fi      
done
