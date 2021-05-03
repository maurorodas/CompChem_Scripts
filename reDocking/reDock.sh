#!/bin/bash

filename='proteins.dat'

for line in `seq 1 5`
do
    pdbcode=`awk -v var=$line '{if(NR==var) print $1}' $filename`
    chain=`awk -v var=$line '{if(NR==var) print $2}' $filename`
    my_query=`echo "$pdbcode$chain"`

    fileProt=`echo "proteina$line"`
    ligand=`echo "ligando$line"`

    if [ -d $fileProt ]
    then
	rm -r $fileProt
	mkdir $fileProt
    else
	mkdir $fileProt
    fi

    if [ $line -eq 1 ]
    then
	protRef=$pdbcode
	chainRef=$chain
	query_ref=`echo "$protRef$chainRef"`
    fi
      
done
