#!/bin/bash

paramsFile='parameters.inp'
# Number of the proteins to perform the redocking calculation
protNumber=`awk '{if(NR==1) print $2}' $paramsFile`
# Name in the pdb file of the ligand to use for the redocking
ligandName=`awk '{if(NR==2) print $2}' $paramsFile`

centerX=`awk '{if(NR==3) print $2}' $paramsFile`
centerY=`awk '{if(NR==3) print $3}' $paramsFile`
centerZ=`awk '{if(NR==3) print $4}' $paramsFile`

sizeX=`awk '{if(NR==4) print $2}' $paramsFile`
sizeY=`awk '{if(NR==4) print $3}' $paramsFile`
sizeZ=`awk '{if(NR==4) print $4}' $paramsFile`

spacingVal=`awk '{if(NR==5) print $2}' $paramsFile`
gaEvals=`awk '{if(NR==6) print $2}' $paramsFile`

adtToolsPath=`awk '{if(NR==7) print $2}' $paramsFile`

# file with proteins and chain to perform a redocking calculation
filename='proteins.dat'
workindDir=`pwd`

# Directory to save all ligands 
if [ -d ligandsDir ]
then
    rm -r ligandsDir
    mkdir ligandsDir
else
    mkdir ligandsDir
fi

if [ -d Results ]
then
    rm -r Results
    mkdir Results
else
    mkdir Results
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
    
    dataFile=`echo "rmsd_"$pdbcode".dat"`

    cp ligandsDir/ligand* $fileProt/.

    cd $fileProt
    pythonsh $adtToolsPath/prepare_receptor4.py -r $protein_h.pdb -A hydrogens -o $protein_h.pdbqt

    ligand=''
    for lig in `seq 1 $protNumber`
    do
	ligand=`echo "ligand$lig"`
	vectorLigs[$lig-1]="Lig"$lig
	
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

	

	for pose in `seq 1 10`
	do
	    poseName=`echo $ligandRandom"_"$pose`
	    mv $poseName.pdbqt poses/.

	    cd poses

	    realPosition=`grep "RANKING" ../$pdbcode$ligand.dlg | awk -v var="$pose" '{ if ( $3 == var ) print NR }'`
	    # bindingEnergy=`grep "RANKING" ../$pdbcode$ligand.dlg | awk -v var="$pose" '{ if ( $3 == var ) print $4 }'`

	    pythonsh $adtToolsPath/pdbqt_to_pdb.py -f $poseName.pdbqt

	    poseOrder=`echo "pose"$realPosition`
	    
	    rm $poseName.pdbqt
	    mv $poseName.pdb $poseOrder.pdb
	    
	    poseRed=`echo $poseOrder"_H"`
	    reduce $poseOrder.pdb > $poseRed.pdb

	    cd ..
	done
	
	for pose in `seq 1 10`
	do
	    cd poses

	    poseRed=`echo "pose"$pose"_H"`
	    
	    rmsd=`pymol ../$ligand_h.pdb $poseRed.pdb -c -d 'align '$ligand_h', '$poseRed'' | grep 'RMSD = ' | awk '{print $4}'`

	    if [ $lig -eq 1 ]
	    then
		echo $pose" "$rmsd >> $workindDir/Results/$dataFile
		echo $pdbcode" "$rmsd >> $workindDir/Results/allRMSDdata.dat
	    else
		sed -i.bak "${pose}s/$/ ${rmsd}/" $workindDir/Results/$dataFile
		echo $pdbcode" "$rmsd >> $workindDir/Results/allRMSDdata.dat
		rm $workindDir/Results/*.bak
	    fi
	    	    
	    cd ..
	    
	done
	
	cd ..
    done

    gnuPlotFile="gPlot"$pdbcode".gp"
    gnuPlotGraph="graph_"$pdbcode".png"
    
    touch $workindDir/Results/$gnuPlotFile

    echo "#!/usr/bin/gnuplot" >> $workindDir/Results/$gnuPlotFile
    echo "" >> $workindDir/Results/$gnuPlotFile
    echo 'set terminal pngcairo enhanced background "#ffffff" fontscale 2.5 dashed size 1920, 1280' >> $workindDir/Results/$gnuPlotFile
    echo "" >> $workindDir/Results/$gnuPlotFile
    echo "set encoding iso_8859_1" >> $workindDir/Results/$gnuPlotFile
    echo "set output '"$workindDir"/Results/"$gnuPlotGraph"'" >> $workindDir/Results/$gnuPlotFile
    echo "" >> $workindDir/Results/$gnuPlotFile
    echo "set xrange [1:10]" >> $workindDir/Results/$gnuPlotFile
    echo "set yrange [0:*]" >> $workindDir/Results/$gnuPlotFile
    echo "" >> $workindDir/Results/$gnuPlotFile
    echo 'set xlabel "Poses" font "Arial, 20"' >> $workindDir/Results/$gnuPlotFile
    echo 'set ylabel "RMSD"  font "Arial, 20"' >> $workindDir/Results/$gnuPlotFile
    echo "" >> $workindDir/Results/$gnuPlotFile
    echo 'set key outside top horizontal font "Arial, 18" maxcols 4' >> $workindDir/Results/$gnuPlotFile
    echo 'set xtics axis nomirror out font "Arial, 20"' >> $workindDir/Results/$gnuPlotFile
    echo 'set ytics axis nomirror out font "Arial, 20"' >> $workindDir/Results/$gnuPlotFile
    echo "set mxtics" >> $workindDir/Results/$gnuPlotFile
    echo "set mytics" >> $workindDir/Results/$gnuPlotFile
    echo "" >> $workindDir/Results/$gnuPlotFile
    echo 'set arrow from 1,2 to 10,2 nohead dt 9 lw 4 lc "red"' >> $workindDir/Results/$gnuPlotFile
    echo 'list = "'${vectorLigs[@]}'"' >> $workindDir/Results/$gnuPlotFile
    echo "item(n) = word(list,n)" >> $workindDir/Results/$gnuPlotFile
    echo 'plot for [i=1:words(list)] "'$workindDir/Results/$dataFile'" using 1:i+1 title item(i) with linespoints lw 2.5 ps 3' >> $workindDir/Results/$gnuPlotFile

    gnuplot $workindDir/Results/$gnuPlotFile
       
    cd ..
    
done

#Sorting RMSD values and extracting the lowest value for each protein and graph them

cd Results

for i in `seq 1 $protNumber`
do
    pdbcode=`awk -v var=$i '{if(NR==var) print $1}' ../$filename`
    dataFile=`echo "rmsd_"$pdbcode".dat"`
    dataLowest=`echo "allProteinsLowestRMSD.dat"`

    vectorProts[$i-1]=$pdbcode

    for j in `seq 1 $protNumber`
    do
	column=$((j + 1))
	lowestRMSD=`sort -n -k $column -t ' ' $dataFile | awk -v var=$column '{if(NR==1) print $var}'`
	
	if [ $i -eq 1 ]
	then
	    echo $j" "$lowestRMSD >> $workindDir/Results/$dataLowest
	else
	    sed -i.bak "${j}s/$/ ${lowestRMSD}/" $workindDir/Results/$dataLowest
	    rm $workindDir/Results/*.bak
	fi
    done
done

echo "#!/usr/bin/gnuplot" >> $workindDir/Results/gnuPlotAllRMSD.gp
echo "" >> $workindDir/Results/gnuPlotAllRMSD.gp
echo 'set terminal pngcairo enhanced background "#ffffff" fontscale 2.5 dashed size 1920, 1280' >> $workindDir/Results/gnuPlotAllRMSD.gp
echo "set encoding iso_8859_1" >> $workindDir/Results/gnuPlotAllRMSD.gp
echo "set output '"$workindDir"/Results/gnuPlotAllRMSD.png'" >> $workindDir/Results/gnuPlotAllRMSD.gp
echo "" >> $workindDir/Results/gnuPlotAllRMSD.gp
echo "set style fill solid 0.5 border -1" >> $workindDir/Results/gnuPlotAllRMSD.gp
echo "set style boxplot outliers pointtype 7" >> $workindDir/Results/gnuPlotAllRMSD.gp
echo "set style data boxplot" >> $workindDir/Results/gnuPlotAllRMSD.gp
echo "set boxwidth  0.5" >> $workindDir/Results/gnuPlotAllRMSD.gp
echo "set pointsize 3" >> $workindDir/Results/gnuPlotAllRMSD.gp
echo "" >> $workindDir/Results/gnuPlotAllRMSD.gp
echo "unset key" >> $workindDir/Results/gnuPlotAllRMSD.gp
echo 'set xtics auto nomirror font "Arial"' >> $workindDir/Results/gnuPlotAllRMSD.gp
echo 'set ytics axis nomirror out font "Arial, 20"' >> $workindDir/Results/gnuPlotAllRMSD.gp
echo "set yrange [*:*]" >> $workindDir/Results/gnuPlotAllRMSD.gp
echo "" >> $workindDir/Results/gnuPlotAllRMSD.gp
echo 'set arrow from 0,2 to '$protNumber',2 nohead dt 9 lw 4 lc "red"' >> $workindDir/Results/gnuPlotAllRMSD.gp
echo 'set ylabel "RMSD" font "Arial, 20"' >> $workindDir/Results/gnuPlotAllRMSD.gp
echo "" >> $workindDir/Results/gnuPlotAllRMSD.gp
echo "plot 'allRMSDdata.dat' using (1):2:(0):1 lc variable lw 2" >> $workindDir/Results/gnuPlotAllRMSD.gp

gnuplot $workindDir/Results/gnuPlotAllRMSD.gp

gnuPlotAllProts="allProteinsLowestRMSD.gp"
gnuPlotGraphAllprots="graph_allProteinsLowestRMSD.png"

touch $workindDir/Results/$gnuPlotAllProts

echo "#!/usr/bin/gnuplot" >> $workindDir/Results/$gnuPlotAllProts
echo "" >> $workindDir/Results/$gnuPlotAllProts
echo 'set terminal pngcairo enhanced background "#ffffff" fontscale 2.5 dashed size 1920, 1280' >> $workindDir/Results/$gnuPlotAllProts
echo "" >> $workindDir/Results/$gnuPlotAllProts
echo "set encoding iso_8859_1" >> $workindDir/Results/$gnuPlotAllProts
echo "set output '"$workindDir"/Results/"$gnuPlotGraphAllprots"'" >> $workindDir/Results/$gnuPlotAllProts
echo "" >> $workindDir/Results/$gnuPlotAllProts
echo "set xrange [1:"$protNumber"]" >> $workindDir/Results/$gnuPlotAllProts
echo "set yrange [0:*]" >> $workindDir/Results/$gnuPlotAllProts
echo "" >> $workindDir/Results/$gnuPlotAllProts
echo 'set xlabel "Ligands" font "Arial, 20"' >> $workindDir/Results/$gnuPlotAllProts
echo 'set ylabel "RMSD"  font "Arial, 20"' >> $workindDir/Results/$gnuPlotAllProts
echo "" >> $workindDir/Results/$gnuPlotAllProts
echo 'set key outside top horizontal font "Arial, 18" maxcols 4' >> $workindDir/Results/$gnuPlotAllProts
echo 'set xtics axis nomirror out font "Arial, 20"' >> $workindDir/Results/$gnuPlotAllProts
echo 'set ytics axis nomirror out font "Arial, 20"' >> $workindDir/Results/$gnuPlotAllProts
echo "set mxtics" >> $workindDir/Results/$gnuPlotAllProts
echo "set mytics" >> $workindDir/Results/$gnuPlotAllProts
echo "" >> $workindDir/Results/$gnuPlotAllProts
echo 'set arrow from 1,2 to '$protNumber',2 nohead dt 9 lw 4 lc "red"' >> $workindDir/Results/$gnuPlotAllProts
echo 'list = "'${vectorProts[@]}'"' >> $workindDir/Results/$gnuPlotAllProts
echo "item(n) = word(list,n)" >> $workindDir/Results/$gnuPlotAllProts
echo 'plot for [i=1:words(list)] "'$workindDir/Results/$dataLowest'" using 1:i+1 title item(i) with linespoints lw 2.5 ps 3' >> $workindDir/Results/$gnuPlotAllProts

gnuplot $workindDir/Results/$gnuPlotAllProts

cd ..
