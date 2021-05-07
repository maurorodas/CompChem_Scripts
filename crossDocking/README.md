# reDock: script for re-Docking calculation
**reDock** is a bash script to perform a simple re-Docking calculation using several crystal structures of the same protein co-cristalized with the same ligand, read the results and graph them. This script is inspired by the tutorial about re-Docking and cross-Docking published by [Bevan & Brown Lab, Public](https://osf.io/82n73/), you can access the tuturial [here](https://osf.io/82n73/wiki/Docking%2C%20Re-docking%2C%20and%20Cross%20Docking/). However, the script was written so that it can be used with different proteins and ligands than those used in the tutorial.

To run the script you only need a list of the proteins and chains where the active site and the ligand of interest are found. This list can be easily obtained from a database such as the **[PDB](https://www.rcsb.org/)**, in which using the advanced search tool, a list of proteins crystallized with the same ligand can be obtained. See **Usage** section for more details.

At this moment the script performs the docking calculation using Lamarckian Genetic Algortihm and Autodock4, you can change this modifing the script.

## Prerequisites

**reDock** is a script that uses several third party software, you need have installed these packages:

1. **[Pymol](https://pymol.org/2/)** to obtain the proteins, align and calculate RMSD
2. **[Reduce](https://ambermd.org/AmberTools.php)** to add hydrogens, **Reduce** is part of **[Ambertools 21](https://ambermd.org/AmberTools.php)**, if you don't want to use Ambertools, you can also get a copy of reduce from its **[official page](http://kinemage.biochem.duke.edu/software/reduce.php)**, you can also use another application to add the hydrogens to the system, but for this you will have to modify the script, feel free to do it
3. **[Autodock4](http://autodock.scripps.edu/)** to performs the docking calculation
4. **[MGLTools](https://ccsb.scripps.edu/mgltools/downloads/)** to prepare all files and to obtain the ligand poses
4. **[GNUPlot](http://www.gnuplot.info/)** to graph the results

## Installation
This is a bash script, you only need to copy the script in your working directory, give it the appropriate execution permissions and run in your bash shell

```bash
chmod +x reDock.sh
./reDock.sh
```

## Usage
**reDock** script needs two files for doing the job, these files must be in the same directory as the script.

The first file need to be named as proteins.dat and contains the pdb codes and chain letter of the proteins, this is an example:

```bash {data-filename="proteins.dat"}
2PRG A
3CS8 A
3DZY D
1FM6 D
1ZGY A
```

**Important Note: *** The first protein is always taken as the reference protein and you need to obtain the box parameters for that protein.

The second file nedd to be named parameter.inp, This file contains the necessary parameters to perform the docking calculation, this is an example

```bash {data-filename="proteins.dat"}
proteinNumber 5
ligand3LetterName BRL
center 50.395 -37.797 19.013
size 41 41 34
spacing 0.375
gaEvals 2500000
adtPath /PathTo/mgltools/MGLToolsPckgs/AutoDockTools/Utilities24
```

The keywords are self-describing, but below is a more detailed description:

- ```proteinNumber``` is the number of proteins in the proteins.dat file
- ```ligand3LetterName``` is three letter code of the common ligand for the proteins
- ```center``` is the xyz center position of the box
- ```size``` is the size of the box
- ```spacing``` by default ADT use and spacing of 0.375, but if you obtain you box using another program you need to change the spacing to the appropriate value, for example, Chimera use a spacing value of 1.0000
- ```gaEvals``` is the number of Genetic Algortihm evaluations
- ```adtPath``` is the path to the Utilities of the AutoDockTools

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.

## License
[GNU General Public License v. 3](https://www.gnu.org/licenses/gpl-3.0.html)

Feel free to use, share and modify it, and if you wish, let me know how it went, or how it can be improved, I still have a lot to learn. greetings from Colombia. [mauricio.rodas@ucaldas.edu.co](mauricio.rodas@ucaldas.edu.co)