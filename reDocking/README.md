# Bash script for Redocking calculation
**reDock** is a bash script to perform a redocking calculation using several proteins and ligands, read the results and graph them.

For now, the script fetch the proteins and align them.

## Prerequisites

**reDock** is a script that uses several third party softwares, you need have installed these packages:

1. Pymol
2. Autodock
3. GNUPlot

## Instalation
Only copy the script in a directory of your choice and run in your bash shell

```bash
./reDock.sh
```

## Usage
**reDock** script needs a file named proteins.dat with the pdb codes and chain letter of the proteins, this is an example:
```bash
2PRG A
3CS8 A
3DZY D
1FM6 D
1ZGY A
```

The first protein is always taken as the reference protein

## License
[GNU General Public License v. 3](https://www.gnu.org/licenses/gpl-3.0.html)

Feel free to use, share and modify it, and if you wish, let me know how it went, or how it can be improved, I still have a lot to learn. greetings from Colombia. [mauricio.rodas@ucaldas.edu.co](mauricio.rodas@ucaldas.edu.co)