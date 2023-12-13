# Dihedral-angle-scan

Dihedral angle is one crucial parameter in constructing potential energy surface 
1. Relaxed Scan
   Scan a series of different dihedral angle, relax the structure before getting the energy
2. Rigid Scan
   Scan a series of different dihedral angle without optimization structure.

Tutorial
1. First, optimize the molecule geometry. Using Chemdraw to create the .mol file (containing chemical structure and preliminary geometry), then using Avogadro to transfer into .xyz file.
2. Use Gaussian to optimize structure.

3. Using create_dihedral.py to create different dihedral angle structures.

4. Next upload all structures to NERSC or cluster with Gaussian. Run the command: for i in {0..180..5}; do mkdir $i; cd $i; mv ../${i}degree.xyz .; cd ..; done.

for i in {0..180..5}; do cd $i; bash createscan.sh ${i}degree.xyz $i; cd ..; done

for i in {0..180..5}; do cd $i; sbatch run_${i}; cd ..; done
