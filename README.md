# Quantum-Espresso-codes

## Some important conversions

1 Angstorm = 1.88971616463 Bohr		\
1 Bohr = 0.529177249 Angstorm = 1 a.u	\
1 Ry = 13.6046981 eV			\
1 eV = 1.60217733 x 10^-19 J		

> https://www.quantum-espresso.org/Doc/pw_user_guide/			\
> https://www.quantum-espresso.org/Doc/INPUT_PW.html			\
> https://www.materialscloud.org/work/tools/qeinputgenerator

## Finding 'scf' for spin-up in LiCab

### What is 'scf'?

In Quantum ESPRESSO, a self-consistent field (SCF) calculation is a basic operation that minimizes a system's energy for a fixed geometry. SCF calculations are an iterative method that involves:        
* Selecting an approximate Hamiltonian        
* Solving the Schrödinger equation to get a more accurate set of orbitals        
* Solving the Schrödinger equation again with the new orbitals until the results converge        

SCF calculations can be performed for varying ionic positions, and the ground state energies can be plotted against those positions. A geometry relaxation is a series of SCF calculations, where the ionic position of the relaxed structure is given.
To perform an SCF calculation for silicon using Quantum ESPRESSO, you need to provide various parameters in an input file. The default value for the calculation parameter is "scf". 
The PWscf (Plane-Wave Self-Consistent Field) package is a core component of the Quantum ESPRESSO distribution that performs many different kinds of self-consistent calculations. These calculations include: 
* Ground-state energy 
* One-electron (Kohn-Sham) orbitals
* Atomic forces
* Stresses
* Polarization
* Effective Screening Medium (ESM) method

```
&control        
	calculation = 'scf',    
	prefix = 'LiCaB'        
	outdir = '/home/lopzz/....../temp'        
	pseudo_dir = '/home/......./pseudo'        
/        
&system        
	ibrav = 2			#ibrav: Bravais lattice index        
	celldm(1) = 11.3553064        
	nat = 3,			#nat: no. of atoms        
	ntyp = 3,			#ntyp: no. of types of atoms        
	ecutwfc = 15			#ecutwfc: K.E, cutoff (Ry) for wavefunction        
	ecutrho = 120			#ecutrho: K.E, cutoff (Ry) for charge density & potential        
 	ocupations = 'smearing'
	degauss = 0.002        
	nspin = 1			#this is particularly for up-spin, for down spin the value will be 2        
	starting_magnetization(1) = 0.5        
	starting_magnetization(2) = 0.5        
	starting_magnetization(3) = 0.5        
/        
&electrons        
	mixing_beta = 0.7		#mixing_beta: mixing factor for self-consistency        
	diagonalization = 'cg'		#cg: conjugate gradient-like band-by-band diagonalization        
/        
&ions    
/        
ATOMIC_SPECIES        
Li 	6.941 Li.pbe-s-rrkjus_psl.0.2.1.UPF      
Ca 	40.078 Ca.pbe-spn-rrkjus_psl.0.2.3.UPF   
B 	10.811 B.pbe-n-rrkjus_psl.0.1.UPF        
ATOMIC_POSITIONS(alat)			#alat: atomic position are in cartesian coordinates in units of cell parameters        
B	0.25 	0.25	0.25        
Ca	0.00	0.00	0.00        
Li	0.50 	0.50	0.50        
k-POINTS (automatic)			#automatic: automatically generated uniform grid of k-points        
6	6	6	1	1	1        
```

## Calculation of cut off energy

### What is cut-off energy?
In Quantum Espresso, the cutoff energy limits the number of plane wave functions used as basis functions to represent the wavefunction. In theory, an infinite number of basis functions would be needed for an exact answer, but this is not computationally feasible.		\
Here are some tips for choosing the cutoff energy in Quantum Espresso:	
* For norm-conserving pseudopotentials, use the default value, or reduce it slightly. Reducing it too much can introduce noise, especially for forces and stress.
* For ultrasoft pseudopotentials, use a larger value than the default, typically 8–12 times ecutwfc.
* For PAW datasets, use 4*ecutwfc, but test first to make sure it works, as it depends on the shape of the augmentation charge. 
You can perform self-consistent calculations with different cutoff energies and sets of K-points to get the right values for your system. You can also run a convergence test to find a setting that's optimized for a model. 

```
#1/bin/sh					
NAME= 'ecut_LiCaB'				
for ecut_LiCaB in 5 10 15 20 25 30 35 40 	
do						
cat > $NAME.$ecut_LiCaB.in << EOF		
&control					
	calculation = 'scf'			
 	prefix = 'LiCaB_ecut'			
  	outdir = '/home/lopazz/Desktop/MSc_Project_2022/temp'
   	pseudo_dir = '/home/lopazz/quantum_espresso/pseudo'	
/								
&system								
	ibrav = 2,						
 	celldm(1) = 11.33,					
  	nat = 3,							
   	ntyp = 3,						
    	ecutwfc = $ecut_LiCab					
     	ecutrho = 320						
/								
&electrons							
	mixing_beta = 0.7					
 	diagonalization = 'cg'					
/
ATOMIC_SPECIES	
Li 	6.941 Li.pbe-s-rrkjus_psl.0.2.1.UPF        
Ca 	40.078 Ca.pbe-spn-rrkjus_psl.0.2.3.UPF     
B 	10.811 B.pbe-n-rrkjus_psl.0.1.UPF        
ATOMIC_POSITIONS(alat)				
B	0.25 	0.25	0.25        
Ca	0.00	0.00	0.00        
Li	0.50 	0.50	0.50        
k-POINTS(automatic)		
10	10	10	1	1	1		
EOF							
~/espresso-5.2.0/bin/pw.x<$NAME.$ecut_LiCaB.in>$NAME.$ecut_LiCaB.out	
done									
```

## k-Points

```
#1/bin/sh					
NAME= 'LiCaB_nk'				
for LiCaB_nk in 4 6 8 10 12 14 16 18
do						
cat > $NAME.$LiCaB_nk.in << EOF		
&control					
	calculation = 'scf'			
 	prefix = 'LiCaB_nk'			
  	outdir = '/home/lopazz/Desktop/MSc_Project_2022/temp'
   	pseudo_dir = '/home/lopazz/quantum_espresso/pseudo'	
/								
&system								
	ibrav = 2,						
 	celldm(1) = 11.33,					
  	nat = 3,							
   	ntyp = 3,						
    	ecutwfc = 15					
     	ecutrho = 120						
/								
&electrons							
	mixing_beta = 0.7					
 	diagonalization = 'cg'					
/
ATOMIC_SPECIES 	
Li 	6.941 Li.pbe-s-rrkjus_psl.0.2.1.UPF        
Ca 	40.078 Ca.pbe-spn-rrkjus_psl.0.2.3.UPF     
B 	10.811 B.pbe-n-rrkjus_psl.0.1.UPF        
ATOMIC_POSITIONS(alat)				
B	0.25 	0.25	0.25        
Ca	0.00	0.00	0.00        
Li	0.50 	0.50	0.50        
k-POINTS(automatic)		
$LiCab_nk	$LiCaB_nk	$LiCaB	1	1	1		
EOF							
~/espresso-5.2.0/bin/pw.x<$NAME.$LiCaB_nk.in>$NAME.$LiCaB_nk.out	
done									
```
Run the file 
```
./LiCaB_nk.sh &
```

## Calculation of lattice constant

```
#1/bin/sh					
NAME= 'LiCaB_alat'				
for LiCaB_alat in 11.0 11.2 11.4 11.6 11.8 12.0 	
do						
cat > $NAME.$LiCaB_alat.in << EOF		
&control					
	calculation = 'scf'			
 	prefix = 'LiCaB_alat'			
  	outdir = '/home/lopazz/Desktop/MSc_Project_2022/temp'
   	pseudo_dir = '/home/lopazz/quantum_espresso/pseudo'	
/								
&system								
	ibrav = 2,						
 	celldm(1) = $LiCaB_alat,					
  	nat = 3,							
   	ntyp = 3,						
    	ecutwfc = 15				
     	ecutrho = 120						
/								
&electrons							
	mixing_beta = 0.7					
 	diagonalization = 'cg'					
/
ATOMIC_SPECIES	
Li 	6.941 Li.pbe-s-rrkjus_psl.0.2.1.UPF        
Ca 	40.078 Ca.pbe-spn-rrkjus_psl.0.2.3.UPF     
B 	10.811 B.pbe-n-rrkjus_psl.0.1.UPF        
ATOMIC_POSITIONS(alat)				
B	0.25 	0.25	0.25        
Ca	0.00	0.00	0.00        
Li	0.50 	0.50	0.50        
k-POINTS(automatic)		
10	10	10	1	1	1		
EOF							
~/espresso-5.2.0/bin/pw.x<$NAME.$LiCaB_alat.in>$NAME.$LiCaB_alat.out	
done
```
Run the file 
```
./LiCaB_alat*.out > LiCaB_alat.dat
```
```
grep! LiCaB_alat.dat
```
```
vi LiCaB_alat.dat
```
```
~/espresso-5.2.0/bin/ev.x
```
The value of lattice constant = 6.3499 Å

## Magnetization and spin-orbit interaction

```
&control        
	calculation = 'scf',    
	prefix = 'LiCaB'        
	outdir = '/home/lopzz/....../temp'        
	pseudo_dir = '/home/......./pseudo'        
/        
&system        
	ibrav = 2			#ibrav: Bravais lattice index        
	celldm(1) = 11.3553064        
	nat = 3,			#nat: no. of atoms        
	ntyp = 3,			#ntyp: no. of types of atoms        
	ecutwfc = 15			#ecutwfc: K.E, cutoff (Ry) for wavefunction        
	ecutrho = 120			#ecutrho: K.E, cutoff (Ry) for charge density & potential        
 	ocupations = 'smearing'
	degauss = 0.002        
	nspin = 1			#this is particularly for up-spin, for down spin the value will be 2        
	starting_magnetization(1) = 0.5        
	starting_magnetization(2) = 0.5        
	starting_magnetization(3) = 0.5        
/        
&electrons        
	mixing_beta = 0.7		#mixing_beta: mixing factor for self-consistency        
	diagonalization = 'cg'		#cg: conjugate gradient-like band-by-band diagonalization        
/        
&ions    
/        
ATOMIC_SPECIES        
Li 	6.941 Li.pbe-s-rrkjus_psl.0.2.1.UPF      
Ca 	40.078 Ca.pbe-spn-rrkjus_psl.0.2.3.UPF   
B 	10.811 B.pbe-n-rrkjus_psl.0.1.UPF        
ATOMIC_POSITIONS(alat)			#alat: atomic position are in cartesian coordinates in units of cell parameters        
B	0.25 	0.25	0.25        
Ca	0.00	0.00	0.00        
Li	0.50 	0.50	0.50        
k-POINTS (automatic)			#automatic: automatically generated uniform grid of k-points        
6	6	6	1	1	1        
```

## Dielectric constant

intersmear - It is the broadening for interband transitions. If it is too small, a number of spurious peaks due to numerical reasons will show up. It is a broadening parameter (in eV). It will be Gaussian or Lorentzian broadening parameter in the case of joint density state calculation or Drude-Lorentz broadening parameter for the dielectric tensor calculation.

intrasmear - It is the broadening parameter for the intraband i.e, metal Drude- like term (in eV)

> https://pranabdas.github.io/espresso/assets/files/eps_man-ab3fac19eb366509dd129c37fbf94ac0.pdf

### Interband

INPUT FILE

OUTPUT FILE


### Intraband


INPUT FILE
```
&inputpp
    calculation = 'eps',
    prefix = 'LiCaB_eps'
    outdir = '/home/lopazz/Desktop/MSc_Project_2022/temp'
 /
 &energy_grid
    smeartype = 'gauss'
    intersmear = 0.0d0
    intrasmear = 0.136d0
    wmax = 30.0d0
    wmin = 0.0d0
    nw = 600
    shift = 0.0d0
 /
```

OUTPUT FILE

```

     Program epsilon v.5.2.0 starts on 23Jun2022 at 18:48:14 

     This program is part of the open-source Quantum ESPRESSO suite
     for quantum simulation of materials; please cite
         "P. Giannozzi et al., J. Phys.:Condens. Matter 21 395502 (2009);
          URL http://www.quantum-espresso.org", 
     in publications or presentations arising from this work. More details at
     http://www.quantum-espresso.org/quote

     Serial version


     Reading input file...
     Broadcasting variables...
     Reading PW restart file...

     Reading data from directory:
     /home/bulumoni/Desktop/Msc_4th_sem_Project_2022/Hirashri/conv/temp/LiCaB_oncv_optics.save

   Info: using nr1, nr2, nr3 values from input

   Info: using nr1, nr2, nr3 values from input

     IMPORTANT: XC functional enforced from input :
     Exchange-correlation      = PBE ( 1  4  3  4 0 0)
     Any further DFT definition will be discarded
     Please, verify this is what you really want


     G-vector sticks info
     --------------------
     sticks:   dense  smooth     PW     G-vecs:    dense   smooth      PW
     Sum        1147    1147    337                25605    25605    4015

     Generating pointlists ...
     new r_m :   0.1786 (alat units)  2.1031 (a.u.) for type    1
     new r_m :   0.1786 (alat units)  2.1031 (a.u.) for type    2
     new r_m :   0.1786 (alat units)  2.1031 (a.u.) for type    3


     Fermi energy [eV] is:  4.48252
     The system is a metal...

     Performing eps calculation...

     Writing output on file...


     epsilon      :     1.23s CPU         1.26s WALL

     calculation  :      0.38s CPU      0.38s WALL (       1 calls)
     dipole_calc  :      0.28s CPU      0.28s WALL (      70 calls)


stop_clock: clock #  1 for      epsilon not running
     epsilon      :     1.23s CPU         1.26s WALL


   This run was terminated on:  18:48:15  23Jun2022            

=------------------------------------------------------------------------------=
   JOB DONE.
=------------------------------------------------------------------------------=

```

### Total or inter+intraband


INPUT FILE
```
&inputpp
    calculation = 'eps',
    prefix = 'LiCaB_eps'
    outdir = '/home/bulumoni/Desktop/Msc_4th_sem_Project_2022/Hirashri/conv/temp/'
 /
 &energy_grid
    smeartype = 'gauss'
    intersmear = 0.136d0
    intrasmear = 0.136d0
    wmax = 30.0d0
    wmin = 0.0d0
    nw = 600
    shift = 0.0d0
 /

```
OUTUT FILE
```

     Program epsilon v.5.2.0 starts on 23Jun2022 at 18:41:29 

     This program is part of the open-source Quantum ESPRESSO suite
     for quantum simulation of materials; please cite
         "P. Giannozzi et al., J. Phys.:Condens. Matter 21 395502 (2009);
          URL http://www.quantum-espresso.org", 
     in publications or presentations arising from this work. More details at
     http://www.quantum-espresso.org/quote

     Serial version


     Reading input file...
     Broadcasting variables...
     Reading PW restart file...

     Reading data from directory:
     /home/bulumoni/Desktop/Msc_4th_sem_Project_2022/Hirashri/conv/temp/LiCaB_oncv_optics.save

   Info: using nr1, nr2, nr3 values from input

   Info: using nr1, nr2, nr3 values from input

     IMPORTANT: XC functional enforced from input :
     Exchange-correlation      = PBE ( 1  4  3  4 0 0)
     Any further DFT definition will be discarded
     Please, verify this is what you really want


     G-vector sticks info
     --------------------
     sticks:   dense  smooth     PW     G-vecs:    dense   smooth      PW
     Sum        1147    1147    337                25605    25605    4015

     Generating pointlists ...
     new r_m :   0.1786 (alat units)  2.1031 (a.u.) for type    1
     new r_m :   0.1786 (alat units)  2.1031 (a.u.) for type    2
     new r_m :   0.1786 (alat units)  2.1031 (a.u.) for type    3


     Fermi energy [eV] is:  4.48252
     The system is a metal...

     Performing eps calculation...

     Writing output on file...


     epsilon      :     1.24s CPU         1.30s WALL

     calculation  :      0.39s CPU      0.41s WALL (       1 calls)
     dipole_calc  :      0.28s CPU      0.28s WALL (      70 calls)


stop_clock: clock #  1 for      epsilon not running
     epsilon      :     1.24s CPU         1.30s WALL


   This run was terminated on:  18:41:30  23Jun2022            

=------------------------------------------------------------------------------=
   JOB DONE.
=------------------------------------------------------------------------------=
```
