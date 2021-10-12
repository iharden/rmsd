# rmsd
Calculates the RMSD between two geometries given as xyz file.

## Description

This command line tool calculates the RMSD-value (in Angstrom) between two structures given as .xyz files. The result is printed to screen. In addition, the coordinates of the aligned geometry are written to file.
The default algorithm for calculating the transformation matrix for molecule alignment is the Kabsch algorithm (https://en.wikipedia.org/wiki/Kabsch_algorithm).
The quaterion-based algorithm (https://doi.org/10.1002/jcc.20110) is implemented as well and can be invoked from the command line. The program is written in C++ and uses OpenMP for parallelization.

## Usage

`rmsd example.xyz example_two.xyz`

`rmsd example.xyz example_two.xyz -q`

`rmsd example.xyz example_two.xyz --quaternion`

`-q` and `--quaternion` are equivalent and trigger the quaternion algorithm for calculating the transformation matrix.

## Requirements

Precompiled binaries for Ubuntu and Windows 10 are avaible and can be used right away. When compiling the source code from scratch, the third-party Eigen package (https://eigen.tuxfamily.org/index.php?title=Main_Page) needs to be downloaded. Eigen is header-only and therefore should be easy to include.

## Contributor

contributed by Ingolf Harden
