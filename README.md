# rmsd
Calculates the RMSD between two geometries given as xyz file.

## Description

This command line tool calculates the RMSD-value between two structures given as .xyz files. The result is printed to screen. In addition, the coordinates of the aligned geometry are written to file.
The default algorithm for calculating the transformation matrix for molecule alignment is the Kabsch algorithm (https://en.wikipedia.org/wiki/Kabsch_algorithm).
The quaternion-based algorithm (https://doi.org/10.1002/jcc.20110) is implemented as well and can be invoked from the command line. The program is written in C++ and uses OpenMP for parallelization.

## Usage

`rmsd example.xyz example_two.xyz`

`rmsd example.xyz example_two.xyz -q`

`rmsd example.xyz example_two.xyz --quaternion`

`-q` and `--quaternion` are equivalent and trigger the quaternion algorithm for calculating the transformation matrix.

## Requirements

Precompiled binaries for Ubuntu and Windows 10 are available and can be used right away. When compiling the source code from scratch, the third-party Eigen package (https://eigen.tuxfamily.org/index.php?title=Main_Page) needs to be downloaded. Eigen is header-only and therefore should be easy to include.

## Contributor

contributed by Ingolf Harden

## License

This program is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but without any warranty; without even the implied warranty of merchantability or fitness for a particular purpose. See the GNU Lesser General Public License for more details.

Unless you explicitly state otherwise, any contribution intentionally submitted for inclusion in this program by you, as defined in the GNU Lesser General Public license, shall be licensed as above, without any additional terms or conditions
