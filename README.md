#rmsd
Calculates the rmsd between two geometries given as xyz file.

Description

This program calculates the RMSD-value (in Angstrom) between two structures given as xyz file. The result is printed to screen. In addition, a xyz file with the aligned geometry of the structure given first is written to file.

Usage

rmsd example.xyz example_two.xyz

rmsd example.xyz example_two.xyz -q
rmsd example.xyz example_two.xyz --quaternion

Tested with ORCA 5.0. Precompiled binaries for Ubuntu and Windows 10 are statically linked and do not have further dependencies. When compiling the source code from scratch, Third-Party libraries (Boost and fmt) are needed. Boost.program_options is not header-only and therefore requires a proper installation of Boost. Note: fmt was used in header_only mode so the FMT_HEADER_ONLY macro needs to be set.
Contributor

contributed by Ingolf Harden
