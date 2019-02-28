polymer_radiolysis_pub
======================
.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.2579871.svg
   :target: https://doi.org/10.5281/zenodo.2579871
   :alt: Zenodo DOI
`polymer_radiolysis_pub <https://github.com/bjodah/polymer_radiolysis_pub>` contains the source code for our model
of radioylsis of aqueous polymer solutions. An archived copy of the Dockerimage is published on Zenodo for reproducibility.

Usage
-----
To compile the source code and run a simulation you may run:

.. code:: bash

   $ curl -LOs https://zenodo.org/record/2579871/files/dockerimage__polymer_radiolysis_pub.tar.xz
   $ openssl sha256 dockerimage__polymer_radiolysis_pub.tar.xz
   SHA256(dockerimage__polymer_radiolysis_pub.tar.xz)= f3560be844775ce2d0e2ffe1a9ca52109a4bf1e3f88b9058f5ab72f3e75d0266
   $ xz --decompress - < dockerimage__polymer_radiolysis_pub.tar.xz | docker image load
   $ ./scripts/run_dockerized_build_and_test.sh

this will compile the source code under ``src/`` and run a small test-simulation.

Licenseing
----------
This code (except files under "external/") is licensed under the GNU General Public
License as published by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version (see file LICESENE for further information).

The files under "external/" have their own licenses which are distributed therein.

Author
------
Björn Ingvar Dahlgren, contact:
 - gmail address: bjodah
 - kth.se address: bda
