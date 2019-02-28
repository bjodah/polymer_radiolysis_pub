The dockerimage was built using `Docker <https://docs.docker.com/install/linux/docker-ce/ubuntu/>`_ version 17.05.0-ce, build 89658be:

.. code:: bash

   $ docker build -t polymer_radiolysis_pub ./environment/

and compressed using pxz:

.. code:: bash

   $ docker image save polymer_radiolysis_pub | pxz -cv > dockerimage__polymer_radiolysis_pub.tar.xz

and the compressed image uploaded to zenodo.
