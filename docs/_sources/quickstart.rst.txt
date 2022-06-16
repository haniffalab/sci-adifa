.. _quickstart:

**********
Quickstart
**********

Follow these instructions to visualise your ``.h5ad`` files in a web portal using the Adfa application.

Prerequisites
-------------

Before you begin, make sure you have installed the following tools:

- `docker`_ (>=20.10)
- `docker-compose`_ (>=1.27)

.. _docker: https://docs.docker.com/get-docker/
.. _docker-compose: https://docs.docker.com/compose/install/

Create a docker-compose file
----------------------------

Create an empty file called ``docker-compose.yml``, and copy the following code into file.

:: 

    version: '2'
    services:
        adifa:
            image: haniffalab/adifa:latest
            container_name: adifa
            entrypoint: sh
            command: ./boot-autodiscover.sh
            ports:
            - 5000:5000
            volumes:
            - ./path/to/your/h5ad:/data
            environment:
            - "DATA_PATH=/data" 

Replace ``/path/to/your/h5ad`` with the path to the directory where your ``.h5ad`` files are located. Save the file.

Start the docker service
------------------------

Run the following command in the same directory:

:: 

    docker-compose up

This command will load the Adifa application, read your ``.h5ad`` files and start the web server. Depending on the size of your data files, the application may take several minutes to load. 

To access the application, open http://localhost:5000 in a browser. 

.. attention::

    To learn more about using Adifa with Docker, including deployment strategies and guidance, visit :ref:`deployment <deployment>`



