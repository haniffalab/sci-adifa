.. _deployment:

**********
Deployment
**********

Docker
======

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

* Replace ``haniffalab/adifa:latest`` with your preferred image tag (see below)
* Replace ``/path/to/your/h5ad`` with the path to the directory where your ``.h5ad`` files are located
* Remove the ``entrypoint`` and ``command`` directives to prevent autodiscovery during the service boot
* Map the service to a different port on the host machine, for example ``- 5000:80`` to serve on port 80

Start the docker service
------------------------

Run the following command in the same directory:

:: 

    docker-compose up

This command will load the Adifa application, read your ``.h5ad`` files and start the web server. Depending on the size of your data files, the application may take several minutes to load. 

To access the application, open http://localhost:5000 in a browser. 

Offical Docker images
*********************

You can pull the official Adifa Docker images from the DockerHub Haniffa Lab 

The ``latest`` tag can be used for the most recent 

::

    docker pull haniffalab/adifa:latest

The ``dev`` tag can be used for the most recent 

::

    docker pull haniffalab/adifa:dev

Building your own custom docker images
**************************************

You can build a custom image from the root of the sc-afida repositry

::

    docker build -t davehorsfall/adifa-custom:latest .


docker-compose
**************

This is then the ``docker-compose.yml`` so we see it.

:: 

    version: '2'
    services:
    # Web server
    nginx: 
        image: nginx:latest
        container_name: nginx
        volumes:
        - ./nginx.conf:/etc/nginx/conf.d/default.conf
        ports:
        - 80:80
        - 443:443
    # Application
    adex:
        image: davehorsfall/cellatlas-api:merge
        container_name: adex
        volumes:
        - ./data:/data
        ports:
        - 5000:5000      
        environment:     
        - "DATA_PATH=/data"
        links:
        - nginx

Nginx
-----

This is then the ``conf.conf`` so we see it.

:: 

    version: '2'
    services:
    # Web server
    nginx: 
        image: nginx:latest
        container_name: nginx