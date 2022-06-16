Configuration
=============

Applications need some kind of configuration to toggle settings and define other such environment-specific things. Adifa handles the config like other Flask applications, with the `Config`_ class.

.. _Config: https://flask.palletsprojects.com/en/2.0.x/config/


Default Configuration
---------------------

The default configuration for Adifa is defined in the ``.env`` file in the `repository root`_, and contains the following values.

::

    SECRET_KEY = 'best-password-here'
    SQLALCHEMY_DATABASE_URI = 'sqlite:///../instance/adifa.sqlite'
    SQLALCHEMY_TRACK_MODIFICATIONS = False
    API_VERSION = 1
    API_PREFIX = '/api'
    API_SERVER = '/'
    DATA_PATH = './instance'

.. _repository root: https://github.com/haniffalab/sci-adifa/blob/main/.env

Overriding the Configuration
----------------------------

To override the deafult configuration you can either:

* use system environment variables (recommended for production)
* creating a file called ``config.py`` in your ``instance`` folder (recommended for development)

Using environment variables
^^^^^^^^^^^^^^^^^^^^^^^^^^^

To override indivisual values in the configuration pass or create a new system environment variable with the same name as the configuration value, 
and assign it the value that is needed. 

Using config.py
^^^^^^^^^^^^^^^

The instance folder was introduced in Flask 0.8, and is designed to not be under version control and be deployment specific. It’s the perfect place
to drop things that either change at runtime or configuration files. To override the configuration, create a new file called ``config.py`` in your
``instance`` folder. This file will not be placed under version control. Copy the contents from ``.env`` into your new ``config.py``, and update the 
values as needed.

Loading new datasets
--------------------

Data folder
^^^^^^^^^^^

Adifa reads files from the data folder to populate a database which it will query to display their information and generate plots.
This data folder is specified through the ``DATA_PATH`` environment variable or through the ``config.py`` file.
You should store your ``.h5ad`` files in this folder.

Autodiscovery
^^^^^^^^^^^^^

The autodiscovery tool manages the whole process to add datasets:
it searches for ``.h5ad`` files within the data folder and reads them to populate the database.
To run this tool use

::

    $ flask autodiscover

.. highlight:: shell


For each file in the data folder you can expect an output like


::

    [2022-05-27 21:57:21,711] INFO in dataset_utils: Inspecting file.h5ad
    [2022-05-27 21:57:22,026] INFO in dataset_utils: Hashing file.h5ad
    [2022-05-27 21:57:23,149] INFO in dataset_utils: Adding file.h5ad

.. highlight:: shell


Whenever you add new files to the data folder you must run the autodiscover tool again for them to be added.
Files that had already been added will have a slightly different output like


::

    [2022-05-27 21:57:21,711] INFO in dataset_utils: Inspecting file.h5ad
    [2022-05-27 21:57:22,026] INFO in dataset_utils: Hashing file.h5ad
    [2022-05-27 21:57:23,149] INFO in dataset_utils: Reprocessed file.h5ad

.. highlight:: shell
