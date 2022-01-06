Installation
============

The basic blog app built in the Flask `tutorial`_.

.. _tutorial: https://flask.palletsprojects.com/tutorial/


Prerequisites
-------------
Python 3.8

Install
-------

**Be sure to use the same version of the code as the version of the docs
you're reading.** You probably want the latest tagged version, but the
default Git version is the main branch. ::

    # clone the repository
    $ git clone git@github.com:haniffalab/sci-adifa.git
    $ cd sci-adifa

Create a virtualenv and activate it::

    $ python -m venv venv
    $ . venv/bin/activate

Or on Windows cmd::

    $ python -m venv venv
    $ venv\Scripts\activate.bat

Install adifa::

    $ pip install -e .
    $ pip install -r requirements.txt


Run
---

::

    $ export FLASK_APP=adifa
    $ export FLASK_ENV=development
    $ flask init-db
    $ flask run

Or on Windows cmd::

    > set FLASK_APP=adifa
    > set FLASK_ENV=development
    > flask init-db
    > flask run

Open http://127.0.0.1:5000 in a browser.


Test
----

::

    $ pip install '.[test]'
    $ pip install pytest coverage
    $ pytest

Run with coverage report::

    $ coverage run -m pytest
    $ coverage report
    $ coverage html  # open htmlcov/index.html in a browser

Sphinx documentation
--------------------

::

    $ cd sphinx
    $ make clean
    $ make html # open _build/html/index.html in a browser

Docker builds
--------------------

::

    $ cd sphinx
    $ make clean
    $ make html # open _build/html/index.html in a browser    
