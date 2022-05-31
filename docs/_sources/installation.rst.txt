Installation
============

Follow these instructions to install a local installation of the Adifa source code.

.. attention::

    This page is for developers. If you want to use Adifa to visualise your data, please visit the :ref:`quickstart <quickstart>` page and follow the instructions there. 

Prerequisites
-------------

Before you begin using Adifa, make sure you have installed the following libraries:

- Python (>=3.8)
- Git

Install
-------------

**Be sure to use the same version of the code as the version of the docs
you're reading.** You probably want the latest tagged version, but the
default Git version is the main branch. 

Clone the repository

::

    $ git clone git@github.com:haniffalab/sci-adifa.git
    $ cd sci-adifa

Create a virtualenv and activate it

::

    $ python -m venv venv
    $ . venv/bin/activate

Or on Windows cmd

::

    $ python -m venv venv
    $ venv\Scripts\activate.bat

Install the requirements

::

    $ pip install -r requirements.txt

Run
---

::

    $ flask init-db
    $ flask autodiscover
    $ flask run

Or on Windows cmd::

    $ flask init-db
    $ flask autodiscover
    $ flask run

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

