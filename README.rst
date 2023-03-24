Tool for computing chi2 maps from .ROOT files
=============================================

.. image:: https://git.km3net.de/spenamartinez/chi2_tool/badges/master/pipeline.svg
    :target: https://git.km3net.de/spenamartinez/chi2_tool/pipelines

.. image:: https://git.km3net.de/spenamartinez/chi2_tool/badges/master/coverage.svg
    :target: https://spenamartinez.pages.km3net.de/chi2_tool/coverage

.. image:: https://git.km3net.de/examples/km3badges/-/raw/master/docs-latest-brightgreen.svg
    :target: https://spenamartinez.pages.km3net.de/chi2_tool


Installation
~~~~~~~~~~~~

It is recommended to first create an isolated virtualenvironment to not interfere
with other Python projects::

  git clone https://git.km3net.de/spenamartinez/chi2_tool
  cd chi2_tool
  python3 -m venv venv
  . venv/bin/activate

Install directly from the Git server via ``pip`` (no cloneing needed)::

  pip install git+https://git.km3net.de/spenamartinez/chi2_tool

Or clone the repository and run::

  make install

To install all the development dependencies, in case you want to contribute or
run the test suite::

  make install-dev
  make test


---

*Created with ``cookiecutter https://git.km3net.de/templates/python-project``*
