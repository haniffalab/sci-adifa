
[![tests](https://github.com/haniffalab/sci-adifa/actions/workflows/test-coverage.yml/badge.svg)](https://github.com/haniffalab/sci-adifa/actions/workflows/test-coverage.yml)
[![build](https://github.com/haniffalab/sci-adifa/actions/workflows/docker-release.yml/badge.svg)](https://github.com/haniffalab/sci-adifa/actions/workflows/docker-release.yml)
[![docker](https://github.com/haniffalab/sci-adifa/actions/workflows/docker-latest.yml/badge.svg)](https://github.com/haniffalab/sci-adifa/actions/workflows/docker-latest.yml)
[![sphinx](https://github.com/haniffalab/sci-adifa/actions/workflows/sphinx-build.yml/badge.svg)](https://github.com/haniffalab/sci-adifa/actions/workflows/sphinx-build.yml)
[![codecov](https://codecov.io/gh/haniffalab/sci-adifa/branch/main/graph/badge.svg?token=RQLL0HKQ5W)](https://codecov.io/gh/haniffalab/sci-adifa)
[![python](https://img.shields.io/badge/python-3.8-blue)](https://python.org)

## SINGLE CELL INSIGHTS

# Adifa - Annotated Data in Flask App

[![docs](https://img.shields.io/badge/Documentation-online-blue)](https://haniffalab.github.io/sci-adifa)
[![doi](https://zenodo.org/badge/DOI/10.5281/zenodo.5824895.svg)](https://doi.org/10.5281/zenodo.5824895)

Adifa is a framework for visualising single-cell gene expression data in a web browser. It is built on [Flask](https://flask.palletsprojects.com/), a micro web framework, and ingests [Annotated Data](https://anndata.readthedocs.io/) objects in the `.h5ad` file format. It includes dimensionality reduction visualisation, heatmap and dotplots, and the ability to explore gene expression and disease markers. The Python-based implementation and usage of [deck.gl](https://deck.gl/) deals efficiently with datasets of up to one million cells.

If you'd like to publish your final publication data objects online please [get in touch]. Discuss usage on [Discourse]. Read the [documentation](https://haniffalab.github.io/sci-adifa). If you’d like to contribute by opening an issue or creating a pull request, please take a look at our [contributing guide](CONTRIBUTING.md). If Adifa is useful for your research, consider [citing the software](https://haniffalab.com/sci-adifa/citing.html). 