<!--
[![tests](https://github.com/haniffalab/sci-adifa/actions/workflows/test-coverage.yml/badge.svg)](https://github.com/haniffalab/sci-adifa/actions/workflows/test-coverage.yml)
[![build](https://github.com/haniffalab/sci-adifa/actions/workflows/docker-release.yml/badge.svg)](https://github.com/haniffalab/sci-adifa/actions/workflows/docker-release.yml)
[![docker](https://github.com/haniffalab/sci-adifa/actions/workflows/docker-latest.yml/badge.svg)](https://github.com/haniffalab/sci-adifa/actions/workflows/docker-latest.yml)
[![sphinx](https://github.com/haniffalab/sci-adifa/actions/workflows/sphinx-build.yml/badge.svg)](https://github.com/haniffalab/sci-adifa/actions/workflows/sphinx-build.yml)
[![codecov](https://codecov.io/gh/haniffalab/sci-adifa/branch/main/graph/badge.svg?token=RQLL0HKQ5W)](https://codecov.io/gh/haniffalab/sci-adifa) -->

# sci-adifa-elmer [![python](https://img.shields.io/badge/python-3.8-blue)](https://python.org)

Authors: Antony Rose and Daniela Basurto-Lozada
***

sci-adifa-elmer is a forked repository from the original Adifa software which this adaptation is built upon. Link to the original software and information below.

This developmental branch builds upon the original Adifa software to add a new window for pseudo-spatial visualisation for droplet based scRNA-seq data that contains locality metadata.

This includes three new visualisations for users to explore their droplet based scRNA-seq data.

Guides and notebooks for mask generation and adapter notebook to process scRNA-seq ready for input into the sci-adifa-elmer will be available at [(https://github.com/arose20/Pseudo-spatial_visualisation_of_scRNAseq_data)](https://github.com/arose20/Pseudo-spatial_visualisation_of_scRNAseq_data)

<br> <br>


## SINGLE CELL INSIGHTS

## Adifa software information

[![docs](https://img.shields.io/badge/Documentation-online-blue)](https://haniffalab.github.io/sci-adifa)
[![doi](https://zenodo.org/badge/DOI/10.5281/zenodo.5824895.svg)](https://doi.org/10.5281/zenodo.5824895)

Adifa is a framework for visualising single-cell gene expression data in a web browser. It is built on [Flask](https://flask.palletsprojects.com/), a micro web framework, and ingests [Annotated Data](https://anndata.readthedocs.io/) objects in the `.h5ad` file format. It includes dimensionality reduction visualisation, heatmaps and dotplots, with the ability to explore gene expression and disease markers. The Python-based implementation and usage of the [deck.gl](https://deck.gl/) framework allows efficient handling of datasets up to one million cells.

Link to original Adifa repository: [https://github.com/haniffalab/sci-adifa](https://github.com/haniffalab/sci-adifa)

Discuss usage on [Discourse]. Read the [documentation](https://haniffalab.github.io/sci-adifa). If you’d like to contribute by opening an issue or creating a pull request, please take a look at our [contributing guide](CONTRIBUTING.md). If Adifa is useful for your research, consider [citing the software](https://haniffalab.com/sci-adifa/citing.html). 
