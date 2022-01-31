from setuptools import setup, find_packages

with open("README.md", "rb") as fh:
    long_description = fh.read().decode()

with open("requirements.txt") as fh:
    requirements = fh.read().splitlines()

with open("requirements-extra.txt") as fh:
    requirements_extra = fh.read().splitlines()

setup(
    name="adifa",
    version="0.0.2",
    packages=find_packages(),
    url="https://haniffalab.github.io/sci-adifa",
    license="MIT",
    author="Haniffa Lab",
    author_email="rse@haniffalab.com",
    description="Adifa is a framework for visualising single-cell gene expression data in a web browser",
    long_description=long_description,
    long_description_content_type="text/markdown",
    install_requires=requirements,
    python_requires=">=3.8",
    include_package_data=True,
    zip_safe=False,
    classifiers=[
        "Development Status :: 1 - Planning",
        "Framework :: Flask",
        "Environment :: Web Environment",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3.8",
        "Natural Language :: English",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    extras_require=dict(prepare=requirements_extra),
)