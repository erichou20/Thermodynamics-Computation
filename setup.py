from setuptools import setup, find_packages

with open("README.md", "r") as f:
    long_description = f.read()

setup(
    name="Thermodynamics-Computation",
    version="0.1.0",
    description= "A package for chemical engineers with thermodynamic computation functions",
    package_dir={"":"Thermodynamics-Computation"},
    package=find_packages(where='Thermodynamics-Computation'),
    long_description=long_description,
    long_description_content_type="text/markdown",
)