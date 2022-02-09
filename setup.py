from glob import glob
from distutils.core import setup

reqs = [
    "astropy",
    "numpy",
    "scipy",
    "matplotlib",
    "mwa_pb_lookup",
    "calplots",
    "casacore",
    "healpy",
    "requests",
    "mysqlconnector",
    "psutil",
]

scripts = glob('gleam_x/*/*py')

setup(
    name="gleam_x",
    version="0.1",
    author="Natasha Hurley-Walker, Paul Hancock, Tim Galvin",
    description="Python scripts to support the processing of GLEAM-X data.",
    url="https://github.com/nhurleywalker/GLEAM-X-pipeline",
    long_description=open("README.md").read(),
    packages=["gleam_x", "gleam_x.bin", "gleam_x.db", "gleam_x.utils"],
    requires=reqs,
    scripts=scripts,
)
