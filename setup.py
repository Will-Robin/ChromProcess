from setuptools import setup

setup(
    name="ChromProcess",
    vesion="0.1.0",
    author=["William E. Robinson"],
    packages=["ChromProcess"],
    install_requires=[
        "numpy >= 1.18",
        "scipy >= 1.6",
        "matplotlib >= 3.3",
        "netCDF4 >= 1.5",
    ],
)
