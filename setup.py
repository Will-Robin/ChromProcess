from setuptools import setup

setup(
    name="ChromProcess",
    vesion="0.0.0",
    author=["William E. Robinson"],
    packages=["ChromProcess"],
    install_requires=[
        "numpy >= 1.21.1",
        "scipy >= 1.7.3",
        "netcdf4 >= 1.5.7",
        "tomli >= 2.0.1",
    ],
)
