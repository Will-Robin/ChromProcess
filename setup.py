from setuptools import setup

setup(
    name="ChromProcess",
    vesion="0.0.1",
    author=["William E. Robinson"],
    packages=["ChromProcess"],
    install_requires=[
        "numpy >= 1.24.3",
        "scipy >= 1.10.1",
        "netcdf4 >= 1.6.2",
        "tomli >= 2.0.1",
    ],
    extras_require={"dev": ["pdoc", "black", "ruff"]},
)
