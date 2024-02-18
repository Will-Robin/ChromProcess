from setuptools import setup

setup(
    name="ChromProcess",
    vesion="0.0.2",
    author=["William E. Robinson"],
    packages=["ChromProcess"],
    install_requires=[
        "numpy >= 1.26.4",
        "scipy >= 1.12.0",
        "netcdf4 >= 1.6.5",
        "tomli >= 2.0.1",
    ],
    extras_require={"dev": ["pdoc", "black", "ruff"]},
)
