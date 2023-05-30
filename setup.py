from setuptools import setup

setup(
    name="ChromProcess",
    vesion="0.0.1",
    author=["William E. Robinson"],
    packages=["ChromProcess"],
    install_requires=[
        "numpy >= 1.23.5",
        "scipy >= 1.9.3",
        "netcdf4 >= 1.6.0",
        "tomli >= 2.0.1",
        "maturin >= 1.0.1"
    ],
    extras_require={"dev": ["pdoc", "black"]},
)
