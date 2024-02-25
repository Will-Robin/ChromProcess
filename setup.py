from setuptools import setup
import ChromProcess

setup(
    name="ChromProcess",
    version=ChromProcess.__version__,
    author=ChromProcess.__author__,
    packages=["ChromProcess"],
    install_requires=[
        "numpy >= 1.26.4",
        "scipy >= 1.12.0",
        "netcdf4 >= 1.6.5",
        "tomli >= 2.0.1",
    ],
    extras_require={"dev": ["pdoc", "black", "ruff"]},
)
