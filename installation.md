# Installation

## 1. Clone The Repository

```
git clone https://github.com/Will-Robin/ChromProcess.git
cd ChromProcess
```

## 2. Create and configure an environment

### Using Conda

Create a virtual environment:

```
conda create -n chromprocess python=3.10
```

Activate the virtual environment:

```
conda activate chromprocess
```

Install ChromProcess:

```
pip install -e .
```

### Using Pip and Venv

Create a virtual environment:

```
python -m venv chromprocess-env
```

Activate it:

#### Linux/macOS

```
source chromprocess-env/bin/activate
```

#### Windows

```
chromprocess-env\Scripts\activate
```

Install ChromProcess:

```
pip install -e .
```

### Using uv

Create a virtual environment:

```
uv sync
```

# Development Installation

To install ChromProcess together with development tools:

## pip

```
pip install -e ".[dev]"
```

## uv

```
uv sync --extra dev
```

This installs `mypy`, `ruff` and `pdoc`.


