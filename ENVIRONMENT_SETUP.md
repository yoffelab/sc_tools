# Environment Setup Guide

This guide explains how to set up and share your Python environment for this project, similar to working with conda environments.

## Using Conda (Recommended for Scientific Computing)

### Creating the environment

```bash
# Create environment from environment.yml
conda env create -f environment.yml

# Activate the environment
conda activate ggo_visium
```

### Sharing with others

To share your exact environment with others:

1. **Export your current environment:**
   ```bash
   conda env export > environment.yml
   ```

2. **Or export without build strings (more portable):**
   ```bash
   conda env export --no-builds > environment.yml
   ```

3. **Share the `environment.yml` file** - others can recreate it with:
   ```bash
   conda env create -f environment.yml
   ```

## Using pip/venv

### Creating a virtual environment

```bash
# Create virtual environment
python -m venv venv

# Activate it (macOS/Linux)
source venv/bin/activate

# Or on Windows
# venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt
```

### Sharing with others

1. **Generate requirements.txt from your current environment:**
   ```bash
   pip freeze > requirements.txt
   ```

2. **Share the `requirements.txt` file** - others can install with:
   ```bash
   pip install -r requirements.txt
   ```

## Using the Environment in Cursor

### Option 1: Select Python Interpreter

1. In Cursor, press `Cmd+Shift+P` (macOS) or `Ctrl+Shift+P` (Windows/Linux)
2. Type "Python: Select Interpreter"
3. Choose your conda environment (e.g., `~/anaconda3/envs/ggo_visium/bin/python`)
   - Or your venv: `./venv/bin/python`

### Option 2: Terminal in Cursor

Cursor uses your system terminal, so you can:

1. Open the integrated terminal (`Ctrl+`` ` or `View > Terminal`)
2. Activate your environment:
   ```bash
   conda activate ggo_visium
   # or
   source venv/bin/activate
   ```
3. Run your scripts normally:
   ```bash
   python scripts/basic_analysis.py
   ```

## Running Scripts

Once your environment is activated (either in terminal or selected in Cursor):

```bash
# From project root
python scripts/basic_analysis.py
python scripts/preprocessing.py
# etc.
```

## Tips

- **Lock versions for reproducibility:** Use `conda env export --from-history` to get only explicitly installed packages
- **Document custom packages:** If you use packages not in PyPI/conda (like `sthd`), document installation instructions
- **Python version:** Specify Python version in `environment.yml` for consistency
- **Platform-specific packages:** Use `environment.yml` for conda (handles binaries better) and `requirements.txt` for pip-only setups

