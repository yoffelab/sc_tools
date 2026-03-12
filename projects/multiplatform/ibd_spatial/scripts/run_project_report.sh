#!/bin/bash
#SBATCH --job-name=ibd_report
#SBATCH --partition=scu-cpu
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --time=02:00:00
#SBATCH --output=/home/fs01/juk4007/elementolab/projects/ibd_spatial/logs/project_report.log

set -eo pipefail
source /etc/profile || true
eval "$(conda shell.bash hook 2>/dev/null)"
conda activate sc_tools

WORKDIR=/home/fs01/juk4007/elementolab/projects/ibd_spatial
SC_TOOLS=/home/fs01/juk4007/elementolab/projects/ibd_spatial/sc_tools_repo

# Install sc_tools if not already available
python -c "import sc_tools" 2>/dev/null || {
    echo "Installing sc_tools..."
    pip install -e "${SC_TOOLS}[benchmark]" --no-deps 2>/dev/null || \
    pip install -e "${SC_TOOLS}" --no-deps 2>/dev/null || true
}

pip install jinja2 plotly 2>/dev/null || true

mkdir -p ${WORKDIR}/figures/QC

python ${WORKDIR}/scripts/generate_project_report.py --milestone all
