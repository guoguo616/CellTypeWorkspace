#!/bin/bash
#!/bin/bash
#SBATCH --job-name=sc_preprocess
#SBATCH --cpus-per-task=5
#SBATCH --time=1:00:00
#SBATCH --mem=150G
#SBATCH --output=/home/platform/project/cell_type_workspace/cell_type_workspace_api/workspace/log/job_log/aknno_%j.out
#SBATCH --partition=highmem

# If the script is run directly, load the required modules
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    module load GCC/11.2.0 OpenMPI/4.1.1 GCCcore/11.2.0
    module load Python/3.9.6
    module load R/4.2.0
    echo 'modules loaded'
fi

SCRIPT_DIR=$(dirname "$0")
export PYTHONPATH=/data2/platform/cell_type_workspace/venv/lib/python3.9/site-packages:$PYTHONPATH


# Check if the correct number of arguments are provided
if [ "$#" -lt 3 ] || [ "$#" -gt 4 ]; then
    echo "Usage: $0 outputdir species tissue_class [override_cluster_json]"
    exit 1
fi

outputdir=$1
species=$2
tissue_class=$3
override_cluster_json=$4

echo '######################## Finding Cell Marker ########################'
echo 'running R script to find marker genes and python script to filter adata'
Rscript /home/platform/project/cell_type_workspace/cell_type_workspace_api/workspace/module/cell_marker.R \
    "$outputdir" "$species" "$tissue_class" "$override_cluster_json"
if [ $? -ne 0 ]; then
    echo "Failed to run the aknno.R"
    exit 1
fi

# Call python script to generate gexf
echo '######################## Generating gexf ########################'
/data2/platform/cell_type_workspace/venv/bin/python3 \
    /home/platform/project/cell_type_workspace/cell_type_workspace_api/workspace/module/write_gexf.py \
    --visium_path "$outputdir/visium_data" \
    --output_path "$outputdir" \
    --is_overridden "$(if [ -n "$override_cluster_json" ]; then echo 'True'; else echo 'False'; fi)"
if [ $? -ne 0 ]; then
    echo "Failed to run the write_gexf.py"
    exit 1
fi

echo 'cell marker done'
