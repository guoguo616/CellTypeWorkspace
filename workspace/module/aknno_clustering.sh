#!/bin/bash
#!/bin/bash
#SBATCH --job-name=sc_preprocess
#SBATCH --cpus-per-task=5
#SBATCH --time=1:00:00
#SBATCH --mem=150G
#SBATCH --output=/home/platform/project/cell_type_workspace/cell_type_workspace_api/workspace/log/job_log/aknno_%j.out
#SBATCH --partition=highmem

module load GCC/11.2.0 OpenMPI/4.1.1 GCCcore/11.2.0
module load Python/3.9.6
module load R/4.2.0
echo 'modules loaded'

SCRIPT_DIR=$(dirname "$0")
export PYTHONPATH=/data2/platform/cell_type_workspace/venv/lib/python3.9/site-packages:$PYTHONPATH


# Check if the correct number of arguments are provided
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 inputfile outputdir species tissue_class"
    exit 1
fi

inputfile=$1
outputdir=$2
species=$3
tissue_class=$4

mkdir -p $outputdir

echo '######################## Unzipping input file ########################'
unzip $inputfile -d $outputdir/visium_data
if [ $? -ne 0 ]; then
    echo "Failed to unzip the input file"
    exit 1
fi

# Check if the extraction result is a single folder and no other files
extracted_items=$(find $outputdir/visium_data -mindepth 1 -maxdepth 1)
if [ $(echo "$extracted_items" | wc -w) -eq 1 ] && [ -d "$extracted_items" ]; then
    echo "Found a single folder in the extracted items"
    # Copy the contents out of the single folder to $outputdir/visium_data
    mv $extracted_items/* $outputdir/visium_data/
    # Remove the now empty folder
    rmdir $extracted_items
fi

echo '######################## Clustering ########################'
Rscript /home/platform/project/cell_type_workspace/cell_type_workspace_api/workspace/module/aknno.R "$outputdir/visium_data" "$outputdir"
if [ $? -ne 0 ]; then
    echo "Failed to run the aknno.R"
    exit 1
fi

echo 'clustering done'

source "$SCRIPT_DIR/cell_marker.sh" "$outputdir" "$species" "$tissue_class"
