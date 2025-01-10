INDIR="/home/carlos/projects/singlecell_Dec2024/cellranger1/"
OUTDIR="/home/carlos/projects/singlecell_Dec2024/secondary_analysis/cellranger1_secondary1/"

PREPROCESS_SCRIPT_1="/home/carlos/projects/singlecell_Dec2024/analysis_code1/preprocess_1.py"
PREPROCESS_SCRIPT_2="/home/carlos/projects/singlecell_Dec2024/analysis_code1/remove_background_soupX.R"
PREPROCESS_SCRIPT_3="/home/carlos/projects/singlecell_Dec2024/analysis_code1/detect_doublets_scDblFinder.R"

samples=("A" "B" "C")

mamba activate singlecell1
#
#
#if [ ! -d $OUTDIR ]; then
#    mkdir -p $OUTDIR
#fi
#
#for sample in A B C; do
#    echo "Processing sample $sample"
#    python3 $PREPROCESS_SCRIPT_1 -s $sample -o $OUTDIR \
#        -f $INDIR/$sample/outs/filtered_feature_bc_matrix.h5 \
#        -r $INDIR/$sample/outs/raw_feature_bc_matrix.h5 -c 5 -m 8
#
#done



mamba deactivate

#for sample in A B C; do
#    echo "Doing SoupX for sample $sample"
#
#    Rscript $PREPROCESS_SCRIPT_2 \
#        --samplename $sample \
#        --rawdata $OUTDIR/$sample/$sample"_rawdata.npz" \
#        --preprocessed $OUTDIR/$sample/$sample"_filtdata.npz" \
#        --cells $OUTDIR/$sample/$sample"_cells.csv" \
#        --genes $OUTDIR/$sample/$sample"_genes.csv" \
#        --clusters $OUTDIR/$sample/$sample"_soupx_groups.csv" \
#        --outdir $OUTDIR/$sample
#done

for sample in A B C; do
    echo "Finding doublets for sample $sample"

    Rscript $PREPROCESS_SCRIPT_3 \
        --samplename $sample \
        --input $OUTDIR/$sample/$sample"_soupX.npz" \
        --outdir $OUTDIR/$sample \
        --seed 123
done