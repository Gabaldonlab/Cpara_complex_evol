# Get sample names
source "scripts/utils/sample_names.sh"
source "scripts/utils/set_paths.sh"

module load bedtools

# CNV path
CNV_PATH="results/cnv_calling"
OUT_PATH="results/cnv_intersect"
ANNOT_PATH="results/annot_cnv"
DAT_PATH="data/aneuploidies"

# Parse CNV outputs to BED format
find $CNV_PATH -type f -name final_CNVcalling.tab -execdir sh -c \
    'echo {}; awk -v OFS="\t" "{print \$1, \$3, \$4, \$2}" {} | grep -v "^chromosome" > {}.bed' \;

# Filter T0 CNVs
# Loop start at 2 because 0 is the T0 and 1 the YPD
for (( i=2; i<${#BP57_STRAINS[@]}; i++ )); do
    strain="${BP57_STRAINS[$i]}"
    t0strain="${BP57_STRAINS[0]}"
    ypdstrain="${BP57_STRAINS[1]}"
    echo -e $strain "-\t" $t0strain "\t" $ypdstrain
    bedtools subtract \
        -a $CNV_PATH/${strain}/final_CNVcalling.tab.bed \
        -b $CNV_PATH/${t0strain}/final_CNVcalling.tab.bed > \
        $OUT_PATH/${strain}/final_CNVcalling_woT0.bed
    bedtools subtract \
        -a $OUT_PATH/${strain}/final_CNVcalling_woT0.bed \
        -b $CNV_PATH/${ypdstrain}/final_CNVcalling.tab.bed > \
        $OUT_PATH/${strain}/final_CNVcalling_woYPD.bed
    ln -frs $OUT_PATH/${strain}/final_CNVcalling_woYPD.bed \
        $OUT_PATH/${strain}/final_CNVcalling.wobg.bed
done

for (( i=2; i<${#DMKU_13A_STRAINS[@]}; i++ )); do
    strain="${DMKU_13A_STRAINS[$i]}"
    t0strain="${DMKU_13A_STRAINS[0]}"
    ypdstrain="${DMKU_13A_STRAINS[1]}"
    echo -e $strain "-\t" $t0strain "\t" $ypdstrain
    bedtools subtract \
        -a $CNV_PATH/${strain}/final_CNVcalling.tab.bed \
        -b $CNV_PATH/${t0strain}/final_CNVcalling.tab.bed > \
        $OUT_PATH/${strain}/final_CNVcalling_woT0.bed
    bedtools subtract \
        -a $OUT_PATH/${strain}/final_CNVcalling_woT0.bed \
        -b $CNV_PATH/${ypdstrain}/final_CNVcalling.tab.bed > \
        $OUT_PATH/${strain}/final_CNVcalling_woYPD.bed
    ln -frs $OUT_PATH/${strain}/final_CNVcalling_woYPD.bed \
        $OUT_PATH/${strain}/final_CNVcalling.wobg.bed
done


for (( i=2; i<${#CANOR_C1_STRAINS[@]}; i++ )); do
    strain="${CANOR_C1_STRAINS[$i]}"
    t0strain="${CANOR_C1_STRAINS[0]}"
    ypdstrain="${CANOR_C1_STRAINS[1]}"
    echo -e $strain "-\t" $t0strain "\t" $ypdstrain
    bedtools subtract \
        -a $CNV_PATH/${strain}/final_CNVcalling.tab.bed \
        -b $CNV_PATH/${t0strain}/final_CNVcalling.tab.bed > \
        $OUT_PATH/${strain}/final_CNVcalling_woT0.bed
    bedtools subtract \
        -a $OUT_PATH/${strain}/final_CNVcalling_woT0.bed \
        -b $CNV_PATH/${ypdstrain}/final_CNVcalling.tab.bed > \
        $OUT_PATH/${strain}/final_CNVcalling_woYPD.bed
    ln -frs $OUT_PATH/${strain}/final_CNVcalling_woYPD.bed \
        $OUT_PATH/${strain}/final_CNVcalling.wobg.bed
done

for (( i=2; i<${#CANOR_C2_STRAINS[@]}; i++ )); do
    strain="${CANOR_C2_STRAINS[$i]}"
    t0strain="${CANOR_C2_STRAINS[0]}"
    ypdstrain="${CANOR_C2_STRAINS[1]}"
    echo -e $strain "-\t" $t0strain "\t" $ypdstrain
    bedtools subtract \
        -a $CNV_PATH/${strain}/final_CNVcalling.tab.bed \
        -b $CNV_PATH/${t0strain}/final_CNVcalling.tab.bed > \
        $OUT_PATH/${strain}/final_CNVcalling_woT0.bed
    bedtools subtract \
        -a $OUT_PATH/${strain}/final_CNVcalling_woT0.bed \
        -b $CNV_PATH/${ypdstrain}/final_CNVcalling.tab.bed > \
        $OUT_PATH/${strain}/final_CNVcalling_woYPD.bed
    ln -frs $OUT_PATH/${strain}/final_CNVcalling_woYPD.bed \
        $OUT_PATH/${strain}/final_CNVcalling.wobg.bed
done

for (( i=2; i<${#CANOR_C3_STRAINS[@]}; i++ )); do
    strain="${CANOR_C3_STRAINS[$i]}"
    t0strain="${CANOR_C3_STRAINS[0]}"
    ypdstrain="${CANOR_C3_STRAINS[1]}"
    echo -e $strain "-\t" $t0strain "\t" $ypdstrain
    bedtools subtract \
        -a $CNV_PATH/${strain}/final_CNVcalling.tab.bed \
        -b $CNV_PATH/${t0strain}/final_CNVcalling.tab.bed > \
        $OUT_PATH/${strain}/final_CNVcalling_woT0.bed
    bedtools subtract \
        -a $OUT_PATH/${strain}/final_CNVcalling_woT0.bed \
        -b $CNV_PATH/${ypdstrain}/final_CNVcalling.tab.bed > \
        $OUT_PATH/${strain}/final_CNVcalling_woYPD.bed
    ln -frs $OUT_PATH/${strain}/final_CNVcalling_woYPD.bed \
        $OUT_PATH/${strain}/final_CNVcalling.wobg.bed
done

for (( i=2; i<${#CANOR_C4_STRAINS[@]}; i++ )); do
    strain="${CANOR_C4_STRAINS[$i]}"
    t0strain="${CANOR_C4_STRAINS[0]}"
    ypdstrain="${CANOR_C4_STRAINS[1]}"
    echo -e $strain "-\t" $t0strain "\t" $ypdstrain
    bedtools subtract \
        -a $CNV_PATH/${strain}/final_CNVcalling.tab.bed \
        -b $CNV_PATH/${t0strain}/final_CNVcalling.tab.bed > \
        $OUT_PATH/${strain}/final_CNVcalling_woT0.bed
    bedtools subtract \
        -a $OUT_PATH/${strain}/final_CNVcalling_woT0.bed \
        -b $CNV_PATH/${ypdstrain}/final_CNVcalling.tab.bed > \
        $OUT_PATH/${strain}/final_CNVcalling_woYPD.bed
    ln -frs $OUT_PATH/${strain}/final_CNVcalling_woYPD.bed \
        $OUT_PATH/${strain}/final_CNVcalling.wobg.bed
done

for (( i=2; i<${#CANOR_PA_STRAINS[@]}; i++ )); do
    strain="${CANOR_PA_STRAINS[$i]}"
    t0strain="${CANOR_PA_STRAINS[0]}"
    ypdstrain="${CANOR_PA_STRAINS[1]}"
    echo -e $strain "-\t" $t0strain "\t" $ypdstrain
    bedtools subtract \
        -a $CNV_PATH/${strain}/final_CNVcalling.tab.bed \
        -b $CNV_PATH/${t0strain}/final_CNVcalling.tab.bed > \
        $OUT_PATH/${strain}/final_CNVcalling_woT0.bed
    bedtools subtract \
        -a $OUT_PATH/${strain}/final_CNVcalling_woT0.bed \
        -b $CNV_PATH/${ypdstrain}/final_CNVcalling.tab.bed > \
        $OUT_PATH/${strain}/final_CNVcalling_woYPD.bed
    ln -frs $OUT_PATH/${strain}/final_CNVcalling_woYPD.bed \
        $OUT_PATH/${strain}/final_CNVcalling.wobg.bed
done

for (( i=2; i<${#CANOR_PB_STRAINS[@]}; i++ )); do
    strain="${CANOR_PB_STRAINS[$i]}"
    t0strain="${CANOR_PB_STRAINS[0]}"
    ypdstrain="${CANOR_PB_STRAINS[1]}"
    echo -e $strain "-\t" $t0strain "\t" $ypdstrain
    bedtools subtract \
        -a $CNV_PATH/${strain}/final_CNVcalling.tab.bed \
        -b $CNV_PATH/${t0strain}/final_CNVcalling.tab.bed > \
        $OUT_PATH/${strain}/final_CNVcalling_woT0.bed
    bedtools subtract \
        -a $OUT_PATH/${strain}/final_CNVcalling_woT0.bed \
        -b $CNV_PATH/${ypdstrain}/final_CNVcalling.tab.bed > \
        $OUT_PATH/${strain}/final_CNVcalling_woYPD.bed
    ln -frs $OUT_PATH/${strain}/final_CNVcalling_woYPD.bed \
        $OUT_PATH/${strain}/final_CNVcalling.wobg.bed
done

# Subtract switch treatments 
for (( i=0; i<${#CANME_SWITCH_1[@]}; i++ )); do
    switched_strain="${CANME_SWITCH_2[$i]}"
    first_strain="${CANME_SWITCH_1[$i]}"
    echo -e $switched_strain "-\t" $first_strain
    bedtools subtract \
        -a $OUT_PATH/${switched_strain}/final_CNVcalling_woYPD.bed \
        -b $CNV_PATH/${first_strain}/final_CNVcalling.tab.bed > \
        $OUT_PATH/${switched_strain}/final_CNVcalling_woSWITCH.bed
    ln -frs $OUT_PATH/${switched_strain}/final_CNVcalling_woSWITCH.bed \
        $OUT_PATH/${switched_strain}/final_CNVcalling.wobg.bed
done

# Subtract switch treatments 
for (( i=0; i<${#CANOR_SWITCH_1[@]}; i++ )); do
    switched_strain="${CANOR_SWITCH_2[$i]}"
    first_strain="${CANOR_SWITCH_1[$i]}"
    echo -e $switched_strain "-\t" $first_strain
    bedtools subtract \
        -a $OUT_PATH/${switched_strain}/final_CNVcalling_woYPD.bed \
        -b $CNV_PATH/${first_strain}/final_CNVcalling.tab.bed > \
        $OUT_PATH/${switched_strain}/final_CNVcalling_woSWITCH.bed
    ln -frs $OUT_PATH/${switched_strain}/final_CNVcalling_woSWITCH.bed \
        $OUT_PATH/${switched_strain}/final_CNVcalling.wobg.bed
done

# Subtract anueploidies from the background-removed tables 
for (( i=0; i<${#SAMPLE_NAMES[@]}; i++ )); do
    strain="${SAMPLE_NAMES[$i]}"
    echo $strain
    if [ -s $DAT_PATH/${strain}_aneuploidies.bed ]
    then
        echo "removing aneups"
        bedtools subtract \
            -a $OUT_PATH/${strain}/final_CNVcalling.wobg.bed \
            -b $DAT_PATH/${strain}_aneuploidies.bed > \
            $OUT_PATH/${strain}/final_CNVcalling_woAneup.bed
        ln -frs $OUT_PATH/${strain}/final_CNVcalling_woAneup.bed \
            $OUT_PATH/${strain}/final_CNVcalling.final.bed
    else
        echo "not removing aneups"
        ln -frs $OUT_PATH/${strain}/final_CNVcalling.wobg.bed \
            $OUT_PATH/${strain}/final_CNVcalling.final.bed
    fi
done