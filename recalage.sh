#!/bin/bash
# VESTIBULOCOCHLEAR TRACTOGRAPHY - Antoine Lachance (2025-2026)

# This script performs ensemble tractography of the CN VIII using fODFs and FA maps
# from TractoFlow. Anatomical ROIs from both FreeSurfer and the ROIs_mean MNI folder are used to
# segment and virtually dissect the vestibulocochlear nerve. Tracking is performed in original space over a
# grid of (step size, theta) combinations, merged, transformed to MNI space, filtered and segmented 
# into three components (cochlear bundle, inferior vestibular bundle, superior vestibular bundle). Final cleaned bundles are saved
# in both original (orig_space) and MNI (mni_space) coordinate spaces.
#FA = 0.2 is used for HCP data; for clinical data, a lower FA threshold (e.g., 0.15) should be considered due to reduced anisotropy.

# EXAMPLE COMMAND
#
#   bash vestibulocochlear_first_order.sh \
#       -s /path/to/subject/folder/sub-01 \
#       -m /path/to/ROIs_mean \
#       -o /path/to/output_dir \
#       -t 10 \
#       -p 0.5 \
#       -e 30 \
#       -f 0.15 \     
#       -n 20000 \
#       -g true
#
# Notes:
#   - If -p (step) and -e (theta) are not provided, the script performs ensemble tracking using:
#       step = {0.1, 0.5, 1.0}
#       theta = {20, 30, 40}
#   - GPU option (-g) accelerates local tracking when supported; omit it to run on CPU.



# Input structure
#
#    [input]
#    ├── sub-01
#    │   ├── freesurfer
#    │   │   └─── aparc.DKTatlas+aseg.mgz
#    │   ├── sub-01__fa.nii.gz
#    │   ├── sub-01__fodf.nii.gz
#    │   └── sub-01__T1w.nii.gz
#    │
#    ├── S2
#    .
#    .


# -m : Path to the MNI-space reference folder.
#       In this project, it should point to the path/to/vestibulocochlear/ROIs_clean/ folder
#       within the vestibulocochlear GitHub project.

usage() { 
    echo "$(basename $0) [-s path/to/subject] [-m path/to/mni] [-o output_dir] [-t nb_threads] [-p step_size] [-e theta_deg] [-f fa_threshold] [-n npv] -g true" 1>&2
    exit 1
}


while getopts "s:m:o:t:p:e:f:n:g:" args; do
    case "${args}" in
        s) subject_dir=${OPTARG} ;;
        m) mni_dir=${OPTARG} ;;
        o) output_dir=${OPTARG} ;;
        t) nb_threads=${OPTARG} ;;
        p) step_size=${OPTARG} ;;
        e) theta=${OPTARG} ;;
        f) fa_threshold=${OPTARG} ;;   
        n) npv=${OPTARG} ;;    
        g) gpu=${OPTARG} ;;
        *) usage;;
    esac
done
shift $((OPTIND-1))

if [ -z "${subject_dir}" ] || [ -z "${mni_dir}" ] || [ -z "${output_dir}" ]; then
    usage
fi




# If the user provides both step size (-p) and theta (-e) values, the script will use those specific parameters and perform a single run
# of local_tracking (i.e., no ensemble tractography will be performed).

if [ -n "${step_size}" ] && [ -n "${theta}" ]; then
    step_list=(${step_size})
    theta_list=(${theta})
else
    step_list=(0.1 0.5 1.0)
    theta_list=(20 30 40)
fi

# Default values if not set
fa_threshold=${fa_threshold:-0.15}
npv=${npv:-20000}


# npv is the total number of seeds per voxel.
# It is divided by the number of step/theta combinations.
# The resulting npv_per_run is the number of seeds actually used in each run.



# number of step/theta combos
n_combos=$(( ${#step_list[@]} * ${#theta_list[@]} ))
# seeds per combo (rounded)
npv_per_run=$(( (npv + n_combos - 1) / n_combos ))  # ceiling division
echo "Using $npv_per_run seeds per voxel per run (based on $npv total)"



#opposite_side=leftright

# Enable GPU if available (highly recommended for X times faster processing)

if [ ! -z "${gpu}" ]; then
    gpu="--use_gpu"
else
    gpu=""
fi

echo "Folder subject: " ${subject_dir}
echo "Folder MNI: " ${mni_dir}
echo "Output folder: " ${output_dir}
echo "Use GPU: " ${gpu}
echo "Number of threads" ${nb_threads}
echo "Tracking grid: steps=${step_list[*]}  thetas=${theta_list[*]}"


nsub=$(basename "${subject_dir}")
subject_parent=$(dirname "${subject_dir}")

    





mkdir -p ${output_dir}/${nsub}/orig_space/{rois,tracking_vestibulocochlear,transfo}
mkdir -p ${output_dir}/${nsub}/mni_space/{rois,tracking_vestibulocochlear}

# Per-combo dirs will be created later, inside the loops.
orig_rois_dir=${output_dir}/${nsub}/orig_space/rois
mni_rois_dir=${output_dir}/${nsub}/mni_space/rois
orig_tracking_dir=${output_dir}/${nsub}/orig_space/tracking_vestibulocochlear
mni_tracking_root=${output_dir}/${nsub}/mni_space/tracking_vestibulocochlear
trials_dir=${orig_tracking_dir}/trials
mkdir -p "${trials_dir}"

echo ""
echo "|------------- PROCESSING CNVIII TRACTOGRAPHY FOR ${nsub} -------------|"
echo ""

echo "|------------- 1) Registrations -------------|"
export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS="${nb_threads}"

echo "Running ANTs registration... (logs saved to ${output_dir}/${nsub}/log.txt)"
antsRegistrationSyN.sh \
    -d 3 \
    -f "${subject_dir}/tractoflow/${nsub}__T1w.nii.gz" \
    -m "${mni_dir}/MNI/mni_masked.nii.gz" \
    -t s \
    -o "${output_dir}/${nsub}/orig_space/transfo/2orig_" \
    > "${output_dir}/${nsub}/log.txt" 2>&1


    ## [ORIG-SPACE] Register all ROIs
for nroi in cn8_left cn8_right; do
    antsApplyTransforms \
    -d 3 \
    -i ${mni_dir}/MNI/${nroi}.nii.gz \
    -r ${subject_dir}/tractoflow/${nsub}__T1w.nii.gz \
    -t ${output_dir}/${nsub}/orig_space/transfo/2orig_1Warp.nii.gz \
    -t ${output_dir}/${nsub}/orig_space/transfo/2orig_0GenericAffine.mat \
    -o ${orig_rois_dir}/${nsub}_${nroi}_orig.nii.gz

    scil_volume_math lower_threshold_eq ${orig_rois_dir}/${nsub}_${nroi}_orig.nii.gz 0.5 ${orig_rois_dir}/${nsub}_${nroi}_orig.nii.gz --data_type int16 -f
done
    
