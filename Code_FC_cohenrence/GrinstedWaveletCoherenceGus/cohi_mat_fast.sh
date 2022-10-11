#!/bin/bash
# example: cohi_mat_fast(inputFileName, saveFileName, TR, lowPass,highPass)
#          TR is the TR value for scans
#          [lowPass, highPass] is the filtering band
    


# Exports all environment variables
#$ -V

# Execute job from current directory
#$ -cwd

# IMPORTANT - deleting the next line will cause your job to fail

echo "cohi_mat_fast('$1','$2',$3,[$4,$5])"
matlab -nodisplay -r "cohi_mat_fast('$1','$2',$3,[$4,$5])"
