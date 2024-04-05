## This code remove facial structure from T1w brain MRI images
#cd your/main/path
code_path=2024_brsc_fc_project/code/
rawdata_path=2024_brsc_fc_data/rawdata/
preproc_path=2024_brsc_fc_data/derivatives/preprocessing/

IDs=(A006 A007 A008 A009 A011 A013 A023 A024 A025 A026 A027 A028 A029 A030 A032 A033 A034 A036 A037 A028 A039 A043 A044 A045 A046 A048 A050 A051 A052 P030 P099)
# 1 Check if the raw files exists
echo "This dataset is composent of ${#IDs[@]} participants"
for ID in "${IDs[@]}"; do
if [ ! -d "$rawdata_path/sub-$ID" ]; then
  echo "sub-"$ID" does not exist."
fi

#Check anat folder and files:
if [ ! -d "$rawdata_path/sub-"$ID"/anat/" ]; then
  echo "sub-"$ID"/anat/ does not exist."
fi
if [ ! -e "$rawdata_path/sub-$ID/anat/sub-"$ID"_T1w.json" ]; then
  echo "sub-"$ID"/anat/sub-"$ID"_T1w.json does not exist."
fi
if [ ! -e "$rawdata_path/sub-$ID/anat/sub-"$ID"_T1w.nii.gz" ]; then
  echo "sub-"$ID"/anat/sub-"$ID"_T1w.nii.gz does not exist."
fi

#Check func folder and files:
if [ ! -d "$rawdata_path/sub-$ID/func/" ]; then
  echo "sub-$ID/func/ does not exist."
fi
if [ ! -e "$rawdata_path/sub-$ID/func/sub-${ID}_task-rest_bold.nii.gz" ]; then
  echo "sub-"$ID"/func/sub-"$ID"_task-rest_bold.nii.gz does not exist."
fi

if [ ! -e "$rawdata_path/sub-$ID/func/sub-${ID}_task-rest_bold.json" ]; then
  echo "sub-"$ID"/func/sub-"$ID"_task-rest_bold.json does not exist."
fi

done

# 2 Deface anat
#we followd the recomandation of https://spine-generic.readthedocs.io/data-acquisition.html



