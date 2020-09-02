#!/usr/bin/env bash

# Script for checking the results of RGB plot clipping

# Define some counts that we expect (the number of sub-folders plus the top folder)
EXPECTED_NUM_FOLDERS=16
EXPECTED_NUM_CLIP_TIF=16

# What folder are we looking in for outputs
if [[ ! "${1}" == "" ]]; then
  TARGET_FOLDER="${1}"
else
  TARGET_FOLDER="./outputs"
fi

# Count the number of folders
FOLDER_LIST=(`find "${TARGET_FOLDER}/" -maxdepth 1 -type d`)
if [[ "${#FOLDER_LIST[@]}" == "${EXPECTED_NUM_FOLDERS}" ]]; then
  echo "Found expected number of folders: ${EXPECTED_NUM_FOLDERS}"
else
  echo "Expected ${EXPECTED_NUM_FOLDERS} folders and found ${#FOLDER_LIST[@]}"
  for i in $(seq 0 $(( ${#FOLDER_LIST[@]} - 1 )))
  do
    echo "$(( ${i} + 1 )): ${FOLDER_LIST[$i]}"
  done
  exit 10
fi

# Check the expected number of image mask files
EXPECTED_CLIP=(`find "${TARGET_FOLDER}/" -type f | grep 'orthomosaic\.tif'`)
if [[ "${#EXPECTED_CLIP[@]}" == "${EXPECTED_NUM_CLIP_TIF}" ]]; then
  echo "Found expected number of orthomosaic.tif files: ${EXPECTED_NUM_CLIP_TIF}"
else
  echo "Expected ${EXPECTED_NUM_CLIP_TIF} orthomosaic.tif files but found ${#EXPECTED_CLIP[@]}"
  for i in $(seq 0 $(( ${#EXPECTED_CLIP[@]} - 1 )))
  do
    echo "$(( ${i} + 1 )): ${EXPECTED_CLIP[$i]}"
  done
  exit 20
fi
