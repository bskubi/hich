#!/bin/bash

# Activate the environment
source activate juicer

# Define the URL and the target path
JAR_URL="https://s3.amazonaws.com/hicfiles.tc4ga.com/public/juicer/juicer_tools_1.22.01.jar"
TARGET_PATH="$CONDA_PREFIX/bin/juicer_tools.jar"

# Download the JAR file
wget -O $TARGET_PATH $JAR_URL

# Make it executable
chmod +x $TARGET_PATH

echo "Juicer tools downloaded to $TARGET_PATH and made executable."
