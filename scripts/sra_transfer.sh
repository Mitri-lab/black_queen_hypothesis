#!/bin/bash -l
#SBATCH --cpus-per-task=1
#SBATCH --mem=2M
#SBATCH --time=03:00:00

HOST=''
USER=''
PASSWD=''
FILE=$1
REMOTE_DIR=''

# Use curl to upload the file
curl -T $FILE ftp://$USER:$PASSWD@$HOST/$REMOTE_DIR/

exit 0