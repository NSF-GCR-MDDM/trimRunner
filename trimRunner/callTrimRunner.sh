#!/bin/bash
# Run them all in parallel
python3 trimRunner.py config/16O_full.jsonc &
python3 trimRunner.py config/24Mg_full.jsonc &
python3 trimRunner.py config/28Si_full.jsonc &
python3 trimRunner.py config/56Fe_full.jsonc &

# Wait for all background jobs to finish
wait
echo "All done."