#!/bin/bash
# Run them all in parallel
python3 trimRunner.py config/218At.jsonc &
python3 trimRunner.py config/218Po.jsonc &
python3 trimRunner.py config/220Rn.jsonc &
python3 trimRunner.py config/226Rn.jsonc &
python3 trimRunner.py config/230Th.jsonc &
python3 trimRunner.py config/234Pa.jsonc &
python3 trimRunner.py config/234Th.jsonc &
python3 trimRunner.py config/234U.jsonc &

# Wait for all background jobs to finish
wait
echo "All done."