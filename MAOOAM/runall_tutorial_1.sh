set -e
python generate_nature_run.py
python generate_nature_run.py freerun
python generate_nature_Mhist.py
python generate_nature_QR.py
python generate_nature_LEs.py
python generate_observations.py
