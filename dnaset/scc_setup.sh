module load python3/3.10.5 pytorch bedtools ucscutils gcc
[ ! -d "env" ] && python -m venv env
source env/bin/activate
pip install --upgrade pip
pip install -r requirements.txt

