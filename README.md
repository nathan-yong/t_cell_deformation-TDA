# Set Up Python Analysis Program (Linux)
### Python Version 3.12

If you don't have `venv`. Note that if you're default python version isn't 3.12. You may have to install python 3.12 and ensure venv is using python 3.12.
```
sudo apt update
sudo apt install python3-venv
```
Go to python directory `t_cell_deformation/python_analysis` 

Start `venv`. If your default version is not 3.12, ensure you are using `python3.12`.
```
python3 -m venv .venv
```
```
source .venv/bin/activate
```
Install dependencies
```
pip install -r requirements.txt
```