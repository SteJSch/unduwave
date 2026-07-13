# unduwave
Undupy is a python wrapper providing an API to use Wave and Undumag with python.

### Installation
```
pip3 install unduwave
```

### Get started
Example for the usage of undupy.

installation:

git clone https://github
install requirements with pip3 -r requirements.txt
if you get error: ImportError: DLL load failes while importing _cext,
go to
https://learn.microsoft.com/en-us/cpp/windows/latest-supported-vc-redist?view=msvc-170#latest-microsoft-visual-c-redistributable-version

and download redistributable package and install 


```Python
from apy import my_lib

lib = my_lib()
res=lib.mult(6)
print(res)
```
