Combination of MachineLearningHEP (of which I'm co-developer) and hipe4ml package for D2H hadronisation analyses (Lc, Ds, and Bs). Moved into a separate repository for testing purposes (it should be seen as a branch of the MLHEP repo). **All credits go to the original authors.**

## Installation requirements
For hipe4ml
```
pip3 install hipe4ml
pip3 install uproot
pip3 install xxhash
pip3 install lz4
```

For MLHEP:
```
pip3 install PyYaml
pip3 install numba
pip3 install --no-binary :all: root_numpy
```
Extra for MLHEP (commented out at the moment):
```
pip3 install keras
pip3 install tensorflow
```

The complete installation requirements for MLHEP are:
```
  install_requires=[ "numpy==1.17.4", "pandas==0.25.3", "scipy==1.2.1", "matplotlib==3.1.2",
                     "seaborn==0.9.0", "uproot==3.11.1", "scikit-learn==0.22.1", "xgboost==0.90",
                     "keras==2.3.1", "tensorflow==1.14.0", "PyYaml==5.1", "pylint==2.4.3",
                     "twisted==19.2.0", "klein==17.10.0", "Jinja2==2.10.3", "numba==0.48.0",
                     "pyarrow==0.13.0", "lz4==2.1.10", "hipe4ml", "xxhash"],
```
Most of them are however already installed indirectly by the hipe4ml package