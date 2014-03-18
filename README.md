NGTS_ZLP
========

The Zero Level Pipeline for the NGTS project

Development
-----------

To help develop this code, preferably inside a virtualenv: run

``` python
python setup.py develop
```

This will effectively create symlinks to the package contents to allow development without installing properly after each change.

To uninstall, run

``` python
python setup.py develop --uninstall
```

The python scripts to run the code will be updated during development.

Installation
------------

To install, run `python setup.py install` or to install directly from git:

```
python setup.py install git+https://github.com/tomlouden/NGTS_ZLP.git
```
