lisn_utils
==========

Overview
--------
The *lisn_utils* package contains several utilities for processing LISN data.

For more information about LISN (Low-Latitude Ionospheric Sensor Network) project visit http://lisn.igp.gob.pe

Installation
------------
Before install the package you must install the following OS system dependences (Use yum or apt-get depending of your LINUX distribution)
* python-devel
* python-setuptools
* gcc
* zlib-devel
* python-matplotlib

The following commands will install these dependencies in a Debian/Ubuntu based system:
```
$ sudo apt-get install python-dev python-setuptools gcc zlib1g-dev python-matplotlib
```

Finally these commands will install the package as well as the rest of the dependencies:
```
$ git clone https://github.com/jucares/lisn_utils.git
$ cd lisn_utils
$ sudo python setup.py install
```

Basic Usage
-----------
OK. Here's some exmaples showing how to read and parse rinex, tec and scintillation files.

Example 1: Parse a rinex file, calculate TEC (total electron content) and plot equivalent TEC vs Time.
```
>>> from lisn_utils import RNXData
>>> rnx = RNXData('path/to/rinex_filename')
>>> rnx.calctec()
>>> rxn.plot('epoch', 'eqTEC')
```
Example 2: Parse a TEC file and plot Slant TEC vs Time.
```
>>> from lisn_utils import TECData
>>> tec = TECData('path/to/tec_filename')
>>> tec.plot('epoch', 'sTEC')
```
Example 3: Parse a Scintillation file and plot S4 index vs Time.
```
>>> from lisn_utils import S4Data
>>> s4 = S4Data('path/to/s4_filename')
>>> s4.plot('epoch', 's4')
```
