#! /usr/bin/python

'''
lisn_utils installer

:copyright:
    Juan C. Espinoza
:license:
    GNU Lesser General Public License, Version 3
    (http://www.gnu.org/copyleft/lesser.html)
'''

import os
from setuptools import setup
from setuptools.extension import Extension
from subprocess import call
from warnings import warn

src  = os.path.join('lisn_utils', 'src')
dst  = os.path.join('lisn_utils', 'bin')
path = os.path.dirname(os.path.abspath(__file__))
lib  = Extension('_crc', sources=[src+'/'+'crc.c',
                                      src+'/'+'crc_wrap.c'])
files = []

print '*** Compiling rnx2crx.c'
dum = call(['gcc', '%s/%s/rnx2crx.c' % (path, src),
            '-o', '%s/%s/rnx2crx' % (path, dst)])
if dum==0:
    files.append(dst+'/'+'rnx2crx')
else:
    warn('*** An error occur compiling rnx2crx. hatanaka compression will not be available')

print 'Compiling crx2rnx.c'
dum = call(['gcc', '%s/crx2rnx.c' % src, '-o', '%s/crx2rnx' % dst])

if dum==0:
    files.append(dst+'/'+'crx2rnx')
else:
    warn('*** An error occur compiling crx2rnx. hatanaka decompression will not be available')

data_files = [('/usr/local/bin/', files)]

if os.path.exists('/etc/lisn.cfg'):
    print '*** Copying old configuration file /etc/lisn.cfg.bk'
    os.rename('/etc/lisn.cfg', '/etc/lisn.cfg.bk')
data_files.append(('/etc/', ['lisn_utils/cfg/lisn.cfg']))

setup(
    name='lisn_utils',
    version='1.0',
    description='Several Utilities for processing LISN data',
    long_description='''
    Utilities for processing Data from the LISN project.

    For more information visit http://lisn.igp.gob.pe.
    ''',
    url='http://lisn.igp.gob.pe',
    author='Juan C. Espinoza',
    author_email='jucar.espinoza@gmail.com',
    license='GNU Lesser General Public License, Version 3 (LGPLv3)',
    platforms='Linux',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: ' + \
        'GNU Library or Lesser General Public License (LGPL)',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Physics',
    ],
    keywords=['GPS', 'RINEX', 'TEC'],
    packages=['lisn_utils'],
    ext_package='lisn_utils',
    ext_modules=[lib],
    zip_safe=False,
    data_files = data_files,
    package_data = {'': ['data/*.dat']},
    install_requires=[
        'numpy',
        'requests',
        ],
)

print '*** Please update User/password in /etc/lisn.cfg'
