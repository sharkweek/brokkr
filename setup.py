#from distutils.core import setup
from setuptools import setup

setup(
    name='mechpy',
    version='0.1',
    description="mechanical engineering toolbox",
    author='Neal Gordon, Andy Perez',
    packages=[
        'composites',
        'bolted_joints'
    ],
    package_dir={'': 'lib'},
    license="MIT",
    long_description=open('README.md').read(),
    url='https://github.com/sharkweek/mechpy',
    keywords=[
        'composites',
        'mechanics',
        'statics',
        'materials'
    ],
    classifiers=[
        "Programming Language :: Python",
        "Programming Language :: Python :: 3.5+",
        "License :: MIT",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering",
        "Development Status :: 3 - Alpha"
    ],
    install_requires=[
        'numpy',
        'matplotlib',
        'scipy',
        'sympy',
        'pint'
    ]
)
