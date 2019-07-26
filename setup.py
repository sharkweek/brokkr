from setuptools import setup

setup(
    name='brokkr',
    version='0.0.1',
    description="mechanical engineering toolbox",
    author='Andy Perez',
    packages=[
        'brokkr/bolted_joints',
        'brokkr/bmath',
        'brokkr/composites',
    ],
    license="MIT",
    long_description=open('README.md').read(),
    url='https://github.com/sharkweek/brokkr',
    keywords=[
        'composites',
        'mechanics',
        'statics',
        'materials',
        'structures',
    ],
    classifiers=[
        "Programming Language :: Python",
        "Programming Language :: Python :: 3.5+",
        "License :: MIT",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering",
        "Development Status :: 3 - Alpha",
    ],
    install_requires=[
        'numpy',
        'pandas',
        'unyt>=2.0.0',
    ],
    python_requires='>=3.5'
)
