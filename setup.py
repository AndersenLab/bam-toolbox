from setuptools import setup
import glob
from bam import __version__
with open('requirements.txt') as f:
    required = f.read().splitlines()

setup(name='bam-toolbox',
      version=__version__,
      packages=['bam','bam.output'],
      description='Tools for working with BAM files',
      url='https://github.com/AndersenLab/bam-toolbox',
      author='Daniel E. Cook',
      author_email='danielecook@gmail.com',
      license='MIT',
      install_requires=required,
      entry_points={
            'console_scripts': [
                  'bam = bam.cli:main'
            ]
      },
      zip_safe=False)
