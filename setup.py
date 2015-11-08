from setuptools import setup


setup(name='bam-toolbox',
      version='0.0.1',
      packages=['bam', 'bam.eav'],
      description='Tools for working with BAM files',
      url='https://github.com/AndersenLab/bam-toolbox',
      author='Daniel Cook',
      author_email='danielecook@gmail.com',
      license='MIT',
      entry_points="""
      [console_scripts]
      bam = bam.bam:main
      """,
      install_requires=["docopt", "clint","pybedtools"],
      zip_safe=False)