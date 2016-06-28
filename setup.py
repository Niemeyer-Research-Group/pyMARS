from setuptools import setup, find_packages

import os
current_dir=os.getcwd()
setup(name='pyMARS',
      version='0.3',
      description='Python-based (chemical kinetic) Model Automatic \
      Reduction Software (MARS), which consists of multiple techniques for \
      reducing the size and complexity of detailed chemical kinetic models.',
      url='http://github.com/kyleniemeyer/pyMARS',
      author='',
      author_email='claytonp@oregonstate.edu',
      license='MIT',
      packages=find_packages(),
      zip_safe=False,
      entry_points={
        'console_scripts': [
                    'pyMARS= pyMARS.__main__:main']},
      install_requires=['argparse', 'h5py' , 'numpy',
                    'Cantera>=2.3.0a2', 'matplotlib']

        )
