from setuptools import setup

def readme():
    with open('README.rst') as f:
        return f.read()

setup(name='maclab',
      version='0.1',
      description='Macroecology tools',
      url='http://github.com/cwarecsiro/maclab',
      author='Chris Ware',
      author_email='chris.ware@csiro.au',
      license='NA',
      packages=['maclab'],
      install_requires=[
          'rasterio',
          'warnings',
          'jenkspy',
          'numpy',
          'operator',
          'functools'
      ],
      test_suite='nose.collector',
      tests_require=['nose'],
      zip_safe=False)