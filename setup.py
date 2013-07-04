from distutils.core import setup
import sys

sys.path.append('gammasens')
import gammasens


setup(name='gammasens',
      version='1.0',
      author='Karl Kosack',
      author_email='karl.kosack@cea.fr',
#      url='',
#      download_url='',
      description='Python-based Sensitivity calculator for Ground-based Gamma-Ray telescopes',
      long_description=gammasens.__doc__,
      package_dir={'': 'gammasens'},
      py_modules=['gammasens','sensitivity','inputs','spectra'],
      provides=['gammasens'],
      keywords='IACT gamma sensitivity',
      license='',
      # classifiers=['Development Status :: 5 - Production/Stable',
      #              'Intended Audience :: Developers',
      #              'Natural Language :: English',
      #              'Operating System :: OS Independent',
      #              'Programming Language :: Python :: 2',
      #              'License :: OSI Approved :: GNU Library or Lesser General Public License (LGPL)',
      #              'License :: OSI Approved :: GNU Affero General Public License v3',
      #              'Topic :: Internet',
      #              'Topic :: Internet :: WWW/HTTP',
      #              'Topic :: Scientific/Engineering :: GIS',
      #             ],
     )
