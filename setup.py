from distutils.core import setup
import sys
from gammasens import __version__

setup(name='gammasens',
      version=str(__version__),
      author='Karl Kosack',
      author_email='karl.kosack@cea.fr',
#      url='',
#      download_url='',
      description='Python-based Sensitivity calculator for Ground-based Gamma-Ray telescopes',
      long_description='',
      packages=['gammasens'],
      keywords='IACT gamma sensitivity',
#      license='',
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
