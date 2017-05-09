from setuptools import setup

def readme():
    with open('README.rst', 'r') as f:
        return f.read()


setup(
    name='fits_corr',
    version='0.1.1a1',
    description='Correlation & Similarity Analysis for FITS data files.',
    long_description=readme(),
    classifiers=[
        # Project Maturity
        'Development Status :: 3 - Alpha',
        
        # Licensing
        'License :: OSI Approved :: MIT License',
        
        # Python Version Support
        'Programming Language :: Python :: 2.7',
        
        # Intended Audience
        'Intended Audience :: Developers',
        
        # Categories
        'Topic :: Scientific/Engineering :: Information Analysis',
        'Topic :: Scientific/Engineering :: Astronomy',
    ],
    url='https://github.com/UCBerkeleySETI/fits-corr,
    author='Pragaash Ponnusamy',
    author_email='pragaash.io@gmail.com',
    license='MIT',
    packages=['fits_corr'],
    entry_points={
        'console_scripts': [
            'fits_corr = fits_corr.__main__:main'
        ]
    },
    install_requires=[
        'numpy>=1.11.3', 
        'scipy>=0.19.0', 
        'scikit-image>=0.12.3',
        'fitsio>=0.9.11',
        'progressbar2>=3.16.0',
    ],
    include_package_data=True,
    zip_safe=False
)