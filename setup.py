import setuptools
import csst

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name='csst',
    version=csst.__version__,
    author='CSST Team',
    author_email='bozhang@nao.cas.cn',
    description='The CSST pipeline',  # short description
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://csst.readthedocs.io',
    project_urls={
        'Source': 'http://github.com/csster/csst',
    },
    packages=setuptools.find_packages(),
    license='MIT',
    classifiers=["Development Status :: 5 - Production/Stable",
                 "Intended Audience :: Science/Research",
                 "License :: OSI Approved :: MIT License",
                 "Operating System :: OS Independent",
                 "Programming Language :: Python :: 3.7",
                 "Topic :: Scientific/Engineering :: Physics",
                 "Topic :: Scientific/Engineering :: Astronomy"],
    package_dir={'csst': 'csst'},
    # include_package_data=True,
    package_data={"": ["LICENSE", "README.md"],
                  "csst": ["msc/deepcr_model/*",
                           "msc/pos_calib_config/*",
                           "msc/flux_calib_config/*",
                           "msc/phot_config/*",
                           ]},
    # install_requires=['sphinx>=4.2.0',
    #                   'numpy>=1.22.0',
    #                   'scipy', 'matplotlib',
    #                   'astropy', 'healpy', 'ccdproc', 'deepCR'],
    python_requires='>=3.7',
)
