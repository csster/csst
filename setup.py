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
    url='http://github.com/csster/csst',
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
    include_package_data=False,
    package_data={"": ["LICENSE", "README.md"]},
    requires=['numpy', 'scipy', 'astropy'],
    python_requires='>=3.9.7',
)
