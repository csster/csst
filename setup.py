import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name='csst',
    version='0.0.1',
    author='Chao Liu',
    author_email='liuchao@nao.cas.cn',
    description='The CSST pipeline',  # short description
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='http://github.com/hypergravity/csst',
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
    requires=['numpy', 'scipy', 'astropy', 'joblib']
)
