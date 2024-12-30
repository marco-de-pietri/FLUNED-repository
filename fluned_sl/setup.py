'''
setup script for fluned_sl
'''

from setuptools import setup, find_packages

# Read the requirements.txt file
with open('requirements.txt', 'r',encoding='utf-8') as f:
    requirements = f.read().splitlines()

setup(
    name="fluned_sl",
    version="0.1.0",
    description="radiological modelling of radioisotopes contained in activated water",
    author="Marco De Pietri",
    author_email="mdepietri@ind.uned.es",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    install_requires=[
            "numpy>=1.21,<2", 
            "iapws",
            ],
    package_data={
        'fluned_sl': ['*.total', ],
    },
    include_package_data=True,
    entry_points={
        'console_scripts': [
            'fluned_sl=fluned_sl.fluned_sl:main',
        ],
    }
)
