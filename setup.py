from setuptools import setup, find_packages



setup(
            name="fluned-1",
            version="0.0.1",
            author="Marco De Pietri",
            author_email="mdepietri@ind.uned.es",
            description="post processing of fluned simuations",
            packages=find_packages(where="src"),
            package_dir={'': 'src'},
            entry_points = {'console_scripts': ['flunedPost = flunedPost.flunedPost:main']},
            install_requires=[
            "h5py",
            "numpy>=1.21,<2", 
            "vtk",
            "pyvista",
            ],
            python_requires=">=3.11",  # Specify supported Python versions
            )
         
