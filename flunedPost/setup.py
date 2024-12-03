import setuptools



setuptools.setup(
            name="flunedPost",
            version="0.0.1",
            author="Marco De Pietri",
            author_email="mdepietri@ind.uned.es",
            description="post processing of fluned simuations",
            packages = setuptools.find_packages(),
            entry_points = {'console_scripts': ['flunedPost = flunedPost.flunedPost:main']},
            install_requires=[
            "h5py",
            "numpy>=1.21,<2", 
            "vtk",
            "pyvista",
            ],
            python_requires=">=3.11",  # Specify supported Python versions
            )
         
