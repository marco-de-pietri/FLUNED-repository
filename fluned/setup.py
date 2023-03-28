import setuptools



setuptools.setup(
            name="fluned",
            version="0.9",
            author="Marco De Pietri",
            author_email="mdepietri@ind.uned.es",
            description="pre processor of fluned simuations",
            packages = setuptools.find_packages(),
            entry_points = {'console_scripts': ['fluned = fluned.fluned:main']})
         
