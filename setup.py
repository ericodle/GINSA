from setuptools import setup, find_packages

setup(
    name='ginsa',
    version='0.1.0',
    packages=find_packages('src'),
    package_dir={'': 'src'},
    entry_points={
        'console_scripts': [
            'GINSA_cli = GINSA.GINSA_cli:main',
            'GINSA_gui = GINSA.GINSA_gui:main',
            'suffix_adder = GINSA.suffix_adder:main',
        ],
    },
    
    install_requires=[
            "anyio==4.2.0",
            "biopython==1.83",
            "certifi==2024.2.2",
            "charset-normalizer==3.3.2",
            "contourpy==1.2.0",
            "cycler==0.12.1",
            "docutils==0.20.1",
            "fonttools==4.47.2",
            "h11==0.14.0",
            "httpcore==1.0.2",
            "httpx==0.26.0",
            "idna==3.6",
            "importlib-metadata==7.0.1",
            "iniconfig==2.0.0",
            "jaraco.classes==3.3.0",
            "keyring==24.3.0",
            "kiwisolver==1.4.5",
            "markdown-it-py==3.0.0",
            "matplotlib==3.8.2",
            "mdurl==0.1.2",
            "more-itertools==10.2.0",
            "nh3==0.2.15",
            "numpy==1.26.3",
            "packaging==23.2",
            "pandas==2.2.0",
            "pillow==10.2.0",
            "pkginfo==1.9.6",
            "pluggy==1.4.0",
            "psutil==5.9.8",
            "Pygments==2.17.2",
            "pyparsing==3.1.1",
            "pytest==8.0.0",
            "python-dateutil==2.8.2",
            "pytz==2024.1",
            "readme-renderer==42.0",
            "requests==2.31.0",
            "requests-toolbelt==1.0.0",
            "rfc3986==2.0.0",
            "rich==13.7.0",
            "setuptools==69.0.3",
            "six==1.16.0",
            "sniffio==1.3.0",
            "tk==0.1.0",
            "twine==4.0.2",
            "tzdata==2023.4",
            "urllib3==2.2.0",
            "wget==3.2",
            "wheel==0.42.0",
            "zipp==3.17.0"
    ],
)
