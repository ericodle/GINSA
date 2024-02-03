from setuptools import setup, find_packages

setup(
    name='ginsa',
    version='0.1.0',
    packages=find_packages('src'),
    package_dir={'': 'src'},
    entry_points={
        'console_scripts': [
            'ginsa-cli = ginsa.ginsa_cli:main',
            'ginsa-gui = ginsa.ginsa_gui:main',
        ],
    },
    
    install_requires=[
          "biopython==1.83",
          "certifi==2024.2.2",
          "charset-normalizer==3.3.2",
          "contourpy==1.2.0",
          "cycler==0.12.1",
          "fonttools==4.47.2",
          "idna==3.6",
          "kiwisolver==1.4.5",
          "matplotlib==3.8.2",
          "numpy==1.26.3",
          "packaging==23.2",
          "pandas==2.2.0",
          "pillow==10.2.0",
          "psutil==5.9.8",
          "pyparsing==3.1.1",
          "python-dateutil==2.8.2",
          "pytz==2024.1",
          "requests==2.31.0",
          "six==1.16.0",
          "tk==0.1.0",
          "tzdata==2023.4",
          "urllib3==2.2.0",
          "wget==3.2"
    ],
)
