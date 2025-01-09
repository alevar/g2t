from setuptools import setup, find_packages

setup(
    name='g2t',
    version='0.1',
    packages=find_packages(),
    install_requires=[],
    entry_points={
        'console_scripts': [
            'g2t=g2t.core:main',
        ],
    },
)
