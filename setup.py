from setuptools import setup

description = 'Python implementations of Gene Set Enrichment Analysis(GSEA).'

with open('README.md') as f:
    long_description = f.read()

setup(
    name='gsea',
    version='0.0.1',
    author='Martin Vuk',
    description=description,
    long_description=long_description,
    license='MIT',
    keywords='bioinformatics',
    setup_requires = ['pytest-runner'],
    tests_require = ['pytest'],
    install_requires=['numpy>=1.12.0'],
    packages=['gsea'],
    scripts= ['bin/gsea'],
    classifiers=[
        'Intended Audience :: Science/Research',
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3.6',
    ],
)
