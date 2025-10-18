from setuptools import setup, find_packages

setup(
    name='coherence-gravity-coupling',
    version='0.2.0',
    packages=find_packages(where='src') + ['.'],
    package_dir={'': 'src', 'cgc': '.'},
    py_modules=['__main__'],
    entry_points={
        'console_scripts': [
            'cgc=__main__:main',
        ],
    },
    install_requires=[
        'numpy>=1.21.0',
        'scipy>=1.7.0',
        'matplotlib>=3.4.0',
        'tqdm>=4.62.0',
        'pyamg>=4.2.0',
    ],
    python_requires='>=3.8',
    author='Coherence-Gravity Research Group',
    description='Framework for simulating coherence-induced modifications to gravitational coupling',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
)
