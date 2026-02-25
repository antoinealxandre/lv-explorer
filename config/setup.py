#!/usr/bin/env python3
"""
Setup script pour LV Explorer

DEPRECATED: project metadata has moved to /pyproject.toml at the repository root.
This file is kept for backwards compatibility only.
"""

from setuptools import setup, find_packages
import os


def read_requirements():
    """Lit le fichier requirements.txt"""
    with open('requirements.txt', 'r') as f:
        return [line.strip() for line in f if line.strip() and not line.startswith('#')]


def read_readme():
    """Lit le README"""
    if os.path.exists('README.md'):
        with open('README.md', 'r', encoding='utf-8') as f:
            return f.read()
    return ''


setup(
    name='lv-explorer',
    version='0.1.0',
    author='Cardiac Analysis Lab',
    author_email='your.email@example.com',
    description='Application unifiée d\'analyse et de simulation ventriculaire',
    long_description=read_readme(),
    long_description_content_type='text/markdown',
    url='https://github.com/yourusername/lv-explorer',
    packages=find_packages(),
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Healthcare Industry',
        'Topic :: Scientific/Engineering :: Medical Science Apps.',
        'Topic :: Scientific/Engineering :: Visualization',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.8',
    install_requires=read_requirements(),
    entry_points={
        'console_scripts': [
            'lv-explorer=lv_explorer.ui.main_window:launch_app',
        ],
    },
    include_package_data=True,
    zip_safe=False,
)