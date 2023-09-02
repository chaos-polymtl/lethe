from setuptools import setup

REQUIRES = []
with open('requirements.txt') as f:
    for line in f:
        line, _, _ = line.partition('#')
        line = line.strip()
        if not line or line.startswith('setuptools'):
            continue
        else:
            REQUIRES.append(line)

setup(
    name='lethe_pyvista_tools',
    version='0.1',
    license='LGPL-2.1',
    packages=['lethe_pyvista_tools'],
    author='Lethe',
    url="https://lethe-cfd.github.io/lethe/index.html",
    author_email='',
    python_requires='>=3.8',
    install_requires=REQUIRES,
    description='Lethe post-processing tools using PyVista'
)