try:
	from setuptools import setup
except ImportError:
	from disutils.core import setup

config = {
	'name': 'SeqDB',
	'version': '0.2',
	'description': 'A simple program to organize sequencing files',
	'packages': ['bin', 'Sequence Files'],
	'install_requires': ['nose', 'BioPython', 'PyQt4'],
	'scripts':'main.py',
	'author': 'Kiet T. Phong',
	'author_email': 'tuankietphong@gmail.com',
	'long_description': 'README.md'
}

setup(**config)
