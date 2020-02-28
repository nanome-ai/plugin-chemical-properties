import pathlib
from setuptools import find_packages, setup

README = (pathlib.Path(__file__).parent / "README.md").read_text()

setup(
	name='nanome-chemical-properties',
	packages=find_packages(),
	version='0.3.4',
	license='MIT',
	description='A Nanome plugin to display chemical properties for a complex using rdkit',
	long_description=README,
    long_description_content_type="text/markdown",
	author='Nanome',
	author_email='hello@nanome.ai',
	url='https://github.com/nanome-ai/plugin-chemical-properties',
	platforms="any",
	keywords=['virtual-reality', 'chemistry', 'python', 'api', 'plugin'],
	install_requires=['nanome', 'cairosvg'],
	entry_points={"console_scripts": ["nanome-chemical-properties = nanome_chemical_properties.ChemicalProperties:main"]},
	classifiers=[
		'Development Status :: 3 - Alpha',

		'Intended Audience :: Science/Research',
		'Topic :: Scientific/Engineering :: Chemistry',

		'License :: OSI Approved :: MIT License',

		'Programming Language :: Python :: 2.7',
		'Programming Language :: Python :: 3.5',
		'Programming Language :: Python :: 3.6',
		'Programming Language :: Python :: 3.7',
	],
	package_data={
        "nanome_chemical_properties": [
            "menus/json/*",
			"icons/*"
        ]
	},
)
