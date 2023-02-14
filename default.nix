{ pgs ? import <nixpkgs> {} }:
	let
	my-python = pgs.python310Packages;
	in
	my-python.buildPythonPackage rec {
		name = "m3l";
		src = ./src;
		propagatedBuildInputs = with my-python; [
			build
			twine
			setuptools
			scipy
			jupyterlab
			];
	}
