with import <nixpkgs> {};
with python310Packages;

buildPythonPackage rec {
  name = "m3l";
  src = ./src;
  propagatedBuildInputs = [ setuptools numpy scipy jupyterlab ];
}
