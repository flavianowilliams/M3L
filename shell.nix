with import <nixpkgs> {};
with pkgs.python310Packages;

buildPythonPackage rec {
  name = "m3l";
  src = ./src;
  propagatedBuildInputs = [ pytest scipy numpy pkgs.libsndfile ];
}
