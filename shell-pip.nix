{ pkgs ? import <nixpkgs> { } }:

let
  pythonEnv = pkgs.python310.withPackages(ps: with ps; [ scipy numpy jupyterlab setuptools build ] );
in

pkgs.mkShell {
  packages = [
    pythonEnv
  ];
}
