let

overlay1 = self: super:
{
  sundials1 = self.callPackage ./CustomSundials { };
};

nixpkgs = builtins.fetchTarball {
    url = "https://github.com/NixOS/nixpkgs/archive/20.03.tar.gz";
    sha256 = "0182ys095dfx02vl2a20j1hz92dx3mfgz2a6fhn31bqlp1wa8hlq";
};

in

{ pkgs ? import nixpkgs { overlays = [ overlay1 ]; } }:

let

haskellPackages = pkgs.haskellPackages;

drv = haskellPackages.callPackage ./default.nix {
  klu = pkgs.suitesparse;
  suitesparseconfig = pkgs.suitesparse;
  sundials_arkode = pkgs.sundials1;
  sundials_cvode = pkgs.sundials1;
  sundials_sunlinsolklu = pkgs.sundials1;
  sundials_sunmatrixsparse = pkgs.sundials1;
};

in

drv
