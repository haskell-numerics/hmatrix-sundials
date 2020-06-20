let

overlay1 = self: super:
{
  sundials1 = self.callPackage ./CustomSundials { };
};

in

{ nixpkgs ? import <nixpkgs> { overlays = [ overlay1 ]; } }:

let

haskellPackages = nixpkgs.haskellPackages;

drv = haskellPackages.callPackage ./default.nix { klu = nixpkgs.suitesparse; suitesparseconfig = nixpkgs.suitesparse; sundials_arkode = nixpkgs.sundials1; sundials_cvode = nixpkgs.sundials1; sundials_sunlinsolklu = nixpkgs.sundials1; sundials_sunmatrixsparse = nixpkgs.sundials1; };

in

drv.env
