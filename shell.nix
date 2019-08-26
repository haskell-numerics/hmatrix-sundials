{ pkgs ? import <nixpkgs> {}, doBenchmark ? false }:

let

hmatrix-sundials = pkgs.haskellPackages.callPackage ./default.nix { };

haskellDeps = ps: with ps; [
  hmatrix-sundials
];

in

  pkgs.stdenv.mkDerivation {
  name = "env";
  buildInputs = [
    (pkgs.haskellPackages.ghcWithPackages haskellDeps)
  ];
}

