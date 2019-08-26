{ pkgs ? import <nixpkgs> {}, doBenchmark ? false }:

let

my-hmatrix-sundials   = pkgs.haskellPackages.callPackage ./default.nix { };

haskellDeps = ps: with ps; [
  my-hmatrix-sundials
];

in

  pkgs.stdenv.mkDerivation {
  name = "env";
  buildInputs = [
    (pkgs.haskellPackages.ghcWithPackages haskellDeps)
  ];
}

