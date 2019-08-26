# { nixpkgs ? import <nixpkgs> {}, compiler ? "default", doBenchmark ? false }:

# let

#   inherit (nixpkgs) pkgs;

#   f = { mkDerivation, base, bytestring, cassava, Chart
#       , Chart-diagrams, clock, containers, deepseq, diagrams-cairo
#       , diagrams-lib, diagrams-rasterific, hmatrix, hspec, inline-c, lens
#       , Naperian, optparse-applicative, plots, split, stdenv
#       , sundials, template-haskell, vector
#       }:
#       mkDerivation {
#         pname = "hmatrix-sundials";
#         version = "0.20.0.0";
#         src = ./.;
#         isLibrary = true;
#         isExecutable = true;
#         libraryHaskellDepends = [
#           base containers deepseq hmatrix inline-c split template-haskell
#           vector
#         ];
#         librarySystemDepends = [ sundials ];
#         executableHaskellDepends = [
#           base Chart Chart-diagrams diagrams-cairo diagrams-lib hmatrix
#           Naperian vector
#         ];
#         testHaskellDepends = [
#           base containers diagrams-lib diagrams-rasterific hmatrix hspec
#           inline-c lens plots split template-haskell vector
#         ];
#         testSystemDepends = [ sundials ];
#         benchmarkHaskellDepends = [
#           base bytestring cassava clock hmatrix optparse-applicative
#         ];
#         homepage = "https://github.com/idontgetoutmuch/hmatrix/tree/sundials";
#         description = "hmatrix interface to sundials";
#         license = stdenv.lib.licenses.bsd3;
#       };

#   haskellPackages = if compiler == "default"
#                        then pkgs.haskellPackages
#                        else pkgs.haskell.packages.${compiler};

#   variant = if doBenchmark then pkgs.haskell.lib.doBenchmark else pkgs.lib.id;

# drv = variant (haskellPackages.callPackage f { });

# in

#   if pkgs.lib.inNixShell then drv.env else drv

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

