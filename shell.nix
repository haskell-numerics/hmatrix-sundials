let

overlay1 = self: super:
{
  sundials1 = self.callPackage ./CustomSundials { };
};

in

{ nixpkgs ? import <nixpkgs> { overlays = [ overlay1 ]; }, compiler ? "default", doBenchmark ? false }:

let

  inherit (nixpkgs) pkgs;

  f = { mkDerivation, aeson, aeson-pretty, base, bytestring
      , cassava, clock, containers, deepseq, filepath, ghc-prim, hmatrix
      , inline-c, katip, klu, mtl, optparse-applicative, split, stdenv
      , suitesparseconfig, sundials_arkode, sundials_cvode
      , sundials_sunlinsolklu, sundials_sunmatrixsparse, tasty
      , tasty-golden, tasty-hunit, template-haskell, text, vector
      }:
      mkDerivation {
        pname = "hmatrix-sundials";
        version = "0.20.0.0";
        src = ./.;
        libraryHaskellDepends = [
          aeson base bytestring containers deepseq ghc-prim hmatrix inline-c
          katip mtl split template-haskell text vector
        ];
        librarySystemDepends = [
          klu suitesparseconfig sundials_arkode sundials_cvode
          sundials_sunlinsolklu sundials_sunmatrixsparse
        ];
        testHaskellDepends = [
          aeson aeson-pretty base bytestring filepath hmatrix katip tasty
          tasty-golden tasty-hunit vector
        ];
        benchmarkHaskellDepends = [
          base bytestring cassava clock hmatrix optparse-applicative
        ];
        description = "hmatrix interface to sundials";
        license = stdenv.lib.licenses.bsd3;
      };

  haskellPackages = if compiler == "default"
                       then pkgs.haskellPackages
                       else pkgs.haskell.packages.${compiler};

  variant = if doBenchmark then pkgs.haskell.lib.doBenchmark else pkgs.lib.id;

  drv = variant (haskellPackages.callPackage f { klu = pkgs.suitesparse; suitesparseconfig = pkgs.suitesparse; sundials_arkode = pkgs.sundials1; sundials_cvode = pkgs.sundials1; sundials_sunlinsolklu = pkgs.sundials1; sundials_sunmatrixsparse = pkgs.sundials1; });

in

  if pkgs.lib.inNixShell then drv.env else drv
