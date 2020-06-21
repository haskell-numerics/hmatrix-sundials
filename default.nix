{ mkDerivation, aeson, aeson-pretty, base, bytestring, cassava
, clock, containers, deepseq, filepath, ghc-prim, hmatrix, inline-c
, katip, klu, mtl, optparse-applicative, split, stdenv
, suitesparseconfig, sundials_arkode, sundials_cvode
, sundials_sunlinsolklu, sundials_sunmatrixsparse, tasty
, tasty-golden, tasty-hunit, template-haskell, text, vector
}:
mkDerivation {
  pname = "hmatrix-sundials";
  version = "0.20.2.0";
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
  homepage = "https://github.com/haskell-numerics/hmatrix-sundials";
  description = "hmatrix interface to sundials";
  license = stdenv.lib.licenses.bsd3;
}
