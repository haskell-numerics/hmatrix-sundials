{ mkDerivation, base, bytestring, cassava, Chart, Chart-diagrams
, clock, containers, deepseq, diagrams-cairo, diagrams-lib
, diagrams-rasterific, hmatrix, hspec, inline-c, lens
, optparse-applicative, plots, split, stdenv, sundials
, template-haskell, vector
}:
mkDerivation {
  pname = "hmatrix-sundials";
  version = "0.20.1.0";
  src = ./.;
  isLibrary = true;
  isExecutable = true;
  libraryHaskellDepends = [
    base containers deepseq hmatrix inline-c split template-haskell
    vector
  ];
  librarySystemDepends = [ sundials ];
  executableHaskellDepends = [
    base Chart Chart-diagrams diagrams-cairo diagrams-lib hmatrix
    vector
  ];
  testHaskellDepends = [
    base containers diagrams-lib diagrams-rasterific hmatrix hspec
    inline-c lens plots split template-haskell vector
  ];
  testSystemDepends = [ sundials ];
  benchmarkHaskellDepends = [
    base bytestring cassava clock hmatrix optparse-applicative
  ];
  homepage = "https://github.com/haskell-numerics/hmatrix-sundials";
  description = "hmatrix interface to sundials";
  license = stdenv.lib.licenses.bsd3;
}
