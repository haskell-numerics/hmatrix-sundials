let

sundialsOverlay = self: super:
{
  sundials1 = self.callPackage ./CustomSundials { };
};

myHaskellPackageOverlay = self: super: {
  myHaskellPackages = super.haskellPackages.override {
    overrides = hself: hsuper: rec {

   tasty-golden =
        let newTastyGoldenSrc = builtins.fetchTarball { url = "https://hackage.haskell.org/package/tasty-golden-2.3.3/tasty-golden-2.3.3.tar.gz";
          sha256 = "0wgcs4pqr30bp801cyrg6g551i7q0vjjmd9gmk5jy44fgdhb7kkl";
          };
            tg = hself.callCabal2nix "tasty-golden" newTastyGoldenSrc {};
          in
          super.haskell.lib.dontCheck tg;
      };
    };
};

nixpkgs = builtins.fetchTarball {
    url = "https://github.com/NixOS/nixpkgs/archive/20.03.tar.gz";
    sha256 = "0182ys095dfx02vl2a20j1hz92dx3mfgz2a6fhn31bqlp1wa8hlq";
};

in

{ pkgs ? import nixpkgs { overlays = [ sundialsOverlay myHaskellPackageOverlay ]; } }:

let

haskellPackages = pkgs.myHaskellPackages;

in

haskellPackages.callPackage ./default.nix {
  klu = pkgs.suitesparse;
  suitesparseconfig = pkgs.suitesparse;
  sundials_arkode = pkgs.sundials1;
  sundials_cvode = pkgs.sundials1;
  sundials_sunlinsolklu = pkgs.sundials1;
  sundials_sunmatrixsparse = pkgs.sundials1;
}
