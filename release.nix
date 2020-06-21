let

sundialsOverlay = self: super:
{
  sundials1 = self.callPackage ./CustomSundials { };
};

myHaskellPackageOverlay = self: super: {
  myHaskellPackages = super.haskellPackages.override {
    overrides = hself: hsuper: rec {

   tasty-golden =
        let newTastyGoldenSrc = builtins.fetchGit {
          url = "https://github.com/feuerbach/tasty-golden.git";
          rev = "500d828f683bd84eae89beaa43f577a64e41faf5";
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
