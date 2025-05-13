{
  description = "Flake for github:formbio/laava";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-24.11";
    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs = {
    self,
    nixpkgs,
    flake-utils,
  }:
    flake-utils.lib.eachDefaultSystem (
      system: let
        pkgs = nixpkgs.legacyPackages.${system};

        pythonEnv = pkgs.python3.withPackages (ps:
          with ps; [
            biopython
            pandas
            pysam
          ]);

        rEnv = pkgs.rWrapper.override {
          packages = with pkgs.rPackages; [
            flextable
            gdtools
            rmarkdown
            systemfonts
            tidyverse
          ];
        };

        formbio-laava = pkgs.stdenv.mkDerivation {
          pname = "formbio-laava";
          version = "4.0.2";
          src = self;
          postPatch = ''
            substituteInPlace src/*.py \
              --replace '/usr/bin/env python3' '${pythonEnv}/bin/python'
            substituteInPlace src/*.R \
              --replace '/usr/bin/env Rscript' '${rEnv}/bin/Rscript'
          '';
          buildPhase = ":";
          installPhase = ''
            mkdir -p $out/opt/laava
            cp -r src/* $out/opt/laava/
            chmod +x $out/opt/laava/*.{py,R,sh}
          '';
        };

        runtimeEnv = with pkgs; [bash busybox formbio-laava pythonEnv rEnv];
      in {
        formatter = pkgs.alejandra;

        devShells.default = pkgs.mkShell {
          packages = with pkgs;
            [
              git
              gnumake
              lazygit
              nextflow
            ]
            ++ runtimeEnv;
        };

        packages = {
          default = formbio-laava;
          formbio-laava = formbio-laava;
          docker = pkgs.dockerTools.buildLayeredImage {
            name = "ghcr.io/formbio/laava";
            tag = "latest";
            contents = runtimeEnv;
            config = {
              WorkingDir = "/opt/laava";
            };
            created = "now";
          };
        };
      }
    );
}
