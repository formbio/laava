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
        runtimeEnv = with pkgs; [nextflow];
      in {
        formatter = pkgs.alejandra;

        devShells.default = pkgs.mkShell {
          packages = with pkgs;
            [
              git
              lazygit
            ]
            ++ runtimeEnv;
        };
      }
    );
}
