#!/bin/bash

PACKAGES_FILE="apt-packages.txt"
ACTION=${1:-"check"}
BASE_IMAGE="continuumio/miniconda3:24.11.1-0"

# Run apt-cache command inside Docker container using base image
run_apt_cache() {
    local cmd="$1"
    docker run --rm --platform=linux/amd64 $BASE_IMAGE bash -c "apt-get update > /dev/null 2>&1 && $cmd"
}

check_versions() {
    echo "Using base image: $BASE_IMAGE"
    
    echo "Current pinned versions:"
    cat $PACKAGES_FILE
    echo -e "\nAvailable versions:"
    while IFS= read -r package; do
        pkg_name=$(echo $package | cut -d'=' -f1)
        echo "=== $pkg_name ==="
        run_apt_cache "apt-cache madison $pkg_name" | head -3
    done < $PACKAGES_FILE
}

update_package() {
    local package=$1
    local version=$2
    
    if [ -z "$version" ]; then
        echo "Getting latest version for $package..."
        # Get latest version
        version=$(run_apt_cache "apt-cache madison $package" | head -1 | awk '{print $3}')
        if [ -z "$version" ]; then
            echo "Error: Could not find version for $package"
            return 1
        fi
    fi
    
    # Update in packages file
    if [[ "$OSTYPE" == "darwin"* ]]; then
        # macOS sed syntax
        sed -i '' "s/^$package=.*/$package=$version/" $PACKAGES_FILE
    else
        # Linux sed syntax
        sed -i "s/^$package=.*/$package=$version/" $PACKAGES_FILE
    fi
    echo "Updated $package to $version"
}

case $ACTION in
    "check")
        check_versions
        ;;
    "update")
        if [ -z "$2" ]; then
            echo "Usage: $0 update <package> [version]"
            exit 1
        fi
        update_package $2 $3
        ;;
    "update-all")
        echo "Updating all packages to latest versions..."
        while IFS= read -r package; do
            pkg_name=$(echo $package | cut -d'=' -f1)
            update_package $pkg_name
        done < $PACKAGES_FILE
        ;;
    *)
        echo "Usage: $0 {check|update <package> [version]|update-all}"
        exit 1
        ;;
esac
