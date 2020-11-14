using Pkg
metadata_packages = [
    "Oscar",
    "LinearAlgebra",
    "Singular",
    "GroebnerBasis",
    "MacroTools",
    "OrderedCollections"
]

for package = metadata_packages
    Pkg.add(package)
end

Pkg.resolve()