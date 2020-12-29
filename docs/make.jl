using Documenter
using CosmoCentral

makedocs(
    sitename = "CosmoCentral",
    format = Documenter.HTML(),
    modules = [CosmoCentral]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
