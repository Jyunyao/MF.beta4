\name{forest_biodiversity_data}
\alias{forest_biodiversity_data}
\title{Biodiversity data for six European forests}

\description{
In addition to row and column names, the original dataset consists of three columns: the “plot” column indicates the name of plots, the “full_species_original” column includes species names, and the column “basal area” (as a proxy of species abundance). See Scherer-Lorenzen et al. (2023) for the original dataset. Because missing values of “basal area” in the original dataset were imputed by the mean of the same species within the country, and some basal areas were added for almost indistinguishable species, the dataset provided with the package is slightly different from the original dataset. Also ,for analysis, original three column names of the dataset turn into “plotID” (the name of plots), “species” (species names) and “abundance” (species abundance). 
}

\usage{
data("forest_biodiversity_data")
}

\format{
a data.frame with 481 rows (the total number of combinations of plot by tree species in 209 plots) and 3 columns (plotID, species name, and the corresponding basal area as species abundance).
}

\references{

Scherer-Lorenzen, M. et al. (2023). The functional significance of tree species diversity in European forests - the FunDivEUROPE dataset [Dataset]. Dryad. https://doi.org/10.5061/dryad.9ghx3ffpz

}


\keyword{datasets}
