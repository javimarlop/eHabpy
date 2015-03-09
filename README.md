eHabitat+
============

[![Gitter](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/javimarlop/eHabpy?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

<a rel="license" href="http://creativecommons.org/licenses/by-sa/3.0/deed.en_US"><img alt="Creative Commons License" style="border-width:0" src="http://i.creativecommons.org/l/by-sa/3.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-sa/3.0/deed.en_US">Creative Commons Attribution-ShareAlike 3.0 Unported License</a>.

[![DOI](https://zenodo.org/badge/4755/javimarlop/eHabpy.png)](http://dx.doi.org/10.5281/zenodo.10612)


**eHabitat+** GRASS GIS 7 scripts and Python 2.7 library for automatic delineation of habitats within protected areas (PA) and calculation of maps of probabilities to find areas presenting similar ecological characteristics to those found in PA within the corresponding ecoregion. A habitat similarity index (HSI) is computed based on the ratio between the extent of similar areas around the PA and the PA extent, as well as some  landscape metrics and indices to characterize similar areas to PA. Processed results are being updated and can be accessed through the [DOPA Explorer](http://ehabitat-wps.jrc.ec.europa.eu/dopa_explorer/).

## Notes

* gdal command line utilities must be installed and accessible in order to run some landscape metrics.

* The Python library can be executed using [Tzar](https://tzar-framework.atlassian.net/wiki/display/TD/Tzar+documentation) distributed computing framework thanks to [doctorluz](https://github.com/doctorluz).

