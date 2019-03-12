# ModelingClimateTrends

R code that accompanies the article

P. F. Craigmile and P. Guttorp,

Modeling and assessing climatic trends

Revision submitted for review

Contact: pfc@stat.osu.edu


## Correction

This code assumes that the uncertainty column in the Berkeley Land and
Ocean summary provides the standard deviation of the global mean temperatures.

According to
<http://berkeleyearth.lbl.gov/auto/Global/Land_and_Ocean_summary.txt>
"Uncertainties represent the 95% confidence % interval for statistical
and spatial undersampling effects as well as ocean biases."

We recommended the uncertainty values by 1.96 to obtain the standard
error for the global mean temperatures.



## Acknowledgments

Craigmile is supported in part by the US National Science Foundation under grants DMS1407604 and SES-1424481, and the National Cancer Institute of the National Institutes of Health under Award Number R21CA212308. The content is solely the responsibility of the authors and does not necessarily represent the oﬃcial views of the National Institutes of Health. Guttorp is supported in part by the US National Science Foundation grants DMS-1106862 (STATMOS) and DMS-1401793 and by NordForsk grant 74456 (eSACP).
