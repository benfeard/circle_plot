# Circle Plots for Cross-Correlations

This is an update version of the code used in DOI:https://doi.org/10.1016/j.bpj.2017.08.026 

Please read the research article for the full methods behind this approach.

## CrossCirclePlot.py

#usage: python CircleDMS.py apoe4.ct e4_average_filtered.dat example_output.ps

This code helps you visualize the cross correlations (both positive and negative) between residues in your protein. Input data should come from molecular dynamics simulations or experimental data.

It is recommended to minimize the number of correlations included to produce a cleaner and easier to understand plot.

## Process_raw_crosscorr.py

If you have used [Wordom](http://wordom.sourceforge.net) to generate the cross-correlation data for your protein, this script can help format it for the Circle Plot.
