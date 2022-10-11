This is the implementation of Discrete-Event Simulation to generate the thermal history of 3D-printed polymer composite parts. The numerical model was published in the following paper.

The citation for the paper is: Bhandari, S., & Lopez-Anido, R. A. (2020). Discrete-event simulation thermal model for extrusion-based additive manufacturing of PLA and ABS. Materials, 13(21), 4985.

The link to the journal paper is: https://doi.org/10.3390/ma13214985

In this implementation, the input files are the gcode file, temperature-dependent material properties in CSV format, and the Input.txt file that defines the model parameters.

The output files are the topology defined by nodes in nodefile.csv, elements (connectivity of nodes) in elementfile.csv, deposition time of element based on toolpath in activation_times.csv, element temperatures in elem_temps.csv, and interpolated node_temperatures in node_temps.csv. An input file for Abaqus FE software is also generated.

The elem_temps.csv file stores the times and temperatures of an element in a given row. The time entries are stored in odd columns (starting at 1) and temperature entries are stored in even columns (starting at 2).

The nodetemps.csv file stores the times and temperatures of a node in a given row. The time entries are stored in odd columns (starting at 1) and temperature entries are stored in even columns (starting at 2).

 Add a custom footer
