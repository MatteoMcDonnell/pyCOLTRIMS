# pyCOLTRIMS

pyCOLTRIMS is a very much in-progress python framework(ish) for doing data processing from data taken on a COLTRIMS. 
This is amazingly specialized for the way we can output data in a .bin file from our sever.

Currently, the iPython notebook is only written to handle double coincidence, but it is easily expandable to triple coincidence.  


The only difference between the two convert files is that the pump-probe file adds another column for stage voltage. 
This means you can use this for any experiment with a stage, not just pump-probe.


## Razib:

The convert files are what you are looking for. The difference between convert.py and convert-pumpprobe.py is that pumpprobe is ... for pump probe. It outputs an extra column. 
The GUI file is for chosing the input bin file and specifying the channels. 

the files that start with "initialProcessing" are the iPython Notebooks I use for analysis, everything else is just misc. 



