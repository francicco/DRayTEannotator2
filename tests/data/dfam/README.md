This small compatibility dataset contains the 80 nt consensus emitted by HMMER
3.4 `hmmemit -c tutorial/MADE1.hmm` (Dfam accession DF0000629.2) and a short
negative control. The source is the official HMMER 3.4 distribution:
https://eddylab.org/software/hmmer/hmmer-3.4.tar.gz

The model is historical; this fixture does not establish biological accuracy
against the complete Dfam 4.0 release. The optional external-tool test accepts a
FamDB database containing this model with appropriate classification metadata.
For a synthetic compatibility database, any added GA cutoff must be identified
as a test cutoff rather than a measured Dfam threshold.
