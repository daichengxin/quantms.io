Introduction to quantms.io
======================================

The ``.qms`` folder will contain multiple files that will be
used to describe the project, the samples, the data acquisition and the
data processing, in addition to the data files. Each of these files will be described in the following
sections:

-  :doc:`metadata`: A json file for metadata about the
   analyzed project
-  :doc:`ae` or :doc:`de`: A csv file based on the `MSstats <https://github.com/Vitek-Lab/MSstats>`_ format for either absolute expression or
   differential expression.

Some general rules for all the files:

-  The files are tab-delimited, json or parquet files
-  Parquet files are compressed and can be read with pandas.

Common data structures and formats
----------------------------------

We have some concepts that are common for some outputs and would be good
to define and explain them here:

Peptidoform
~~~~~~~~~~~

-  **Peptidoform**: A peptidoform is a peptide sequence with
   modifications. For example, the peptide sequence ``PEPTIDM`` with a
   modification of ``Oxidation`` would be ``PEPTIDM[Oxidation]``. The
   peptidoform show be written using the `Proforma
   specification <https://github.com/HUPO-PSI/ProForma>`__. This concept
   is used in the following outputs:

   -  :doc:`psm`
   -  :doc:`feature`
   -  :doc:`peptide`

Modifications
~~~~~~~~~~~~~

-  **Modifications**: A modification is a chemical change in the peptide
   sequence. Modifications can be annotated as part of the Proforma
   notation inside the peptide or as a separate column. When annotating
   the modification as a separate column, the format should be as close
   as possible to the `mzTab format
   notation <https://github.com/HUPO-PSI/mzTab/tree/master/specification_document-releases/1_0-Proteomics-Release>`__.
   The modifications will encode the following information on each
   peptide or psm:

   -  Modification name or accession: For example, ``Oxidation`` or
      ``UNIMOD:35``. Modifications SHOULD be reported using UNIMOD. If a
      modification is not defined in UNIMOD, a CHEMMOD definition must
      be used like ``CHEMMOD:-18.0913``, where the number is the mass
      shift in Daltons.
   -  Position: The position of the modification in the peptide
      sequence. Terminal modifications in proteins and peptides MUST be
      reported with the position set to 0 (N-terminal) or the amino acid
      length +1 (C-terminal) respectively. For example, ``1`` or
      ``1,2,3``.
   -  Localization Probability: The probability of the modification
      being in the reported position.

Those three properties can be combined in one string as:

::

   {position}({Probabilistic Score:0.9})|{position2}|..-{modification accession or name}

For example:

::

   1(Probabilistic Score:0.8)|2(Probabilistic Score:0.9)|3-UNIMOD:35`. 

This concept is used in the following outputs:

-  :doc:`psm`
-  :doc:`feature`
-  :doc:`peptide`
