"""PDC (Proteomic Data Commons) metadata access for QPX.

Self-contained GraphQL helpers used by the CDAP sample/run builder to recover
the channel -> aliquot -> biological-sample mapping that CDAP ``.psm`` files do
not carry. Network access only; no pridepy dependency (pridepy remains the
download layer for ``pdc2qpx``).
"""
