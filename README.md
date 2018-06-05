# IUPipe
---
A Python pipeline for IUPred 1.0
---
This tool is a pipeline for batch prediction of intrinsic disorder on collections of protein sequences in FASTA format.
The intrinsic disorder prediction is done by IUPred 1.0 <sup>1</sup> (Dosztanyi *et al.*, 2005), a simple yet robust 
method often used as a basic benchmark. The main advantages of IUPred 1.0 are:

- It is computationally cheap.
- It is not a binary score.

---
<sup>1</sup> As of May 2018, a new Python3-based version of IUPred, named IUPred2A, has been released. 
An adapted version of this tool will be created at some point.

## Bibliography

Dosztanyi, Z., V. Csizmok, P. Tompa, and I. Simon. 2005. “IUPred: Web Server for the Prediction of Intrinsically Unstructured Regions of Proteins Based on Estimated Energy Content.” Bioinformatics 21 (16): 3433–34. https://doi.org/10.1093/bioinformatics/bti541.

Mészáros, Bálint, Gábor Erdős, and Zsuzsanna Dosztányi. 2018. “IUPred2A: Context-Dependent Prediction of Protein Disorder as a Function of Redox State and Protein Binding.” Nucleic Acids Research, no. June (June): 1–9. https://doi.org/10.1093/nar/gky384.
