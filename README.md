# RDKit tethered minimization and use for rdock
## Introduction
This is a brief example on how to perform tethered minimization on a common 3D substructure of a set of ligands versus a reference ligand. The result of this minimization can be used as is for modeling tasks (ligand alignments in 3D for instance) or like we do in this case here to define output sd files allowing to run tethered docking using rdock. 
Tethered docking has several advantages over free docking, notably the main one: it avoids the "what is the right pose" problem to a certain extent. 
[More information on how we use this in structure based drug design on my blog](https://www.discngine.com/blog/2019/6/6/tethered-minimization-of-small-molecules-with-rdkit-towards-tethered-docking-on-proteins-with-rdock)

## Prerequisits

In order to use this script you need python3 and rdkit. You can install rdkit using conda or use a docker image like the one here: https://hub.docker.com/r/informaticsmatters/rdkit

## Example Usage: 
```
python tetheredMinimization.py sample/reference_pu8_1uyd.sdf sample/ligands_resorcinol.sdf tethered.sdf free.sdf
```
