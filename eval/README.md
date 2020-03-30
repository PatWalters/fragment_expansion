```
./gen_restricted_confs.py --smi data/x0161.smi --fix data/x0161_lig.sdf --out data/x0161_confs.sdf

./score_confs.py --prot data/x0161_prot.pdb --lig data/x0161_confs.sdf --ref data/x0161_lig.sdf --out data/x0161_out.sdf
```
