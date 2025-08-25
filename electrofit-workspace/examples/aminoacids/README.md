# RESP Charge Fitting with Antechamber and Gaussian: A Tutorial

## Setting Up the Environment

Activate the conda environment and start tleap.

```bash
conda activate AmberTools23
tleap
```

You should see something like:

```
Welcome to LEaP: ...
>
```

## Source the Desired Force Field

```bash
source leaprc.protein.ff14SB
```

## Build Your Sequence

```bash
x = sequence { ACE ALA NME }
```

## Save Files

```bash
savepdb x ala_capped.pdb
saveMol2 x ala_capped.mol2 1
```

The `saveMol2 ... 1` syntax ensures LEaP exports partial charges in the Mol2 file.

You can now open the file and edit it in Avogadro, to make changes.

---

## Optional: Convert MOL2 to a Readable Format for Antechamber

If `antechamber` gives an error, convert the `.mol2` file to `.pdb` and back:

```bash
obabel ala_capped.mol2 -O ala_capped_tmp.pdb
obabel ala_capped_tmp.pdb -O ala_capped_2.mol2
```

This may alter charges. Keep the original `.mol2` file for later reference.

---

## Generate Gaussian Input Using Antechamber

```bash
antechamber -i ala_capped_2.mol2 -fi mol2 -nc 0 -at gaff2 \
            -o gau_input.gcrt -fo gcrt -gv 1 -ge output.gesp
```

## Run Gaussian

```bash
g16 gau_input.gcrt
```

## Generate ESP File

```bash
espgen -i output.gesp -o output_file.esp
```

This file (`output_file.esp`) is needed for the RESP fitting.

---

## Generate RESP Input Files

```bash
antechamber -i gau_input.gcrt.log -fi gout -o output_file.prepi -fo prepi -c resp -s 2
```

This will generate various files, but the most relevant ones are:

- `ANTECHAMBER_RESP1.IN`
- `ANTECHAMBER_RESP2.IN`

---

## Editing `ANTECHAMBER_RESP1.IN`

Modify `ANTECHAMBER_RESP1.IN` to fix ACE and NME charges.

### Original

```
Resp charges for organic molecule

 &cntrl

 nmol = 1,
 ihfree = 1,
 ioutopt = 1,

 &end
    1.0
Resp charges for organic molecule
    0   22
    1    0
    6    0
    1    0
```

### Modified (Fixing ACE and NME)

```
Resp charges for organic molecule

 &cntrl

 nmol = 1,
 ihfree = 1,
 ioutopt = 1,
 iqopt = 2,  # <---- Added line to enforce fixed charges
 qwt = 0.00050,

 &end
    1.0
Resp charges for organic molecule
    0   22
    1   -1
    6   -1
    1   -1
    1   -1
    6   -1
    8   -1
    7    0
    1    0
    6    0
    1    0
    6    0
    1    0
    1    0
    1    0
    6    0
    8    0
    7   -1
    1   -1
    6   -1
    1   -1
    1   -1
    1   -1
```

---

## Preparing `.qin` File for Fixed Charges

Create a `.qin` file (`ala.qin`) using the charges out of the first mol2 file for the residues ACE and NME:

```
 0.112300 -0.366200  0.112300  0.112300  0.597200 -0.567900  0.000000  0.000000 
 0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000 
-0.415700  0.271900 -0.149000  0.097600  0.097600  0.097600
```

These charges are known and widely used (e.g., <https://carlosramosg.com/amber-custom-residue-parameterization>).

- Each line contains 8 numbers.
- Use 6 decimal places.
- Include a space before each line if the first charge is positive.

---

## Running First-Stage RESP

```bash
resp -O -i ANTECHAMBER_RESP1.IN -o resp_output1.OUT \
     -p resp_pch1.pch -t resp_chg1.chg -e output_file.esp -q ala.qin
```

This generates `resp_chg1.chg` with the first-stage fitted charges.

---

## Running Second-Stage RESP (if needed)

If your molecule contains methyl/methylene groups, perform a second RESP fit.

### Before Editing `ANTECHAMBER_RESP2.IN`

```
Resp charges for organic molecule

 &cntrl

 nmol = 1,
 ihfree = 1,
 ioutopt = 1,
 iqopt = 2,
 qwt =  0.00100,

 &end
    1.0
Resp charges for organic molecule
    0   22
    1    0
    6    0
    1    1
    1    1
    6  -99
    8  -99
    7  -99
    1  -99
    6  -99
    1  -99
    6    0
    1    0
    1   12
    1   12
    6  -99
    8  -99
    7  -99
    1  -99
    6    0
    1    0
    1   20
    1   20
```

### After Editing (Fixing ACE/NME Again)

```
Resp charges for organic molecule

 &cntrl

 nmol = 1,
 ihfree = 1,
 ioutopt = 1,
 iqopt = 2,
 qwt =  0.00100,

 &end
    1.0
Resp charges for organic molecule
    0   22
    1   -1
    6   -1
    1   -1
    1   -1
    6   -1
    8   -1
    7  -99
    1  -99
    6  -99
    1  -99
    6    0
    1    0
    1   12
    1   12
    6  -99
    8  -99
    7   -1
    1   -1
    6   -1
    1   -1
    1   -1
    1   -1
```

---

## Run Second-Stage RESP

```bash
resp -O -i ANTECHAMBER_RESP2.IN -o resp_output2.OUT \
     -p resp_pch2.pch -t resp_chg2.chg -e output_file.esp -q resp_chg1.chg
```

---

## Final Output

The file `resp_chg2.chg` contains the final RESP-fitted charges.

**Yey!** You now have RESP-derived charges for your capped amino acid.
