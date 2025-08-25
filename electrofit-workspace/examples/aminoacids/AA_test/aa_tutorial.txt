    conda activate AmberTools23
    tleap

You should see something like:
    Welcome to LEaP: ...
    > 

Source the Desired Force Field:
    source leaprc.protein.ff14SB

Build your sequence:
    x = sequence { ACE ALA NME }

Save:
    savepdb x ala_capped.pdb
    saveMol2 x ala_capped.mol2 1

The saveMol2 ... 1 syntax ensures LEaP exports partial charges in the Mol2 file.

Now you can open the file and edit it in Avogadro, to make changes.
Depending on weather you can use command (1) (look below) rigth away or antechmaber gives you an error do the following: 
    
    Optional:

    Next we need to convert the mol2 file into a file that can be read be antechamber.
    To do so, we convert the file into pdb and than back into a mol2 file (ala_capped_2). 
    Note, this will change the charges. But we do need the initial charges of the ACE and NME cap later.
    So dont delete the initial mol2 file or look up the charges from e.g., https://carlosramosg.com/amber-custom-residue-parameterization.

Now we execute antechamber to generate gaussian input (1):
    antechamber -i ala_capped_2.mol2 -fi mol2 -nc 0 -at gaff2 -o gau_input.gcrt -fo gcrt -gv 1 -ge output.gesp 

Then we run gaussian: 
    rung16 gau_input.gcrt

Then we generate an ESP from the gesp: 
    espgen -i output.gesp -o output_file.esp

We need this "output_file.esp" later...

Then we generate RESP input. I do this with: 
    antechamber -i gau_input.gcrt.log -fi gout -o output_file.prepi -fo prepi -c resp -s 2
This will give you along with the desired input files charges (QOUT, qout, ...) and some other files. 
Ignore them their are wrong. Interesting is only ANTECHAMBER_RESP1.IN, ANTECHAMBER_RESP2.IN for now.

Open the file ANTECHAMBER_RESP1, and edit it.

Belwo "Resp charges for organic molecule" we need to replace the numbers in the second column for the first 6 (excluding the firs row) atoms (ACE) and the last 6 atoms (NME) with "-1". 
The "-1" tells antechamber to hold the charges for these atoms constant during the fit. 
Also add one "iqopt = 2," below  "ioutopt = 1," this will say antechamber to use the charges you give as input.
Meaning, you can feed the correct charges for ACE and NME, which you do not want to be changed during the fit. 

In the end it can look like this:
"""
    Resp charges for organic molecule

    &cntrl

    nmol = 1,
    ihfree = 1,
    ioutopt = 1,
    iqopt = 2, <---- new line added 
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
"""

This is done for alanin, but would look similar for your molecule.

Now we only need to prepare a .qin file that will tell antechamber the charges to hold constant. 
The file looks like this:
"""
 0.112300 -0.366200  0.112300  0.112300  0.597200 -0.567900  0.000000  0.000000 
 0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000 
-0.415700  0.271900 -0.149000  0.097600  0.097600  0.097600
"""
You have always 8 numbers with 6 digits per line. One space in front of each line if the first charge in that line is positive, else no space.
Between charges two spaces if no "-" sign and one space if a "-" sign.
I called the file ala.qin here. 

After this you can do the first stage of RESP fitting:

    resp -O -i ANTECHAMBER_RESP1.IN -o resp_output1.OUT -p resp_pch1.pch -t resp_chg1.chg -e output_file.esp 

Wich will give you your first set of charges. If you dont have methyl or methylene groups you can stop here (groups within ACE and NME does not count, since we constrain them to be constant).
If you have some, you want to do second stage RESP.

Here agian you need to set the frist and last six numbers in the second column under "Resp charges for organic molecule".
Most likely  "iqopt = 2," is already set to two, because we will feed the output charges from the first stage RESP fit (resp_chg1.chg) as input.
Before:
"""
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

"""
After:
"""
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
"""

Second stage: 
    resp -O -i ANTECHAMBER_RESP2.IN -o resp_output2.OUT -p resp_pch2.pch -t resp_chg2.chg -e output_file.esp -q resp_chg1.chg

Yey, you made it. The file resp_chg2.chg contains your final charges. 