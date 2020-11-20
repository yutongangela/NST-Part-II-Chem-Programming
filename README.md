# NST-Part-II-Chem-Programming

## Exercise-1
This is my first programming exercise. 
In the first exercise, I wrote a programme that determines the Hückel energies and degeneracies of polyene and 3 Platonic solids.

### User Input
Running the programming you will see "Type of the molecule is 1-polyene; 2-sp2 Platonic solid (enter the number)".

Simply follow the prompt, type `1` for polyene and `2` for Platonic solid.
If other characters are used, you will see "Invalid input. Please try again". And it's time to try again:)

If you typed `1` for polyene, you will then see the prompt asking for `cyclic` or `linear`, as well as for the number of atoms. Wrong input, including any thing other than `cyclic`, `linear` or positive integer `n`, will return "This is not a valid polyene molecule". Try again. 

If you typed `2` for Platonic solids, you will then see a prompt asking you to choose from `tetrahedron` or `cube` or `dodecahedraon`. Apologise for not including information of other 2 Platonic solids. The programme will be updated in the future.


### Results 
The programme returns the energies in a list, followed by degeneracies. 


## Exercise-2
In the second exercise, I wrote a script to calculate symmetric stretching and bending frequencies of H2O and H2S. 
Before running the script, you need to:
1. Download the `H2Ooutfiles.tar.bz2` and `H2Soutfiles.tar.bz2` files to your machine
2. Decompress by 
```bash
bunzip2 filename
```
3. Extract from tar archive with 
```bash
tar -xvf filename
```
4. Get path to that folder by typing `pwd`, for example
```console
(base) dhcp-10-241-174-73:H2Ooutfiles yutong$ pwd
/Users/yutong/Downloads/H2Ooutfiles
```
### User Input
Enter whether you want to get information about H2O or H2S. All other inputs are invalid.

Copy and paste path to the folder containing that molecule, for example `/Users/yutong/Downloads/H2Ooutfiles`, press Enter.

### Results
The programme with return information about the molecule plus 3 figures, if your input is correct.

For H2S, the results are:
```bash
H2S has the following information:
----------------Section-1-----------------
Energy at equilibrium geometry: -398.675628032 Hartrees
H-X bond length at equilibrium geometry:  1.35 angstroms
X-H-X bond angle at equilibrium geometry:  94 degrees
------------Section-2-----------------
The fitting of quardratic equation gives k_r = 839 and k_theta = 8.812612070232378e-19 
---------------Section-3---------------
symmetric streching frequency:2668.8362781341502 cm-1
bending frequency: 1281.1750214775745 cm-1
---------------END-------------------
```



If you have any question, please drop me an email (yc439@cam.ac.uk) Thank you!
