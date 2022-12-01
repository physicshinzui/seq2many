# seq2many

## Given a reference sequence, seq2many generates mutant sequences specified.
seq2many is possible to carry out
- single-point mutation (you specify)
- multiple mutations (Give an input file for the list of mutations)
- deep mutational scanning (20 mutations per one position)

Requirements: 
- Python > 3.8.5 


## Example
Suppose the following reference sequence (Let this file name be `ref.seq`):
```
>6M17_E
CPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGD
EVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQA
GSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELL
```

### Single-point mutation (The amino acid at position 0 mutates into K)
`python seq2many.py -i ref.seq -m single -p 0 -ia K`

### Multiple mutation 
`python seq2many.py -i ref.seq -m multi -ml mutlist.inp`

where
```mutlist.inp
A:1   # The amino acid at the position 1 muates into A.
L:2   # ... at the position 2 mutates into L.
K:5   # ... at the position 5 mutates into K.
```

### Deep mutational scanning
`python seq2many.py -i ref.seq -m deep`

This command generates many sequence files.
