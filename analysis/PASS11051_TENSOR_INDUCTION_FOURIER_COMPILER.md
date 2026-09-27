# Pass 11051 — the 27D compiler is nine parallel qutrit Fourier blocks

Producer: `analysis/w33_pass11051_tensor_induction_fourier_compiler.py`
Certificate: `data/w33_pass11051_tensor_induction_fourier_compiler.json`

The key identity is now constructive:

`V tensor Ind_L(theta)`
` ~= Ind_L(Res_L(V) tensor theta)`
` ~= Ind_L(Reg L)`
` ~= Reg(H27)`.

A cyclic vector in the induced-times-Schrodinger carrier generates an orthogonal regular orbit. The resulting 27x27 synthesis matrix has exactly 81 nonzero `mu3` entries.

Its support graph splits into **nine disjoint 3x3 components**. After dephasing every component is exactly the qutrit Fourier exponent matrix

`[[0,0,0],[0,1,2],[0,2,1]]`.

So the normalized compiler is nine parallel `F3` mixers plus permutations and phase gauges, not a generic 27x27 unitary.
