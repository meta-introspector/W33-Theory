# Pass 11170 — perfect qutrit gates exist for 1, 2, 3 and 5 qutrits, never for 4; an explicit perfect five-qutrit gate

Producer: `analysis/w33_pass11170_perfect_gate_ladder.py`
Regression: `tests/test_w33_pass11170_perfect_gate_ladder.py`

**Four qutrits: impossible (known).** A perfect four-qutrit gate of any kind would have an AME(8,3) Choi state. The
shadow inequalities exclude AME(8,3): Huber, Eltschka, Siewert and Gühne, arXiv:1708.06298. The Huber–Wyderka table of
AME states lists AME(8,3) as "No (Shadow)"; we read its data file. So a register of four qutrits can never scramble
perfectly in one tick. This is a corollary of published work, not a new result.

**Independent check on the F₉-linear slice.** F₉-linear perfect gates on m qutrits are the unitary m×m matrices over F₉
with every entry and every 2×2 minor nonzero (Pass 11164). We ran an exhaustive clique search:

| m | rows (unit vectors, all entries nonzero) | perfect row sets |
|---|---|---|
| 2 | 16 | 32 (= 64/2!) |
| 3 | 192 | 2048 (= 12 288/3!) |
| 4 | 1280 | **0** |

The m = 4 zero means there is no Hermitian self-dual [8,4,5]₉ code (Pass 10970's [8,4,5]₇ is a different code, over F₇): such a code has a generator [I | ηK] with K unitary
and superregular.

**Five qutrits: an explicit perfect gate.**
* Among circulant matrices over F₉ there are 24 perfect ones for m = 3 (including Pass 11164's P² + i𝟙), 0 for m = 4,
  and **80 for m = 5**.
* The first, with first row (i, i, 1−i, −1, 1−i), embeds as S ∈ Sp(10,3). Its Choi state has entropy 5 trits on **all
  252** 5|5 cuts, so it is AME(10,3) and the gate is perfect.
* AME(10,3) is known: [[10,0,6]]₃ comes from Grassl–Gulliver's circulant Hermitian self-dual [10,5,6]₉ code, via
  Grassl–Rötteler arXiv:1502.05267 and Danielsen arXiv:1106.2428. Our circulants are presumably instances of it.

**The ladder.** Perfect qutrit ticks exist for registers of 1, 2, 3 and 5 qutrits, and not for 4.

**Correction to Pass 11164 (over-read).** Pass 11164 called the closed form K = P² + i𝟙 new. It is not. With
η = 1 + i (norm −1), [I | ηK] generates a Hermitian self-dual [6,3,4]₉ double-circulant code. That is Grassl–Gulliver's
construction (Des. Codes Cryptogr. 52 (2009) 57–81), the code behind [[6,0,4]]₃ in Grassl–Rötteler. Only its reading as
a perfect three-qutrit tick in W(3,3)'s Clifford group is ours. The 11164 note, the paper and the ledger are corrected
accordingly.
