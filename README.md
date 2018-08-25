# Generation of secure elliptic curves

This was my TIPE for the year 2017/2018. By publishing it I hoped to make high level implementations of common algorithms related to elliptic curves more easily available, including the CM method and the Pollard-rho attack on discrete logarithms. If it can be an example to anyone that'd be great too.

# What's in it

The productions are all there: slides, MCOT, DOT, abstract...  There are also many useless documents I produced during the year, but the important content is certainly the code in the `algo` folder (and the `attacks` folder, though it is less relevant to this TIPE).  Nearly all functions in there are documented, and example usage can be found in `tests.sage`.

The distinction between the `articles` and `resources` directories is blurry, but surely there is one.
 
# Usage

To generate a secure elliptic curve of a given security `n` with cofactor smaller than `k`, load the files `rd_approach.sage` and/or `cm_approach.sage` in the sage interpreter and call either `gen_random` or `gen_cm` with parameters `2^{2n}` and `k`. See the documentation of each for details.

# Reproducibility of results

To generate graphs, you first need to generate performance files using the `bench_*` functions in `tests.sage`. Modify them to suit your needs, but keep in mind that generating curves of large order can take a _very_ long time, as the implementations given here are not optimized.

Then you can call the `plot_*` functions to get the actual graphs, and don't forget to modify the `bits` array according to the computed benchmarks.

# Memento

I thought after a while that this subject was far too deep for me to understand and be productive in. And I was certainly right for the maths part, of which I understood only the tiny portion shown in my slides. For this reason I focused on programming, something I could do. It paid off, especially once I realised I could actually implement the CM method, allowing me to compare two approaches. This was at the very end of the school year, when I was really questionning the interest of my work. Fake it 'till you make it.
-- C. Nolan