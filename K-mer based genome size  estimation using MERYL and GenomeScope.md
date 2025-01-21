Many genome assemblers require a genome size estimation. This can be achieved in a number of different ways, the easiest of which is to simply look it up [(e.g, for animals)](https://www.genomesize.com/) for your organism or a different, but closely related species. An alternative to this is to do it using K-mers.

To do this, I use [meryl](https://github.com/marbl/meryl) to count K-mers and a genome profiling tool called [GenomeScope](https://github.com/tbenavi1/genomescope2.0)  to get useful estimates like heterozygosity and estimated genome size. 

I will assume that the programs are both installed, but if you're unsure how to do this, you can follow the instructions on each tool's github repository. In a nutshell, meryl can be installed using a `make` file, and my installation of GenomeScope is from `conda` .

# Create a K-mer database and a histogram from it
Meryl consists of several commands. The one we will be using is `count`. There are several other required inputs:
1. `k` is the k-mer of length *k* . Here, we use 21.
2. `output` is the output directory name
3. sequence data, here we simply call them positionally `$r1` and `$r2`

The `threads` and `memory` arguments are optional, but if you're working on an HPC, you're likely to get an angry email from your network administrator (automated or otherwise) if you don't watch your resource usage. Here, I will limit the `meryl` to a modest 8 CPUs and 32 Gb of RAM. How long these steps take will depend on the size of your genome and the resources you allocated.

After creating a K-mer database, we will use the `histogram` tool to create a histogram of the K-mers, which we output here with the extension `.hist`. 

```
#set path your meryl installation
PATH=PATH/TO/MERYL/INSTALLATION:$PATH

#set variables names
r1=PATH/TO/READS.fastq.gz
r2=PATH/TO/READS2.fastq.gz
out=some_informative_but_short_name

#count K-mers using meryl
meryl count k=21 $r1 $r2 output ${out}.meryl threads=8 memory=32


#use the two k-mer databases we created and create a histogram
meryl histogram ${out}.meryl > ${out}.hist threads=25 memory=75
```

Now that we have a K-mer histogram, we will pass the `.hist` file to GenomeScope, which is an R script that the authors of the program wrote to model genome characteristics such as size, heterozygosity, and error.

The GenomeScope R script requires several input arguments:
1. `-i` is the input `.hist` file that we created in the previous step
2. `-k` is the K-mer length. Since our database that we created is based on `K=21`, here we set it to 21
3. `-o` is the output directory name
4. `-n` is the output file prefix
5. `-p`  (optional)  is the ploidy of the organism
6. `fitted_hist` outputs a fitted histogram

```
gs=PATH/TO/GenomeScope2/genomescope.R

Rscript $gs -i ${out}_gt1.hist -k 21 -o ${out}_gt1_gs -n ${out}_gt1 -p 2
```

This should be done in a matter of seconds. In the output directory, you should find several `.png` images of the figures that GenomeScope created based on the non-linear model fit to your data, like the one below.
![[images/Bf05_gt1_transformed_linear_plot.png]]

The key things to take away here are:
1. `len`, which is the estimated genome size in bp
2. `aa` and `ab`, which are homozygosity and heterozygosity, respectively
3. `kcov` is the haploid peak
