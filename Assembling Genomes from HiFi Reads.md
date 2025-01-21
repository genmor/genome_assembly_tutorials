The main advantage of PacBio's HiFi sequencing platform are its read lengths and accuracy. If the marketing material is anything to go by, you can expect to get reads up to 25 kbp long and 99.9% accurate base calls. For eukaryotic genomes, these qualities are great because these read lengths are long enough to span most (if not all) repetitive elements and the accuracy more or less eliminates the need for assembly polishing.

Numerous programs have been written, or re-written to handle these read types. All of the programs I've used so far have been fairly easy to install, use, and troubleshoot. All of the programs appear to have an active community and developers surrounding it, so support should fine for a little while.

The programs I will be showing are shown below in a table. Note that this is strictly for genome assembly—a some programs also do haplotig/duplicate contig purging and phasing (with appropriate data).

| Assembler                                          | Does it Purge? | Does it Phase? |
| -------------------------------------------------- | -------------- | -------------- |
| [canu](https://github.com/marbl/canu)              | No             | No             |
| [hifiasm](https://github.com/chhylp123/hifiasm)    | Yes            | Yes            |
| [Flye](https://github.com/mikolmogorov/Flye)       | No             | No             |
| [IPA](https://github.com/PacificBiosciences/pbipa) | Yes            | Yes            |


Installation of these programs are generally straightforward, so I will leave it to the reader to go through the documentation of the program they choose for installation. Furthermore, running these programs are quite easy—the only concern that I would note is that these are generally resource-hungry programs (i.e., high memory and CPU needs). Thus, it's important to keep in mind how much you is available to you (whether on a high performance cluster or on a beefy personal machine). 

## Canu

```
spp=whatever_you_want

canu \
 -p ${spp}_canu -d ${spp}_canu_out \
 genomeSize=1.3g \
 useGrid=false \
 maxMemory=180 maxThreads=40 \
 -pacbio-hifi /path/to/where/your/PacBio_HiFi_reads/live/*.fastq.gz
```

The main arguments to be aware of are `genomeSize`, `useGrid`, `maxMemory`, `maxThreads`, and the `-pacbio-hifi` flag. `genomeSize` is the estimated or expected haploid genome size of the your organism. You might get that information from a database, or a flow cytometry experiment, or using a K-mer based estimation. `useGrid`, `maxMemory`, and `maxThreads` are resource-related arguments. How many CPU threads and how much RAM are available to you? Do you have a job scheduler (like SLURM)? The high performance computing system at the University of Calgary uses SLURM to allocate resources, so I could set `useGrid=true`, but I find it simpler to request resources directly and limit the the amount of RAM and CPU threads (and also a good way to avoid getting angry emails from system administrators). The `-pacbio-hifi` flag tells Canu that the reads are HiFi reads (it's also able to handle PacBio CLR and Nanopore reads). The last argument is the file path to the fastq files. Here, I use the wildcard `*` with the extension `fastq.gz` because it's easier to reference them this way. There are two additional flags `-p` and `-d`—these are the file prefix and the output directories respectively. 

The main outputs of Canu are:
1. prefix.report
2. prefix.contigs.fasta
The prefix is derived from the `-p` argument. The report file contains a wealth of summary information that should be reported in manuscript or a write-up of some sort. The `contigs.fasta` file is the main assembly. Canu also outputs an `prefix.unassembled.fasta` which contains reads that couldn't be incorporated into the main assembly.

## hifiasm

```
spp=whatever_you_want

hifiasm /path/to/where/your/PacBio_HiFi_reads/live/*.fastq.gz \
 -o ${spp}_hifiasm \
 --hg-size 1.3g \
 --primary \
 -t 40 
```

hifiasm tries to be an all-in-one assembler: it assembles HiFi reads to contigs, purges duplicate contigs and haplotigs (by default), and if provided with Hi-C data, can also output a fully phased assembly. I haven't tried the latter mode, so I won't cover it. The main arguments that I use here are `--hg-size` and `--primary`. `--hg-size` is the estimated/expected size of the haploid genome, and `--primary` outputs a primary and an alternate genome. `-t` is the number of CPUs the program will use, and `-o` is the output prefix. 

There are additional flags that you can set (that I didn't) that might help with the output. I haven't messed with these too much because I find the defaults to be good enough. These are `--hom-cov` and `-l`, which respectively tweak the homozygous coverage peak and the level of purging. `--hom-cov` is by default inferred automatically, but if you find that outputted assembly is significantly larger than you expected, you might set it manually. `-l` changes how much hifiasm will purge:
1. `-l 0` will disable purging
2. `-l 1` will purge "only contained haplotigs" (`-s 0.75`)
3. `-l 2` will purge "all types of haplotigs" (`-s 0.75`)
4. `-l 3` will purge "all types of haplotigs in the most aggressive way" (`-s 0.55`)
Purging defaults to `-l 3`. Purging is done based on sequence similarity, which can be controlled using `-s` and specifying a value 0–1. The default similarity threshold are given in parentheses. I would advise against going below `-s 0.55`. 

hifiasm outputs assembly graphs in `.gfa` format, rather than `.fasta`. The main outputs in the mode specified above are:
1. prefix.p_ctg.gfa
2. prefix.a_ctg.gfa
Here, `p_ctg.gfa` and `a_ctg.gfa` the primary and alternate assemblies. Our main focus is on will be on the primary assembly, as the alternate assembly contains haplotigs and duplicate contigs.

To convert a `.gfa` file to a `.fasta`, you can do the following:

```
awk '/^S/{print ">"$2"\n"$3}' p_ctg.gfa | fold > p_ctg.fasta
```

You can specify each `.gfa` individually and convert them one at a time. However, it's likely you'll want to convert all of the `.gfa` files into a `.fasta`, so I just convert all them in a `for` loop like this:

```
for i in *.gfa 
do
awk '/^S/{print ">"$2"\n"$3}' $i | fold > ${i%.gfa}.fasta
done
```


## Flye

```
spp=whatever_you_want

flye --pacbio-hifi /path/to/where/your/PacBio_HiFi_reads/live/*.fastq.gz \
--genome-size 1.3g \
-o ${spp}_flye \
-t 40 \
```

Just like with Canu or hifiasm, the main argument to pay attention here are `--genome-size` (the estimated genome size) and `-t` (the number of CPU threads that will be used). `-o` is the output directory. Flye offers a number additional functionalities because it can handle several other data types (Nanopore and PacBio CLR), but I won't go into them here.

The outputs from Flye are pretty straightforward:
1. `assembly.fasta`
2. `assembly_graph.gfa`
3. `assembly_info.txt`
In the event that the Flye fails for some reason, you can try to resume the using the `--resume` and `--resume-from` options.

## IPA

```
ipa local -i /path/to/fofn \
--genome-size 1300000000 --nthreads 10 --njobs 4
```

IPA, or the "Improved Phased Assembler" is the "in-house" assembler from PacBio. It, like hifiasm, tries to be an all-in-one assembler and performs read assembly, contig/haplotig purging, and phasing. IPA is designed to work with a job scheduler (like SLURM) or locally. Again, I find it easier to control the number of requested resources by running the program locally, so we specify `local`. Along similar vein, I will specify `--nthreads` and `--njobs` which splits the tasks that IPA does across 4 jobs, each using 10 CPU threads. Here, I used an `fofn` (a file of file names) as the input argument for where the reads are located instead of supplying a path. The contents of this file should look like the following:
```
/path/to/pacbio_hifi_reads/read1.fastq.gz
/path/to/pacbio_hifi_reads/read2.fastq.gz
```
It's simply a list of files paths, each on their own line. The path should be the full path to the `fastq.gz` file you have for the specimen. 

As mentioned before, IPA performs duplicate contig and haplotig purging by default. If you want to turn this off, you can do so by setting the flag `--no-purge-dups`.

The main outputs of IPA are:
1. `p_ctg.fasta`
2. `a_ctg.fasta`
Similar to the other assemblers, `p_ctg.fasta` and `a_ctg.fasta` are primary and alternate assemblies.

One thing to note about the `fasta` output from IPA is that the header might contain `\`s, which will cause other bioinformatics programs to throw an error. You can avoid this by removing these and replacing them with some other, more innocuous character like `_` using the following command:
```
sed -i 's!/!_!g' p_ctg.fasta
```


