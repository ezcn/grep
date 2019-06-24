# Singularity tutorial
#### Marco Chierici, FBK

### Introduction
[Singularity](https://www.sylabs.io/docs/) is a container solution similar to [Docker](https://www.docker.com/). It was created by the need for reproducible and portable scientific workflows.

A software container allows a user to pack a software application _and all of its dependencies_ into a single package, which can be copied, moved, and shared without further overhead. Containers are also _portable_: they can be simply copied like you would do with regular files. A container is also _agnostic_ about the operating system where it is run, because everything is self-contained.

### Reference material
For a more in-depth guide, please refer to this great NIH tutorial: https://github.com/NIH-HPC/Singularity-Tutorial

### Installation
On KORE, Singularity v2.6.6 is already installed system-wide.

If you want to create your own containers, you need to install Singularity on your machine and `scp` the created containers to KORE. You cannot create containers directly on KORE since Internet access is required to do so.

To install Singularity in Ubuntu or CentOS, follow [these instructions](https://github.com/NIH-HPC/Singularity-Tutorial/tree/master/00-installation).

  - *Hint*: to avoid potential issues, install a version matching KORE's (2.6.6).

### Monolithic vs single-app 

Monolithic containers:
  * they contain everything you need to work on your data (_i.e.,_ it's like having a complete Ubuntu system)
  * you just have one image
  * it is more difficult to maintain them (if you need to update or add an app, you need to rebuild the complete image from scratch)
 
Single-app containers:
  * "one app - one container" phylosophy (_i.e.,_ it's like having the executable binary of the app)
  * you have multiple images, one for each app
  * it is easier to maintain them (you just quilckly update the single image, which is smaller than a monolithic container)

On KORE, we started experimenting with a monolithic container, but then switched to single-app containers for increased efficiency and flexibility.

### Basics: interacting with containers

A Singularity container is like an isolated, stand-alone operating system. So, you can `ssh` into it and interactively execute commands. This may be good for troubleshooting or testing purposes:

```sh
$ singularity shell image.img
Singularity image.img:~> cat /etc/os-release
NAME="Ubuntu"
VERSION="18.04.1 LTS (Bionic Beaver)"
ID=ubuntu
ID_LIKE=debian
PRETTY_NAME="Ubuntu 18.04.1 LTS"
(...)
```

You can also execute a command directly (non-interactively):

```sh
$ singularity exec image.img cat /etc/os-release
NAME="Ubuntu"
VERSION="18.04.1 LTS (Bionic Beaver)"
(...)
```

Finally, it is possible to run a predefined script, which must be set up during the creation of the container - see next section):

```sh
$ singularity run bwa-0.7.17.img
Program: bwa (alignment via Burrows-Wheeler transformation)
Version: 0.7.17-r1188
Contact: Heng Li <lh3@sanger.ac.uk>

Usage:   bwa <command> [options]

Command: index         index sequences in the FASTA format
         mem           BWA-MEM algorithm
(...)
```


### Basics: container creation

To create a container, you need to have a _recipe_ (usually a file simply named `Singularity`) with the instructions for building it: for example, the operating system (ie, Ubuntu 18.04) and which software and dependencies are to be installed.

Singularity is compatible with Docker containers, so if a Docker recipe already exists in a repository, you can use it as the starting point for your Singularity recipe.

A great source of Docker recipes is [Quay.io](https://quay.io/).

Here's an example of a Singularity recipe:

```
Bootstrap: docker
From: quay.io/biocontainers/bwa:0.7.17--ha92aebf_3

%help
    BWA 0.7.17 container

%labels
    Maintainer Marco Chierici

%runscript
    exec bwa "$@"

%post
    mkdir -p /data /projects /work /scratch
```

We'll see the sections shortly. Suppose you save the previous recipe in a `Singularity` file. In order to actually create the container, we issue the following command:

```sh
$ sudo singularity build bwa-0.7.17.img Singularity
```

This will create a container named `bwa-0.7.17.img` from the recipe.

Let's go back to the recipe file. As you see, the file contains a header followed by several sections. The `Boostrap` line tells Singularity to start from a Docker container; the `From` line specifies where to pull the container from. The `%runscript` section tells Singularity which command to run when the container is called with `singularity run`. In this case, we are instructing singularity to call `bwa "$@"`, where `$@` captures all command line options. This minimizes the transition from a "local" call such as

```sh
$ bwa mem -t 16 -o out.sam ref.fa input1_R1.fq input1_R2.fq
```

to a "containerized" call such as

```sh
$ singularity run bwa-0.7.17.img mem -t 16 -o out.sam ref.fa input1_R1.fq input1_R2.fq
```

### Basics: interacting with files

This is probably the nastier part of working with containers. For security reasons, a container cannot access the whole local filesystem. In order to make it interact with local files and directories, you need to "expose" them and make them available through specific "mount points" (directories inside the container).

This is accomplished runtime with the following syntax:

```sh
$ singularity <command> -B <src>:<dest>[:<options>]
```

where `<src>` is the local directory, `<dest>` is the container directory on which `<src>` will be mounted, and `<options>` are optional mount options (ro: read only; rw: read/write).

**Example:** say you need to read/write files from the local directory `/mpbastudies3/181113_Vincenza-Colonna`, making it accessible on the container directory `/data`.

```sh
$ singularity exec -B /mpbastudies3/181113_Vincenza-Colonna:/data bwa-0.7.17.img ls -l /data
total 2481504
-rw-rw-r--    1 chierici mpba          2616 Nov 13  2018 1808KHX-0023_6samples_md5sum.txt
-rw-rw-r--    1 chierici mpba     2538233002 Nov 13  2018 181113_Vincenza-Colonna_1808KHX-0023_6sample_HiSeqX.zip
-rw-rw-r--    1 chierici mpba         10639 Nov 13  2018 181113_Vincenza-Colonna_1808KHX-0023_6sample_HiSeqX_stats.xlsx
drwxrwxr-x    2 chierici mpba           238 Nov 13  2018 AS006-Av-L
(...)
```

The final result is the same you'd get as you executed `ls -l /mpbastudies3/181113_Vincenza-Colonna` on the local system.

### Basics: piping commands

A bioinformatics pipeline is usually the result of many commands, one piped to the following one. A basic example is aligning with `bwa` and converting the alignments to BAM format on the fly:

```sh
bwa mem -t 16 ref.fa input1_R1.fq input1_R2.fq | samtools view -bo out.bam -
```

Now, suppose you have one container for bwa (`bwa-0.7.17.img`) and another container for samtools (`samtools-1.9.img`). Then, the piping is achieved by:

```sh
singularity run bwa-0.7.17.img mem -t 16 ref.fa input1_R1.fq input1_R2.fq | \
    singularity run samtools-1.9.img view -bo out.bam -
```

What if the input/output data are in the directory `/mpbastudies3/181113_Vincenza-Colonna/test` ?

```sh
singularity run -B /mpbastudies3/181113_Vincenza-Colonna/test:/data bwa-0.7.17.img mem -t 16 /data/ref.fa /data/input1_R1.fq /data/input1_R2.fq | \
    singularity run -B /mpbastudies3/181113_Vincenza-Colonna/test:/data samtools-1.9.img view -bo /data/out.bam -
```

This becomes much more readable if we define Shell variables as shortcuts for the singularity calls (mind the capitalization):

```sh
BWA="singularity run -B /mpbastudies3/181113_Vincenza-Colonna/test:/data bwa-0.7.17.img"
SAMTOOLS="singularity run -B /mpbastudies3/181113_Vincenza-Colonna/test:/data samtools-1.9.img"
$BWA mem -t 16 /data/ref.fa /data/input1_R1.fq /data/input1_R2.fq | \
    $SAMTOOLS view -bo /data/out.bam -
```

_(to be continued...)_