---
layout: post
title: "Introduction to RNA-Seq - updated"
modified:
categories: courses
excerpt:
tags: [ngs, rna, introduction, galaxy, qc, alignment, de]
image:
  feature:
date: 2015-10-22
---

<section id="table-of-contents" class="toc">
  <header>
    <h3>Overview</h3>
  </header>
<div id="drawer" markdown="1">
*  Auto generated table of contents
{:toc}
</div>
</section><!-- /#table-of-contents -->

We will be retracing most of the steps required to get from an Illumina FASTQ sequence file received from a sequencing facility as part of an RNA sequencing analysis all the way to transcript assembly and identifying  differentially expressed genes in two different ways. We will focus on the following  steps: initial quality control, read mapping, transcript isoform determination (reference-guided) followed by quantification. This session borrows material from [Jeremy's RNA-seq sample history](https://main.g2.bx.psu.edu/u/jeremy/p/galaxy-rna-seq-analysis-exercise), a [Galaxy workflow](https://usegalaxy.org/u/mejia-guerra/w/basic-rna-seq-analysis---differential-expression-functional-genomics-workshop-2012) and the [UC Davis Bioinformatics Core workshop](http://training.bioinformatics.ucdavis.edu/docs/2012/05/RNA/index.html). 


#### Access to Galaxy
{:.no_toc}

Before tackling this module you should have a [working knowledge of Galaxy](http://bioinformatics.sph.harvard.edu/ngs-workshops/courses/introduction-to-galaxy/), understand how to find relevant tools, annotate your history or get new data sets from the data library.  As a result we will minimize the screenshots and pointers in this document, focusing just on _what_ to do while guiding you through the process in class. As always, do ask questions and most importantly collaborate with your fellow students.

In order for you to be able to access Galaxy on your assigned dedicated machine on the Cloud, you have been given a web or IP address in the form of A.B.C.D where A, B, C and D are numbers separated by dots. You will need it in order to access Galaxy from the web browser on your laptop. 


#### Changes to the default Galaxy installation
{:.no_toc}

The default Galaxy instance relies on [TopHat/Cufflinks](http://tophat.cbcb.umd.edu/manual.html) for RNA-seq analysis, part of a well-tested workflow system. If you are not interested in discovering splice variants -- which can be prone to false positives -- a simple alignment with TopHat followed by analysis of differential expression with [DESeq](http://bioconductor.org/packages/release/bioc/html/DESeq.html) / [DSS](http://www.bioconductor.org/packages/release/bioc/html/DSS.html) (genes) or [DEXSeq](http://bioconductor.org/packages/release/bioc/html/DEXSeq.html) (exons) is likely to be simpler. For this workshop we added a more recent version ([TopHat2](http://genomebiology.com/2013/14/4/r36)) from the [Galaxy ToolShed](http://toolshed.g2.bx.psu.edu/) to our Galaxy instance, and will be running [edgeR](http://www.bioconductor.org/packages/release/bioc/html/edgeR.html) for the differential expression analysis.

> If you want to work through this session outside of this course, take a look at our section on [how to set up your own Galaxy server](http://bioinformatics.sph.harvard.edu/ngs-workshops/posts/running-your-own-galaxy-instance/) or download the [complete introduction dataset](https://dl.dropbox.com/u/407047/Blog/Data/RNA-Seq_Data_files.tgz) and use one of the public Galaxy instances.

## Quality control

The dataset we are using is part of a larger study described in [Kenny PJ et al, Cell Rep 2014](http://www.ncbi.nlm.nih.gov/pubmed/25464849). The authors are investigating interactions between various genes involved in Fragile X syndrome, a disease in which there is aberrant production of the FMRP protein. FMRP has been linked to the microRNA pathway, as it has been shown to be involved in miRNA mediated translational suppresion. **The authors sought to show that FMRP associates with the RNA helicase MOV10 for regulation of brain mRNAs.**

From this study we are using the [RNA-Seq](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE50499) data which is publicly available in the [SRA](http://www.ncbi.nlm.nih.gov/sra). The RNA was extracted from HEK293F cells that were transfected with a MOV10 transgene and normal control cells. Using this data, we will evaluate transcriptional patterns associated with MOV10 overexpression. 

The libraries for this dataset are stranded and were generated using the dUTP method. Sequencing was carried out on the Illumina HiSeq-2500 for 100bp single end reads. The full dataset was sequenced to ~40 million reads per sample, but for this workshop we will be looking at a small subset on chr1 (~300,000 reads/sample) to keep things manageable and allow algorithms to finish within a few minutes. We will switch over to using all samples for the differential expression analysis by retrieving them from the `Data Libraries` in the `Shared Data` menu as needed.For each group we have three replicates as described in the figure below.


![Experimental_design](../../images/screenshots/exp_design.png)


### Exploring the FASTQ files

Open up Galaxy from your machine as before, start up a new blank history and give it a name. In the top menu, move your mouse cursor over `Shared data` and select `Data libraries`. This allows you to import files into Galaxy that have been shared with you.

Browse to the `RNA-Seq` entry, open the`Sequence and reference data` folder and mouse over the arrow to the right of the filename `Mov10oe1.fastq`. In the context menu select `Import this dataset into selected histories`. Move back to the `Analyze Data` menu.

This is the kind of data you would expect to receive from your sequencing core facility: a large number of short nucleotide reads in [FASTQ format](https://dl.dropbox.com/u/407047/Blog/Documents/Literature/QC/Nucleic%20Acids%20Res%202009%20Cock.pdf) (PDF). 

> Take a look at the file contents using Galaxy's preview and file view functionality. What additional information does FASTQ contain over FASTA? How many reads are in your file? Pick the first sequence in your FASTQ file and note the identifier. What is the quality encoding character of it's first and last nucleotide, respectively?

## Quality Controls

The FASTQ file contains output reads from the sequencer that need to be mapped to a reference genome for us to understand where those reads came from on the sequenced genome. However, before we can delve into read mapping, we first need to make sure that our preliminary data is of sufficiently high quality. This involves several steps:

1. Obtaining summary quality statistics on the reads and review diagnostic graphs 
2. Eliminate sequencing artifacts
3. Filter out genetic contaminants (primers, vectors, adaptors)
4. Filter out low-quality reads
5. Recalculate quality statistics and review diagnostic plots on filtered data

Iterate through steps 2-5 until the data is of sufficient quality before proceeding to mapping.

### Obtain Quality Statistics

Under the `NGS:QC and manipulation` tool heading, select the `Compute quality statistics` tool (from the `FASTX toolkit` section) and apply it to your FASTQ file.

> Explore the results. What information does the data matrix contain?

Galaxy contains a number of tools to convert this matrix into graphs (e.g., `Draw Quality Score Boxplot` and `Draw Nucleotides Distribution Chart`), but in most cases you will stick to using a standard QC tool such as [FASTQC](http://www.bioinformatics.bbsrc.ac.uk/projects/fastqc/ ) with provides many more diagnostic plots. Find `FASTQC` in the Galaxy tool kit and give it a try.

> Do you observe similar quality statistics as before? What additional information do you get with this tool? What can you say about the per-base GC content? Does this match the human genome? Why is there a shift in the GC distribution from the theoretical distribution? Finally, what is the meaning of the k-mer content tab? This publication from [Schroeder et al](https://dl.dropbox.com/u/407047/Blog/Documents/Literature/QC/PLoS%20ONE%202010%20Schr%C3%B6der.pdf) (PDF) has some background information.

Compare your results with a [pre-generated FastQC plot](https://dl.dropboxusercontent.com/u/74036176/Mov10oe_1-fastqc_report.html) obtained from the full FASTQ data set (i.e., using all reads). 

> Do you note any differences in the output or format? Why is that?

Some of the diagnostic plots from a [Core conference call](http://bioinfo-core.org/index.php/9th_Discussion-28_October_2010) might be useful if you are curious to see what different modes of failure frequently look like.

Depending on your FASTQ results you may need to test for contamination (either by vector or other sourcers) and/or remove adapter sequences from your read prior to alignment. It is good to remove these sequences to increase your mapping efficiency. 

### Screen for contamination

Not all required tools are available through Galaxy just yet. For example, depending on the source of your data you probably want to screen your reads for vector or adapter contamination using tools such as [SeqTrim](https://dl.dropbox.com/u/407047/Blog/Documents/Literature/QC/BMC%20Bioinformatics%202010%20Falgueras.pdf) (PDF). Other tools check for sample cross-contamination ([ContEst](https://dl.dropbox.com/u/407047/Blog/Documents/Literature/QC/Bioinformatics%202011%20Cibulskis.pdf) (PDF)) or contamination with data from other species ([FASTQ Screen](http://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/)).  For now we will just explore how you can spot potential contamination from the FastQC results alone. 

> Where is the most overrepresented sequence in your FastQC results from the Mov10oe1 sample coming from?

Compare your own FastQC results with a [pre-generated FastQC summary](http://dl.dropbox.com/u/4253254/Training/HSCI_Workshop/D19VHACXX_s1_0_illumina12index_12_SL17357_fastqc/fastqc_report.html) that highlights a potential problem.

> Study the report. Is there any vector contamination present? Is there any adapter contamination present?

As an emergency solution to remove contaminating sequences, Galaxy comes with a way to `Clip adapter sequences`. If you notice a particularly strong contamination with a vector or primer you can use your FASTQ file as the `Library to clip`, keep the default `Minimum sequence length`, and enter the sequence you spotted in the FASTQC report as the `Adapter` to clip. You probably want to keep both clipped and unclipped sequences in this case; keeping only reads that contained an adapter or barcode can be useful if you expect this to be a feature of all valid reads. For command-line use tools such as [cutadapt](https://code.google.com/p/cutadapt/) or [AlienTrimmer](ftp://ftp.pasteur.fr/pub/gensoft/projects/AlienTrimmer/) are more versatile and less labor-intensive to use.
 
Normally, you would take the filtered output and feed it back into Galaxy, but our dataset has very little contamination so we will move on with the original data set. 

### Filter by quality

As we do not plan on using our data to identify genomic variation individual sequencing errors become less of a problem as long as the sequence can be still mapped to the reference genome. Still, you want to filter sequences where the overall quality is just too low. 

After reviewing the quality diagnostics from your FASTQC report, choose an appropriate lower nucleotide quality that you think would be good for this data (for instance, a quality score below 20 means there is more than a 1/100 chance of an incorrect base call). Use this cutoff to filter out low-quality reads from your data with the `Filter By Quality` tool. After you have filtered your data generate another `FASTQC` quality report.

> Do you have an improvement in read quality? What has changed and why? What percentage of reads did you retain at your cutoff?


### Read Trimming

After the general filtering which removes sequence reads of overall low quality it can be a good idea to trim the end of your reads, cutting away nucleotides at the end that are still of low-quality

> Is this necessary in your case? 

Use the `Trim sequences` tool to trim if required. Make sure to name your finalized FASTQ file that has been filtered and improved. You will be using this file in the next section for mapping. Again, tools such as [sickle](https://github.com/najoshi/sickle) provide more flexibility by allowing for adaptive, window-based trimming, but are not yet part of the default Galaxy toolkit -- though they are available through the ToolShed.

It might also be helpful to create a new history at this time and copy over your final sequence file. An easy way to do this is through the `History` menu and the `Copy datasets` function which allows you to pick one or more items in your current history, copying them over to a new or existing one. 



If for some reason any of the previous steps did not work you can retrieve a quality-controlled, trimmed and filtered FASTQ file from the `Data Library` in the `Snapshots` folder (`mov10oe1_filtered.fastq`).

## Read alignment

Next we are going to look at the steps we need to take once we have a clean, filtered FASTQ file that is ready for alignment. Use the filtered FASTQ file that you prepared in the previous section. The alignment process consists of the following steps:

* Choose an appropriate reference genome to map your reads against
* Choose an appropriate gene annotation model to guide your alignment
* Perform the read alignment

### Load the filtered FASTQ file and Reference Genome

To avoid excessive runtimes during later steps we will not align the reads against the whole human genome, just _chr1_. We have retrieved the sequence for this chromosome from UCSC GoldenPath (hg19) and added it to the `Data Library`. This means we won't use the pre-built genome indices, but submit our own reference sequence. We are also using only a subset of reads from a single sample for the next few steps. 

1. Start up a new blank history and give it a name.
2. Either load your filtered FASTQ files from the previous tutorial, or download a filtered FASTQ file from the Shared Data Libraries in Galaxy. 
3. Find the chromosome 1 genomic sequence in `RNA-Seq`, `Sequence and Reference Data` as `chr1.fa` and import it into your new history. 

Using it generates a slight time overhead since the reference genome needs to be indexed prior to each alignment run, but this is acceptable for now. 

### Read Alignment: TopHat

Aligning RNA-seq reads to a reference genome (`hg19`) introduces unique issues. While alignment and sequencing errors are not a big problem -- you mostly care about reads being mapped to the correct gene -- sequence data was generated from _spliced_ RNA. In other words, while the majority of reads will map to exons a significant share will map to _splice junction boundaries_ not present as such in the reference genome, and few standard NGS aligners can handle alignment gaps that can easily be 10kb or more. This is where dedicated mappers such as TopHat come in.

#### How does TopHat work?

Straight from the [TopHat manual](http://tophat.cbcb.umd.edu/manual.html):

_"TopHat finds splice junctions without a reference annotation. By first mapping RNA-Seq reads to the genome, TopHat identifies potential exons, since many RNA-Seq reads will contiguously align to the genome. Using this initial mapping information, TopHat builds a database of possible splice junctions, and then maps the reads against these junctions to confirm them.
Short read sequencing machines can currently produce reads 100bp or longer, but many exons are shorter than this, and so would be missed in the initial mapping. TopHat solves this problem by splitting all input reads into smaller segments, and then mapping them independently. The segment alignments are "glued" back together in a final step of the program to produce the end-to-end read alignments._

_TopHat generates its database of possible splice junctions from three sources of evidence. The first source is pairings of "coverage islands", which are distinct regions of piled up reads in the initial mapping. Neighboring islands are often spliced together in the transcriptome, so TopHat looks for ways to join these with an intron. The second source is only used when TopHat is run with paired end reads. When reads in a pair come from different exons of a transcript, they will generally be mapped far apart in the genome coordinate space. When this happens, TopHat tries to "close" the gap between them by looking for subsequences of the genomic interval between mates with a total length about equal to the expected distance between mates. The "introns" in this subsequence are added to the database. The third, and strongest, source of evidence for a splice junction is when two segments from the same read are mapped far apart, or when an internal segment fails to map."_

Although TopHat can find splice junctions without a reference annotation, we will give it some guidance by supplying it with one. We will use the UCSC chr19 gene annotations as a gene annotation model (`chr1-hg19_genes.gtf` from the `Sequence and reference data` folder in the `RNA-Seq` library). 

* Go ahead and import `chr1-hg19_genes.gtf` to your history

We will use the `Tophat2` tool from the `NGS:RNA Analysis` folder to map your trimmed reads. Remember to use `chr1.fa` from your history instead of the built-in index. To use the reference gene annotation model, do the following: 

1. Under `TopHat settings to use:` select `Full Parameter List` 
2. Scroll down to `Use Own Junctions` and select `Yes`
3. Scroll down again to `Use Gene Annotation Model` and select `Yes`
4. Scroll down to `Gene Model Annotations` and confirm that `chr1-hg19_genes.gtf` is selected

Finally, to run the alignment, click `Execute`. TopHat will generate 4 files, describing:   

* accepted hits - a BAM file containing the read alignments which TopHat thinks map to a gene  
* predicted splice junctions in BED format
* predicted deletions in BED format
* predicted insertions in BED format

It can take a while for TopHat to finish the alignment, even for such a small data set. We _strongly_ recommend that you take the time to work through the [TopHat documentation](http://tophat.cbcb.umd.edu/manual.html) to understand the data sets that are being generated by the alignment process, what different parameters are available, and how you can modify the output.


### Convert BAM files to SAM files

Most aligners produce output in ["BAM" format](http://samtools.sourceforge.net/)  by now (the binary version of the Sequence Alignment/Map (SAM) format, see [publication from Heng Li](https://dl.dropbox.com/u/407047/Blog/Documents/Literature/Exome%20Seq/Bioinformatics%202009%20Li-3.pdf) (PDF) for more details). In order to take a look at the data structure we need to convert the `accepted hits` BAM file into the text-based SAM format.

* Apply the `BAM-to-SAM` tool under the `NGS:SAM Tools` heading to convert your output BAM file to a readable SAM format. Make sure to include the header.

> Explore the SAM format. What key information is contained?

### Evaluating TopHat results

Explore the TopHat output from the previous session. We can get some basic information from the BAM file by calculating the index statistics using the `BAM Index Statistics` tool under the `NGS:Picard` tool heading (run it on the `Accepted Hits` output).

> What kinds of statistics are reported? What do they tell you about the alignment?

Also take a look at the other reports.

> How many splice junctions were generated? What is the distribution of reads covering splice junctions? The score column in the `Splice Junctions` BED formatted report can help with this (the score here is the number of alignments spanning the junction).


### Explore alignment in a genome browser

Take a look at your alignment in IGV (Integrative Genomics Viewer). To do so, first expand the alignment result (the `accepted hits` BAM file) in your history and click on the `web current` link beside `display with IGV`. This will download IGV to your machine. It should start automatically after a few security warnings (if not, find the `igv.jlnp` file and click it yourself).

You should see a window that looks like this:

[![IGV Startup](../../images/screenshots/IGV_RNA_1.png)](../../images/screenshots/IGV_RNA_1.png)

Click in the search bar and enter the gene name MOV10:
 
[![Highlight](../../images/screenshots/IGV_RNA_2.png)](../../images/screenshots/IGV_RNA_2.png)

You should now be able to see the individual reads, you can zoom in further using the controls on the upper right and scroll around by clicking and dragging in the highlighted alignment track. Note that IGV automatically translated the gene symbol into the matching genomic coordinates for your genome build:

[![Sample Gene](../../images/screenshots/IGV_RNA_3.png)](../../images/screenshots/IGV_RNA_3.png)

> Can you find any reads that span an intron (i.e. the read is derived from processed RNA)?

Next, add the TopHat generated `splice junctions` (BED format). As Galaxy does not support automatic export of BED files to IGV, you will have to add it to IGV yourself. To do so, _download_ the file from Galaxy and take note of where you saved it. Now move over to IGV and click the `File` menu and then `Load from file`. Select the file you downloaded and open it in IGV. 

> You can do the same to import the `chr1-hg19_genes.gtf` UCSC gene model into IGV, though the IGV RefSeq track should suffice to get a sense of the known gene models. 

If you have moved around in the browser navigate back to MOV10. You should now see the the putative splice junctions (highlighted track)…

[![Junctions](../../images/screenshots/IGV_RNA_4.png)](../../images/screenshots/IGV_RNA_4.png)

... in addition to the mapped reads and gene model (which can serve as an internal control). 

> Have a look at the genes PPM1J and PTPN22 as well. 
> For a different view of the predicted splice junctions, right click on the junctions track and select `Collapsed`.
> Can you find cases where the gene model predicts a splice junction but TopHat could not predict any? 
> Likewise, can you find regions where TopHat suggests a splice junction not supported by the RefSeq annotation?

As before, if any of the previous steps did not work out you can retrieve intermediate results from the `Snapshot` library under `Shared Data`. We have deposited the TopHat2-aligned reads in BAM format (`mov10oe1_accepted_hits.bam`).

If you are just interested in differential gene expression a very simple option is to feed the TopHat2-aligned reads directly into `CuffDiff`, the final tool of this session. The [UC Davis RNA-Seq course](http://training.bioinformatics.ucdavis.edu/docs/2012/05/RNA/index.html) has basic information on how to do so, although you most like would want to switch to other tools such as DESeq or edgeR (see below). 

Here, we first explore differential expression after transcript assembly which has the drawback of generating more false positives, but can also detect novel transcripts before moving on to a comparison of read counts with edgeR.


## Assemble transcripts

A number of different approaches exist to quantify gene expression counts and assess significant differences between sample groups (check the resources section for pointers). Most require some sort of unification between samples, for example either generating or providing gene models between which to compare the mapped read counts.

We will start by assembling the reads mapped to exons and splice junctions into complete transcripts using [Cufflinks](http://cufflinks.cbcb.umd.edu/manual.html#cufflinks). Similar to TopHat the algorithm comes with a large number of parameters; we recommend working through the manual before applying it to your own data. Like TopHat, CuffLinks can also run naively (i.e., assembling transcripts _de novo_ rather than using known gene models), but it is less prone to false positives when guided by a reference annotation. We will supply it with the same annotation we used for the TopHat alignment:

1. Select `Cufflinks` from the "NGS:RNA Analysis" section
2. Select your "accepted hits" output from TopHat2 as the `SAM or BAM file of aligned RNA-Seq reads:`
3. Under `Use Reference Annotation` select `Use reference annotation as guide` and  make sure `chr1-hg19_genes.gtf` is selected 
4. Under `Perform Bias Correction`, select `Yes`
5. Under `Reference sequence data` select `History` and confirm that the `chr1.fa` from your history is selected under `Using reference file`
6. Under `Use multi-read correct`, select `Yes`

While CuffLinks is running to browse through its manual or explore a paper from the Pachter group on [biases in RNA-Seq expression estimates](http://genomebiology.com/content/pdf/gb-2011-12-3-r22.pdf) (PDF) which explains some of the corrections applied by CuffLinks. The algorithm will generate three files describing the assembled transcripts, the transcript expression and the gene expression.

* Open the "assembled transcripts" output and explore its contents to find transcripts with multiple exons
* Add the assembled transcripts to the IGV browser view and compare it to the gene model. As before, you cannot send the file directly to IGV, so you will have to go the `Download` and `Load From File` route again.

> Try to assess how well Cufflinks did in generating 'complete' transcript models. `Expand` the track to do so if needed.


## Merge transcripts from different samples

In your own analyses, you will have multiple conditions and (hopefully!) multiple replicates for those conditions. Each sample will have its own set of assembled  transcripts or gene model, which can make comparing the samples a challenge. To do so, you will need to 'merge' (unify) the gene models between different samples using [Cuffmerge](http://cufflinks.cbcb.umd.edu/manual.html#cuffmerge).

Cuffmerge helps analyze the transcribed fragments in an assembly by comparing assembled transcripts to a reference annotation using [Cuffcompare](http://cufflinks.cbcb.umd.edu/manual.html#cuffcompare), filtering out likely artifacts and merging the finally merging the multiple sample derived assemblies.

### Bringing in replicates and another sample
 
Up to now we have been working with a single sample to save time. Now we will bring the other three ENCODE data sets into our analysis. We have generated the necessary files -- the second H1hESC replicate and two replicates for the CD20 cell line -- for you and placed them in the Shared Data Library `RNA-Seq`; they can be found in the `Additional sample data` folder. Open up this folder and import the following assembled transcripts for each of the other 5 ENCODE samples to your history:

* `Ctrl1_transcripts.gtf`
* `Ctrl2_transcripts.gtf`
* `Ctrl3_transcripts.gtf`
* `Mov10oe2_transcripts.gtf`
* `Mov10oe3_transcripts.gtf`

Run `Cuffmerge` from the `RNA Analysis` tool menu on these "assembled transcripts" gene models along with the "assembled transcript" gene model generated for the `Mov10oe1` sample by Cufflinks in the previous step. You will need a total of six input files, use the `Add new Additional GTF Input Files` button to add them all.

Again, use the `chr1-hg19_genes.gtf` gene model as a reference annotation. Use sequence data `chr1.fa` from your history as the `source for the reference list`; this will add additional annotation to your output.

> Time permitting, explore the resulting GTF file of merged transcripts in IGV. What transcripts got removed or merged?


## Differential expression with CuffDiff

As the last step, [Cuffdiff](http://cufflinks.cbcb.umd.edu/manual.html#cuffdiff) will use the TopHat2-aligned reads (the `Accepted Hits` file) together with the unified transcript model (the 'merged transcripts' file generated by Cuffmerge) to find transcripts exhibiting differential expression between your samples. 

> Note that doing so without biological replicates is not advised as the results would be _very_ unreliable.

For this step, you will need the `accepted hits` Tophat2 alignment output for all four samples. Again, we have prepared these for you and placed them in the Shared Data Library `RNA-Seq` in the `Additional sample data` folder. Open up the folder with the other additional samples data and import the following "bam" files for each of the other 5 samples to your history. 

* `Ctrl1_accepted_hits.bam`
* `Ctrl2_accepted_hits.bam`
* `Ctrl3_accepted_hits.bam`
* `Mov10oe2_accepted_hits.bam`
* `Mov10oe3_accepted_hits.bam`

You now have all data required to identify differentially expressed transcripts and genes. Use the imported "accepted hits" BAM files together with the data you generated for your own `Mov10oe1` replicate and feed them into `Cuffdiff`:

1. Under `Transcripts` select your Cuffmerge generated unified gene model. 
2. Using `Add new Condition` make sure you have entry fields for two conditions, one for the Mov10 samples and the other for the Ctrl samples.
3. Using `Add new Replicate` add three replicates for each group, and under the `Add replicate` headings select the appropriate "accepted hits" BAM files for each of your groups and replicates. The file names can be difficult to differentiate so make sure you add the right samples to each group (ie., use the history step numbers instead of the filenames to identify them).
4. Under `Library normalization method` pick `geometric`
6. Under `Use multi-read correct`, select `Yes` 
7. Under `Perform Bias Correction`, select `Yes`
8. Under `Reference sequence data` select `History` and confirm that the `chr1.fa` from your history is selected under `Using reference file`

Finally, click `Execute` and use the time while Cuffdiff is running to go over the output file documentation. Once the comparison has finished find the "transcript FPKM tracking" entry in your history and explore its contents.

 > Find an entry for a _novel_ isoform and an entry for a transcript matching the annotation you provided (hint: look at the 'Class Code' section of the Cufflinks manual). List the FPKM value and overlapping (or nearest) gene.
 
Next, look at the "gene differential expression testing" entry in your history. Sort the results by `q-value` and filter to select transcripts that are differentially expressed (use Galaxy here, but now might also be a good time to export the spreadsheet to any other tool you prefer.  

> * What are the p- and q-values for the top differentially expressed transcript? Load the "accepted hits" bam files into IGV and look at the gene; also explore the novel isoform transcript you looked up before.
> * Time permitting, send the two additional Mov10 replicates and the three Ctrl BAM files to the genomic viewer. Explore one of the genes found to be differentially expressed and see if the aligned reads support the CuffDiff call.

## Differential expression with edgeR

If you are not interested in isoform recovery but just want to know which genes are differentially expressed frameworks such as [edgeR](http://www.bioconductor.org/packages/release/bioc/html/edgeR.html) tend to have a lower false positive rate and allow for describing the experimental setup with generalized linear models. edgeR will _not_ work for experiments without replicates. We strongly recommend working through the user guide before applying edgeR to any of your experimental data sets.

### Summarizing count data

Tools like edgeR, DESeq and others need a simple count matrix that summarizes the number of reads mapping to features of interest; [htseq-count](http://www-huber.embl.de/users/anders/HTSeq/doc/count.html) is one such tool. 

* Under `NGS:RNA Analysis` find `Count reads in features with htseq-count`. 
* Similar to the TopHat/Cufflinks suite you need to again define your sample groups and replicates. Create two groups with two replicates each and enter the `accepted hits` files
* htseq-count needs to know what the features are that you want read summaries for. Set the `chr1-hg19_genes.gtf` as input to the `feature file` selector
* Set the `Mode to handle reads overlapping more than one feature` to `Union` and `Is the data from a strand-specific assay?` to `Yes but reverse`

Leave the other values at their default, but do take the time to read up on how edgeR handles reads overlapping more than one feature while the matrix is being generated. 

> Take a look at the replicate variance.

### Differential gene expression

The resulting htseq-count matrix is quite small since the default parameters ensures  features (genes) without mapped reads in any sample are dropped. This means edgeR will finish quickly as it only needs to run a couple of comparisons. 

* Find `edgeR` in the `RNA Analysis` tool panel
* Provide it with the matrix generated by `htseq-count`, leaving all other values at default and run the algorithm

> * Compare the `norm expr matrix` with the matrix generated by htseq-count. You should be able to see the read count normalization
> * From the `analysis plots` section take a look at the `MDS plot` to check how good the replicates are, and at the `tagwise dispersion plot` to see how genes with lower read counts tend to have higher variation

Explore some of the genes reported to be differentially expressed in IGV.

## Data visualization

Some frameworks allow for a more interaction exploration of RNA-Seq-based count data, but frequently give you less control over the actual comparison being run. [Degust](http://www.vicbioinformatics.com/degust/index.html) is one such tool that is worth exploring:

* Degust expects input files with only one header line. In Galaxy, find `Remove beginning` in the `Text manipulation` set. Use the non-normalized count matrix as input and remove the first (header) line, downloading the result to your desktop
* Navigate to [Degust](http://www.vicbioinformatics.com/degust/index.html) and start by clicking on `Upload your counts file` and you should be re-directed to another page. Here, you will `Choose your file` by navigating to the directory where you downloaded your counts matrix and then click `Upload`
* At the configuration page, give your project a name and select `TAB separated` as your Format type. Check all columns in the `Info columns` selector.

We need to specify to Degust which samples belong to which groups. Click on `Add condition`, and give the condition a name (i.e., Ctrl) and from the drop-down menu select the appropriate samples. Do the same for Mov10 Add the `Gene Link` column if you want to, then `Save changes` will save the configuration and press `View` to generate results. You should be directed to a window that looks similar to the screenshot below:

[![Degust](../../images/screenshots/mov10_degust.png)](../../images/screenshots/mov10_degust.png)

* At the very top you will find an MA plot showing expression for two conditions. Each dot represents a gene and is color-coded by FDR (red = more significant). You can click and drag on the plot to select genes within a rectangle - the heatmap and table below will be filtered to only those genes. 
* The heatmap shows the log fold-change for each gene shown in the MA plot above. Each vertical strip corresponds to a gene. Hovering the mouse over the heatmap will show the corresponding gene in the plot above

> In the top left hand corner, you will find a box where you can apply various filters on the FDR and logFC which will dynamically filter the plots and the gene table below. Experiment with these. 

The gene table should update automatically and show you the genes filtered to show the same genes as plotted above. Anything displayed in the gene list table can be downloaded as a CSV file, and double-clicking on a gene name should bring you to the matching NCBI page. Specific genes may be found using the search box


## Functional enrichment

A large number of methods and external tools are available to make sense of a gene list, most of which revolve around deteching functional enrichment in Gene Ontology, curated pathways or networks. A post on [Getting Genetics Done](http://gettinggeneticsdone.blogspot.com/2012/03/pathway-analysis-for-high-throughput.html) should get you started, and we will expand this section for future workshops. For now we will explore two standard approaches that cover the basics, giving you an idea whether your experimental setup worked.

Since there's relatively little we can get in terms of functional enrichment given that you are only looking at ~300,000 bp mapping to chr1, we have provided the [gene list](hhttps://dl.dropboxusercontent.com/u/204381225/Galaxy_nanocourses/Mov10_full_gene_list.txt) for the entire dataset.

### Gene Ontology

[g:Profiler](http://biit.cs.ut.ee/gprofiler_beta/) (beta) is a great framework to obtain a quick overview of gene lists in the 20-500 genes range. It can be used with default parameters, but it frequently makes sense to switch off the `Significant only` option and setting a custom pvalue threshold (the `User p-value` under `Advanced options`). It supports ranked and unranked gene lists and can also handle multiple gene lists at the same time, although [ToppCluster](http://toppcluster.cchmc.org/) tends to work better for these situations. 

* Cut and paste the gene symbols into the g:Profiler `Query` box; take your time explore some of the additional options
* Click on `g:Profile!`

[![g:Profiler](../../images/screenshots/mov10_gprofiler.png)](../../images/screenshots/mov10_gprofiler.png)

> Explore the results. The output is a spreadsheet with enrichment of GO terms, pathway / disease enrichments, miRNA and TF-targets as well as a basic visualization. It is also worthwhile to capture the gene names and descriptions table.

To provide additional graphic representations you can fall back to [Revigo](http://revigo.irb.hr/) which takes the previously generated GO terms and pvalues. The output is a quick scatterplot or treemap (see below) organized by semantic similarity, and more importantly a new GO list that highlights redundant terms which can be filtered out in the data summary provided back to researchers.

> Time permitting, extract the GO identifiers obtained from g:Profile and import them into Revigo. You probably want to switch the g:Profiler output from `Graphical (PNG)` to `Excel spreadhseet (XLS)` to be able to select GO identifiers and -- optionally -- p-values. 

[![Revigo](../../images/screenshots/mov10_revigo.png)](../../images/screenshots/mov10_revigo.png)

### Network-based

It is a bit more difficult to provide a workflow for any network-related tasks as these depend heavily on the biological question being asked. Network-based approaches can be useful in a number of cases including, but not limited to:

* Providing functional information for comparatively small gene lists (10-30 genes) by pulling in related genes, identifying interaction partners, etc.
* Comparing two or more gene lists beyond the basic gene list overlap, for example by exploring interactions between two gene lists, or expanding one gene list and again checking for interactions
* When trying to identify additional candidate genes, e.g. in cases where the primary analysis only revealed a small number of potential markers, or researchers are interested in identifiying putative new targets based on a small screen or literature review
* For using all the available data (i.e., a complete list of ranked genes or proteins) to identify subnetworks of interest

There are almost too many network-based tools out there, but [GeneMania](http://genemania.org/) has become one of the standard systems to annotate small to mid-sized gene lists, with or without adding functionally related proteins. The reporting function is extremely useful and exhaustive.

Navigate to the GeneMania website and cut and paste your identifiers again, making sure you pick `H. sapiens` as species. Open the advanced options and decide what kind of interaction information you would like to use to connect the genes in your list:

[![GeneMania](../../images/screenshots/mov10_genemania.png)](../../images/screenshots/mov10_genemania.png)

GeneMania enables you to add 'related' or linker genes to the results. Those are genes that were not present in your input list, but connect your genes in some way:

[![Linker](../../images/screenshots/mov10_genemania_linker.png)](../../images/screenshots/mov10_genemania_linker.png)

Pick a number and click on `Go`. 

> Explore the reporting functionality as well as how the enrichment terms change depending on whether you added linker genes or not, or in relation to the interaction information used.

That's it! You went from a number of FASTQ files (your primary data) all the way to a list of enriched biological functions in the matter of a few hours. Next, we recommend that you take the time to explore the reading and resources section and start working with your own data. Please do [contact us](mailto:bioinformatics@hsph.harvard.edu) with questions, comments and feedback.

